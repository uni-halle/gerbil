/*********************************************************************************
Copyright (c) 2016 Marius Erbert, Steffen Rechner

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies
of the Software, and to permit persons to whom the Software is furnished to do
so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS
FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR
COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER
IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*********************************************************************************/

#ifndef KMCHASHTABLE_H_
#define KMCHASHTABLE_H_

#include "FailureBuffer.h"
#include "SyncQueue.h"
#include "TempFile.h"
#include "KmerDistributer.h"
#include <algorithm>
#include <chrono>

namespace gerbil {
namespace cpu {

#define HT_HISTO_SIZE 512

#ifdef HT_JUMPS
#define IF_HT_JUMPS(x) x
#else
#define IF_HT_JUMPS(x)
#endif

#ifdef HT_FAILS
#define IF_HT_FAILS(x) x
#else
#define IF_HT_FAILS(x)
#endif

#ifdef HT_HISTO
#define IF_HT_HISTO(x) x
#else
#define IF_HT_HISTO(x)
#endif

/**
 * A HashTable for Counting of Kmers.
 */
template<unsigned K>
class HasherTask {

	//stat
	IF_HT_FAILS(atomic<uint64> _fKMersNumber;)
	IF_HT_JUMPS(atomic<uint64> _jumps;)

	uint64 _maxPartSize;		// number of entries in hash table
	uint64 _maxSizeUsage;		// ?

	SyncSwapQueueMPSC<KmcBundle>* _kmcSyncSwapQueue;	// output queue

	// distributes kmers to various hash tables
	KmerDistributer* distributor;

	TempFile* _tempFiles;		// array of temporary files?
	std::string _tempPath;		// path to temporary files

	uint_cv _thresholdMin;		// minimal number of occurrences to be output

	byte _threadsNumber;		// number of threads
	std::thread** _threads;		// array of threads

	// for statistic
	std::atomic<uint64> _kMersNumber;
	std::atomic<uint64> _uKMersNumber;
	std::atomic<uint64> _btUKMersNumber;

	IF_HT_HISTO(
			atomic<uint64> _histo[HT_HISTO_SIZE];
	)

	uint32 getMaxSteps() const;

public:
	HasherTask(const byte &threadsNumber, KmerDistributer* distributor,
			SyncSwapQueueMPSC<KmcBundle>* kmcSyncSwapQueues,
			TempFile* tempFiles, const uint_cv &thresholdMin,
			const uint64 &maxSize, std::string tempPath);

	~HasherTask();

	inline uint64 getKMersNumber() const {
		return _kMersNumber.load();
	}
	inline uint64 getUKMersNumber() const {
		return _uKMersNumber.load();
	}
	inline uint64 getBtUKMersNumber() const {
		return _btUKMersNumber.load();
	}

	void print();
	void printStat();

	void hash(SyncSwapQueueMPSC<cpu::KMerBundle<K>>** kMerQueue);
	void join();
};

template<unsigned K>
HasherTask<K>::HasherTask(const byte &threadsNumber,
		KmerDistributer* distributor,
		SyncSwapQueueMPSC<KmcBundle>* kmcSyncSwapQueue, TempFile* tempFiles,
		const uint_cv &thresholdMin, const uint64 &maxSize,
		std::string tempPath
) :
		_kmcSyncSwapQueue(kmcSyncSwapQueue), _tempFiles(tempFiles),
		_thresholdMin(thresholdMin), _threads(NULL), _tempPath(tempPath), _maxSizeUsage(0),
		_kMersNumber(0), _uKMersNumber(0), _btUKMersNumber(0), _threadsNumber(threadsNumber),
		distributor(distributor)
{
	_maxPartSize = maxSize / _threadsNumber;

	//stat
	IF_HT_FAILS(_fKMersNumber.store(0);)
	IF_HT_JUMPS(_jumps.store(0);)
	IF_HT_HISTO(
			for(uint i(0); i < HT_HISTO_SIZE; ++i)
			_histo[i].store(0);
	)
}

template<unsigned K>
HasherTask<K>::~HasherTask() {
}

template<unsigned K>
inline void HasherTask<K>::printStat() {
puts("");
IF_HT_FAILS(printf("fkmers         : %12lu\n", _fKMersNumber.load());)
IF_HT_FAILS(printf("fkmers/kmer    : %15.2f\n", (double)_fKMersNumber.load() / _kMersNumber);)
IF_HT_JUMPS(printf("jumps          : %12lu\n", _jumps.load());)
IF_HT_JUMPS(printf("jumps/kmer     : %15.2f\n", (double)_jumps.load() / _kMersNumber);)

double r = (double) (sizeof(KMer<K> ) + sizeof(uint_cv)) / 1024 / 1024;
uint64 maxSize = _maxPartSize * _threadsNumber;
printf("size unused    : %12lu of %lu\n", maxSize - _maxSizeUsage, maxSize);
printf("memory unused  : %.0f MB of %.0f MB\n", (maxSize - _maxSizeUsage) * r,
		maxSize * r);
}

#define KMCHT_fill() {																\
	if(HYBRID_COUNTER && useSort) {													\
		kmb->copyAndInc(nextKey);													\
		kmb->clear();																\
	} else {																		\
	while(kmb->next(nkMer)) {														\
		hPos = nkMer->getHash() % partSize;											\
		i = 0;																		\
		KMer<K>* curKey;															\
		while(true) {																\
			if((curKey = keys + hPos)->isEmpty()) {									\
				values[hPos] = 1;													\
				curKey->set(*nkMer);												\
				++binUKMers;														\
				break;																\
			}																		\
			if(curKey->isEqual(*nkMer)) {											\
				++values[hPos];														\
				break;																\
			}																		\
			if(++i > maxSteps) {													\
				outBuffer->addKMer(*nkMer);											\
				++binFKMers;														\
				break;																\
			}																		\
			hPos += i * i;															\
			hPos %= partSize;														\
		}																			\
		IF_HT_JUMPS(jumps += i;)													\
	}																				\
	kmb->clear();																	\
}}

#define KMCHT_extract() {																	\
	if(HYBRID_COUNTER && useSort) {															\
		std::sort(keys, nextKey); 															\
		KMer<K>* endKey_p = nextKey;														\
		KMer<K> key;																		\
		key.clear();																		\
		uint64 val = 0;																		\
		for(KMer<K>* curKey = keys; curKey < endKey_p; ++curKey) {							\
			if(curKey->isEqual(key)) {														\
				++val;																		\
			} else {																		\
				if(val) {/* key not empty*/													\
					++binUKMers;															\
					binKMers += val;														\
					IF_HT_HISTO(++histo[val < HT_HISTO_SIZE ? val : 0];)					\
					if(_thresholdMin > val)													\
						++btUKMersNumber;													\
					else if(!curKmcBundle->add<K>(key, val)) {								\
							IF_MESS_SUPERREADER(sw.hold();)									\
							_kmcSyncSwapQueue->swapPush(curKmcBundle);						\
							IF_MESS_SUPERREADER(sw.proceed();)								\
							curKmcBundle->add<K>(key, val);									\
					}																		\
				}																	\
				key.set(*curKey);													\
				val = 1;															\
			}																		\
			curKey->clear();														\
		}																			\
		if(val) {/* key not empty*/\
			++binUKMers;															\
			binKMers += val;														\
			IF_HT_HISTO(++histo[val < HT_HISTO_SIZE ? val : 0];)					\
			if(_thresholdMin > val)													\
				++btUKMersNumber;													\
			else if(!curKmcBundle->add<K>(key, val)) {								\
					IF_MESS_SUPERREADER(sw.hold();)									\
					_kmcSyncSwapQueue->swapPush(curKmcBundle);						\
					IF_MESS_SUPERREADER(sw.proceed();)								\
					curKmcBundle->add<K>(key, val);									\
			}																		\
		}																			\
	} else {																			\
		uint_cv val;																	\
		KMer<K>* endKey = keys + partSize;												\
		KMer<K>* curKey = keys;															\
		for(; curKey < endKey; ++curKey) {												\
			if(!curKey->isEmpty()) {													\
				binKMers += (val = *(values + (curKey - keys)));						\
				IF_HT_HISTO(++histo[val < HT_HISTO_SIZE ? val : 0];)					\
				if(_thresholdMin > val)													\
					++btUKMersNumber;													\
				else if(!curKmcBundle->add<K>(*curKey, val)) {							\
						IF_MESS_SUPERREADER(sw.hold();)									\
						_kmcSyncSwapQueue->swapPush(curKmcBundle);						\
						IF_MESS_SUPERREADER(sw.proceed();)								\
						curKmcBundle->add<K>(*curKey, val);								\
				}																		\
				curKey->clear();														\
			}																			\
		}																				\
	}																					\
}

// compute next size for hashtable (one thread)
#define KMCHT_setPartSize() {																	\
	useSort = _tempFiles[curTempFileId].getKMersNumber() <= _maxPartSize;						\
	if(HYBRID_COUNTER && useSort)																\
		nextKey = keys;																			\
	else {																						\
		partSize = apprUKMersNumber / FILL;														\
		if(partSize < 128)																		\
			partSize = 128;																		\
		if(partSize > _maxPartSize)																\
			partSize = _maxPartSize;															\
		if(maxPartSizeUsage < partSize)	{														\
			keysEnd = keys + partSize;															\
			for(KMer<K>* p = keys + maxPartSizeUsage; p < keysEnd; ++p)							\
				p->clear();																		\
			maxPartSizeUsage = partSize;														\
		}																						\
	}																							\
}

template<unsigned K>
void HasherTask<K>::hash(SyncSwapQueueMPSC<cpu::KMerBundle<K>>** kMerQueues) {
_threads = new std::thread*[_threadsNumber];
for (uint8_t i(0); i < _threadsNumber; ++i) {
	_threads[i] =
			new std::thread(
					[this](const byte tId, SyncSwapQueueMPSC<KMerBundle<K>>** kMerQueues) {
						IF_MESS_HASHER(
								StopWatch sw(CLOCK_THREAD_CPUTIME_ID);
								sw.start();
						)

						KMerBundle<K>* kmb = new KMerBundle<K>();
						KMerBundle<K>* skmb = new KMerBundle<K>();
						KmcBundle* curKmcBundle = new KmcBundle();
						const uint64 maxSteps(getMaxSteps());

						FailureBuffer<K>* inBuffer = new FailureBuffer<K>(FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, _tempPath, tId, 0);
						FailureBuffer<K>* outBuffer = new FailureBuffer<K>(FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, _tempPath, tId, 1);

						SyncSwapQueueMPSC<KMerBundle<K>>* kMerQueue = kMerQueues[tId];
						KMer<K>* nkMer;
						KMer<K>* keys = new KMer<K>[_maxPartSize];
						uint_cv* values = new uint_cv[_maxPartSize];

						KMer<K>* keysEnd;

						uint64 i;
						uint64 hPos;

						uint64 binKMers = 0;
						uint64 binUKMers = 0;
						uint64 binFKMers = 0;

						uint64 kMersNumber = 0;
						uint64 uKMersNumber = 0;
						uint64 fKMersNumber = 0;
						uint64 btUKMersNumber = 0;

						uint64 maxPartSizeUsage = 0;

						IF_HT_JUMPS(uint64 jumps = 0;)
						IF_HT_HISTO(
								uint64 histo[HT_HISTO_SIZE];
								for(uint i(0); i < HT_HISTO_SIZE; ++i)
								histo[i] = 0;
						)

						bool notEmpty = kMerQueue->swapPop(kmb);
						uint_tfn curTempFileId = notEmpty ? kmb->getTempFileId() : 0;
						uint64 partSize;
						float ratio  = distributor->getSplitRatio(false, tId, curTempFileId);
						//printf("cpu hasher thread %i: my ratio is %f\n", tId, ratio);
						uint64 apprUKMersNumber = _tempFiles[curTempFileId].approximateUniqueKmers(FILL) * ratio;
						bool useSort;
						KMer<K>* nextKey;

						KMCHT_setPartSize();

						// for time measurement
						typedef std::chrono::microseconds ms;
						std::chrono::steady_clock::time_point start, stop;
						ms duration(0);

						while(notEmpty) {

							//printf("cpu hasher thread %i: not empty\n", tId);

							if(curTempFileId == kmb->getTempFileId()) {
							
								// insert kmers into hash table
								start = std::chrono::steady_clock::now();
								KMCHT_fill();
								stop = std::chrono::steady_clock::now();
								duration += std::chrono::duration_cast<ms>(stop-start);
								
								IF_MESS_SUPERREADER(sw.hold();)
								if(kMerQueue->swapPop(kmb)) {
									IF_MESS_SUPERREADER(sw.proceed();)
									continue;
								}
								else
								notEmpty = false;
								IF_MESS_SUPERREADER(sw.proceed();)
							}

							/* New File! Extract all kmers from Table */

							// measure time for extracting all kmers
							start = std::chrono::steady_clock::now();
							KMCHT_extract();

							// handle failures
							std::swap(skmb, kmb);
							bool firstFails = true;
							while(!outBuffer->isEmpty()) {
								std::swap(inBuffer, outBuffer);
								IF_DEB_DEV(
										if(firstFails)
										putchar('.');
										else
										putchar('!');
								)
								// fill failures
								if(firstFails && binFKMers > 1024) {
									apprUKMersNumber = binFKMers * ((double)binUKMers / binKMers) * 1.2; // +20% space
									firstFails = false;
								}
								else
								apprUKMersNumber = binFKMers;

								fKMersNumber += binFKMers;
								KMCHT_setPartSize();
								//printf("%1u, %4u: %9lu, %9lu, %9lu => %9lu ==> %9lu\n", tId, curTempFileId, binKMers, binUKMers, binFKMers, apprUKMersNumber, partSize);

								binFKMers = 0;

								while(inBuffer->getNextKMerBundle(kmb))
								KMCHT_fill();

								KMCHT_extract();

								inBuffer->clear();
							}
							std::swap(skmb, kmb);
							
							// stop timing
							stop = std::chrono::steady_clock::now();
							duration += std::chrono::duration_cast<ms>(stop-start);

							// report throughput to the kmer distributor
							float throughput = duration.count() == 0 ? 0 : binKMers / duration.count();
							this->distributor->updateThroughput(false, tId, throughput);

							kMersNumber += binKMers;
							uKMersNumber += binUKMers;
							fKMersNumber += binFKMers;

							_tempFiles[curTempFileId].incUKMersNumber(binUKMers);

							// preparation next bin
							if(notEmpty) {
								binKMers = 0;
								binUKMers = 0;
								binFKMers = 0;
								ratio  = distributor->getSplitRatio(false, tId, curTempFileId);
								//printf("cpu hasher thread %i: my ratio is %f\n", tId, ratio);
								curTempFileId = kmb->getTempFileId();
								apprUKMersNumber = _tempFiles[curTempFileId].approximateUniqueKmers(kMersNumber ? (double)uKMersNumber / kMersNumber : FILL) * ratio;
								KMCHT_setPartSize();
								// reset timer
								duration = ms(0);
							}
						}
						if(!curKmcBundle->isEmpty())
						_kmcSyncSwapQueue->swapPush(curKmcBundle);

						IF_HT_HISTO(
								for(uint i(0); i < HT_HISTO_SIZE; ++i)
								_histo[i] += histo[i];
						)

						delete kmb;
						delete skmb;
						delete curKmcBundle;
						delete[] keys;
						delete[] values;
						delete inBuffer;
						delete outBuffer;

						_kMersNumber += kMersNumber;
						_uKMersNumber += uKMersNumber;
						_btUKMersNumber += btUKMersNumber;

						//printf("Hasher[%2u]: kmers: %12lu    ukmers: %12lu\n", tId, kMersNumber, uKMersNumber);

						__sync_add_and_fetch(&_maxSizeUsage, maxPartSizeUsage);

						IF_HT_FAILS(_fKMersNumber += fKMersNumber;)
						IF_HT_JUMPS(_jumps += jumps;)
						IF_MESS_SUPERREADER(
								sw.stop();
								printf("Hasher[%2u]: %.3f s\n", tId, sw.get_s());
						)

					}, i, kMerQueues);
}
}

template<unsigned K>
inline void HasherTask<K>::join() {
for (uint i = 0; i < _threadsNumber; ++i) {
	_threads[i]->join();
	delete _threads[i];
}
delete[] _threads;
IF_HT_HISTO(
	printf("-------------------------------------\n");
	printf("count\tnumber\n");
	for(uint i(1); i < HT_HISTO_SIZE; ++i)
	printf("%u\t%lu\n", i, _histo[i].load());
	printf(">=%u\t%lu\n", HT_HISTO_SIZE, _histo[0].load());
	printf("-------------------------------------\n");
)
}

// private
template<unsigned K>
uint32 HasherTask<K>::getMaxSteps() const {
uint32 ms = 5;
uint64 mp = _maxPartSize;
while (mp >>= 1)
++ms;
return ms;
}

}
}

#endif /* KMCHASHTABLE_H_ */

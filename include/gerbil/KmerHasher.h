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

#ifndef KMERHASHER_H_
#define KMERHASHER_H_

#include "ThreadBarrier.h"
#include "Bundle.h"
#include "CpuHasher.h"
#include "GpuHasher.h"
#include "KmerDistributer.h"

namespace gerbil {

/*
 * This class does the counting of the kmers. It reads SuperBundles and stores
 * counting results as KmcBundles. The class starts a number of threads that
 * use hash tables for counting.
 */
class KmerHasher {

private:

	// input queue for super bundles
	SyncSwapQueueSPMC<SuperBundle>* _superBundleQueue;

	// output queue for kmer count bundles
	SyncSwapQueueMPSC<KmcBundle> _kmcSyncSwapQueue;

	// threading stuff
	std::thread** _processSplitterThreads;
	std::thread* _processThread;
	uint8 _processSplitterThreadsNumber;
	uint8 _numCPUHasher;
	uint8 _numGPUHasher;
	Barrier* barrier;								// a thread barrier

	// The kmer distributor observes throughput of the hash tables and
	// determines the ratio of kmers of a file for each hash table
	KmerDistributer *distributor;

	// auxilliary variables
	const uint32_t _k;
	const uint_tfn _tempFilesNumber;	// number of temporary files
	TempFile* _tempFiles;				// array of temporary files
	const uint32 _thresholdMin;			// minimal occurance of kmers to be output
	uint_tfn* _tempFilesOrder;			// processing order of temporary files
	std::string _tempFolder;			// directory of temporary files ?
	const bool _norm;					// whether to use normalized kmers

	// for synchronization
	uint32 _syncSplitterCounter;
	std::mutex _mtx;
	std::condition_variable _cv_barrier;

	// for statistic
	uint64 _kMersNumberCPU, _kMersNumberGPU;
	uint64 _uKMersNumberCPU, _uKMersNumberGPU;
	uint64 _btUKMersNumberCPU, _btUKMersNumberGPU;

	uint64 _maxKmcHashtableSize;
	uint32 _kMerBundlesNumber;

	uint64 _histogram[HISTOGRAM_SIZE];

	/**
	 * Spawns Splitter Threads.
	 */
	template<uint32_t K>
	void processSplit(SyncSwapQueueMPSC<cpu::KMerBundle<K>>** cpuKMerQueues,
			SyncSwapQueueMPSC<gpu::KMerBundle<K>>** gpuKMerQueues) {

		// spawn new threads that do the splitting
		for (uint8_t tId = 0; tId < _processSplitterThreadsNumber; ++tId) {
			_processSplitterThreads[tId] =

					new std::thread(
							[this](const uint8_t &tId, SyncSwapQueueMPSC<cpu::KMerBundle<K>>** cpuKMerQueues, SyncSwapQueueMPSC<gpu::KMerBundle<K>>** gpuKMerQueues) {

								if(_norm) {
									// use normalized kmers
									processThreadSplit<K, true>(tId, cpuKMerQueues, gpuKMerQueues);
								}
								else {
									// do not use normalized kmers
									processThreadSplit<K, false>(tId, cpuKMerQueues, gpuKMerQueues);
								}
							}, tId, cpuKMerQueues, gpuKMerQueues);
		}
	}

	template<uint32_t K>
	void process_template() {

		_processThread =
				new std::thread([this]() {
#ifdef DEB
						printf("KmerHasher start...\n");
						StopWatch swr;
						StopWatch swt(CLOCK_THREAD_CPUTIME_ID);
						swr.start();
						swt.start();
#endif

						// create Hash Tables (CPU-side)
						cpu::HasherTask<K> cpuHasher(_numCPUHasher,
								distributor,
								&_kmcSyncSwapQueue, _tempFiles, _thresholdMin,
								_maxKmcHashtableSize, _tempFolder);

						// create  Hash Tables (GPU-side)
						gpu::HasherTask<K> gpuHasher(_numGPUHasher,
								distributor,
								&_kmcSyncSwapQueue, _tempFiles, _tempFolder, _thresholdMin);

						// Create Queues for kmer bundles
						SyncSwapQueueMPSC<cpu::KMerBundle<K>>** cpuKMerQueues;
						SyncSwapQueueMPSC<gpu::KMerBundle<K>>** gpuKMerQueues;

						cpuKMerQueues = new SyncSwapQueueMPSC<cpu::KMerBundle<K>>*[_numCPUHasher];
						gpuKMerQueues = new SyncSwapQueueMPSC<gpu::KMerBundle<K>>*[_numGPUHasher];

						const uint32_t total_threads = _numCPUHasher + _numGPUHasher;

						for (uint8 i = 0; i < _numCPUHasher; ++i) {
							// TODO: Calibrate the argument here
							cpuKMerQueues[i] = new SyncSwapQueueMPSC<cpu::KMerBundle<K>>(
									_kMerBundlesNumber / total_threads);
						}

						for (uint8 i = 0; i < _numGPUHasher; ++i) {
							// TODO: Calibrate the argument here
							gpuKMerQueues[i] = new SyncSwapQueueMPSC<gpu::KMerBundle<K>>(
									_kMerBundlesNumber / total_threads);
						}

						// start splitting
						processSplit<K>(cpuKMerQueues, gpuKMerQueues);

						// start hashing
						cpuHasher.hash(cpuKMerQueues);
						gpuHasher.hash(gpuKMerQueues);

						// wait for splitters to finalize
						joinSplitterThreads();

						// wait for queues to finalize
						for(uint8_t i=0; i < _numCPUHasher; i++)
						cpuKMerQueues[i]->finalize();

						for(uint8_t i=0; i<_numGPUHasher; i++)
						gpuKMerQueues[i]->finalize();

						// wait for hashtables to finalize
						cpuHasher.join();
						gpuHasher.join();

						// wait for kmc queue to finalize
						_kmcSyncSwapQueue.finalize();

						// clean up
						for (uint8 i = 0; i < _numCPUHasher; ++i)
						delete cpuKMerQueues[i];

						for (uint8 i = 0; i < _numGPUHasher; ++i)
						delete gpuKMerQueues[i];

						delete[] cpuKMerQueues;
						delete[] gpuKMerQueues;

						// statistical / debug output
						cpuHasher.printStat();
						for(uint i = 0; i < HISTOGRAM_SIZE; ++i)
							_histogram[i] = cpuHasher.getHistogramEntry(i) + gpuHasher.getHistogramEntry(i);

						_kMersNumberCPU = cpuHasher.getKMersNumber();
						_uKMersNumberCPU = cpuHasher.getUKMersNumber();
						_btUKMersNumberCPU = cpuHasher.getBtUKMersNumber();
						_kMersNumberGPU = gpuHasher.getKMersNumber();
						_uKMersNumberGPU = gpuHasher.getUKMersNumber();
						_btUKMersNumberGPU = gpuHasher.getBtUKMersNumber();
					});

	}

	/**
	 * * This macro distributes the kmer to one of the kmer bundle queues.
	 */
#define ADD_KMER_TO_BUNDLE(kmer, curTempFileId)													\
	{																							\
		/* 																						\
		 * determine id of hash table according to current										\
		 * split ratio.																			\
		 */																						\
		const uint32_t h = distributor->distributeKMer<K>(kmer, curTempFileId);					\
																								\
		/* assign kmer to bundle cpu or gpu */													\
		if(h >= _numCPUHasher) {																\
			/* add to gpu */																	\
			if (!gpuKMerBundles[h-_numCPUHasher]->add(kmer)) {									\
				gpuKMerBundles[h-_numCPUHasher]->setTempFileId(curTempFileId);					\
				gpuKMerQueues[h-_numCPUHasher]->swapPush(gpuKMerBundles[h-_numCPUHasher]);		\
				gpuKMerBundles[h-_numCPUHasher]->add(kmer);										\
			}																					\
		}																						\
		else {																					\
			/* add to cpu */																	\
			if (!cpuKMerBundles[h]->add(kmer)) {												\
				cpuKMerBundles[h]->setTempFileId(curTempFileId);								\
				cpuKMerQueues[h]->swapPush(cpuKMerBundles[h]);									\
				cpuKMerBundles[h]->add(kmer);													\
			}																					\
		}																						\
	}

	/**
	 * The Task of each Splitter Thread:
	 * Split SuperMers into kmers and store them into parameter queues.
	 */
	template<uint32_t K, bool NORM>
	void processThreadSplit(const uint8_t &threadId,
			SyncSwapQueueMPSC<cpu::KMerBundle<K>>** cpuKMerQueues,
			SyncSwapQueueMPSC<gpu::KMerBundle<K>>** gpuKMerQueues) {
#ifdef DEB_MESS_SUPERSPLITTER
		// do timing stuff
		StopWatch sw(CLOCK_THREAD_CPUTIME_ID);
		sw.start();
#endif

		uint64 temp_sbc(0);
		uint64 temp_kmerc(0);
		uint64 temp_kmerc2(0);

		// temporary variables
		SuperBundle* sb = new SuperBundle();// the current super bundle in use
		KMer<K> kMer, iKMer;		// the current k-mer/ inverse k-mer in use
		const KMer<K>* nKMer;					// the current normalized k-mer
		uint16 l;								// size of current super s-mer
		byte* b;						// byte representation of current s-mer.
		size_t h;				// current hash value of k-mer (see distributor)
		byte nextBase;							// next base of current s-mer
		uint_tfn curTempFileId;	// id of temp file, note: each temp file has its own id

		// Kmer bundles
		cpu::KMerBundle<K>** cpuKMerBundles =
				new cpu::KMerBundle<K>*[_numCPUHasher];
		gpu::KMerBundle<K>** gpuKMerBundles =
				new gpu::KMerBundle<K>*[_numGPUHasher];

		// init Kmer bundles
		for (uint8_t i = 0; i < _numCPUHasher; i++)
			cpuKMerBundles[i] = new cpu::KMerBundle<K>();
		for (uint8_t i = 0; i < _numGPUHasher; i++)
			gpuKMerBundles[i] = new gpu::KMerBundle<K>();

#ifdef DEB_MESS_SUPERSPLITTER
		sw.hold();	// timing
#endif

		// take next superbundle out of superbundle queue
		_superBundleQueue->swapPop(sb);

		// timing
#ifdef DEB_MESS_SUPERSPLITTER
		sw.proceed();
#endif

		// for each temporary file
		for (uint rdyNmbs(0); rdyNmbs < _tempFilesNumber; ++rdyNmbs) {

			// not clear what is going on here?
			curTempFileId = _tempFilesOrder[rdyNmbs];

			// get number of kmers in this file
			uint64_t kmersInFile = _tempFiles[curTempFileId].getKMersNumber();

			/* Update the distributor that controls the ratio of
			 *  kmers distributed the various hash tables. */
			distributor->updateFileInformation(curTempFileId, kmersInFile);
			this->barrier->sync();

			// Some sort of progress bar
			if (threadId == 0)
				putchar('*');

			// if super bundle is not yet empty and belongs the the current file
			if (!sb->isEmpty() && sb->tempFileId == curTempFileId) {

#ifdef DEB_MESS_SUPERSPLITTER
				// timing stuff
				sw.hold();
#endif

				// while superbundle is not yet empty and next supermer belongs
				// to the current file: read next superbundle, jump over if next binId
				while ((!sb->isEmpty() || _superBundleQueue->swapPop(sb))
						&& sb->tempFileId == curTempFileId) {

					// timing
#ifdef DEB_MESS_SUPERSPLITTER
					sw.proceed();
#endif

					// while ?
					while (sb->next(b, l)) {
						//c += l - K + 1;

						// extract first kmer out of supermer and add to kmer bundle
						if (NORM) {
							// normalize first
							KMer<K>::set(b, kMer, iKMer);
							nKMer = &kMer.getNormalized(iKMer);
							ADD_KMER_TO_BUNDLE((*nKMer), curTempFileId)

						} else {
							// just add to bundle
							kMer.set(b);
							ADD_KMER_TO_BUNDLE(kMer, curTempFileId);
						}

						// extract all other kmers from supermer
						for (uint i = K; i < l; ++i) {
							nextBase = (*(b + (i >> 2))
									>> (6 - ((i & 0x3) << 1))) & 0x3;
							kMer.next(nextBase);

							// add to bundle
							if (NORM) {
								iKMer.nextInv(nextBase);
								nKMer = &kMer.getNormalized(iKMer);
								ADD_KMER_TO_BUNDLE((*nKMer), curTempFileId);
							} else {
								ADD_KMER_TO_BUNDLE(kMer, curTempFileId);
							}
						}
					}

					// clean up?
					sb->clear();

					// timing
#ifdef DEB_MESS_SUPERSPLITTER
					sw.hold();
#endif
				}

				// set current file id to kmer bundles and push them to the queues
				for (uint8_t i(0); i < _numCPUHasher; ++i) {
					if (!cpuKMerBundles[i]->isEmpty()) {
						cpuKMerBundles[i]->setTempFileId(curTempFileId);
						cpuKMerQueues[i]->swapPush(cpuKMerBundles[i]);
					}
				}

				// for gpu: set current file id to kmer bundles and push them to the queues
				for (uint8_t i(0); i < _numGPUHasher; ++i) {
					if (!gpuKMerBundles[i]->isEmpty()) {
						gpuKMerBundles[i]->setTempFileId(curTempFileId);
						gpuKMerQueues[i]->swapPush(gpuKMerBundles[i]);
					}
				}

				// timing
#ifdef DEB_MESS_SUPERSPLITTER
				sw.proceed();
#endif
			}

			//MEMORY_BARRIER
			memoryBarrier();

			this->barrier->sync();
		}

		// clean up
		delete sb;

		// free kmer bundles
		for (uint8_t i = 0; i < _numCPUHasher; i++)
			delete cpuKMerBundles[i];
		for (uint8_t i = 0; i < _numGPUHasher; i++)
			delete gpuKMerBundles[i];

		delete[] cpuKMerBundles;
		delete[] gpuKMerBundles;

// timing
#ifdef DEB_MESS_SUPERSPLITTER
		sw.stop();
		printf("superplitter[%2u]: %.3f s\n", threadId, sw.get_s());
#endif
	}

	/**
	 * Wait for all splitter threads to finalize.
	 */
	void joinSplitterThreads() {
		for (uint8_t i = 0; i < _processSplitterThreadsNumber; i++) {
			_processSplitterThreads[i]->join();
			delete _processSplitterThreads[i];
		}
	}

public:

	/**
	 * Constructor.
	 */
	KmerHasher(const uint32_t &k, const uint32 &kmcBundlesNumber,
			SyncSwapQueueSPMC<SuperBundle>* const superBundleQueue,
			const uint8 &processSplitterThreadsNumber,
			const uint8 &processHasherThreadsNumber,
			const uint8 &processGPUHasherThreadsNumber, TempFile* tempFiles,
			const uint_tfn &tempFilesNumber, const uint32 &thresholdMin,
			const bool &norm, std::string pTempFolder,
			const uint64 &maxKmcHashtableSize, const uint32 &kMerBundlesNumber,
			uint_tfn* tempFilesOrder) :
			_kmcSyncSwapQueue(kmcBundlesNumber), _processSplitterThreadsNumber(
					processSplitterThreadsNumber), _numCPUHasher(
					processHasherThreadsNumber), _numGPUHasher(
					processGPUHasherThreadsNumber), barrier(nullptr), _k(k), _tempFilesNumber(
					tempFilesNumber), _thresholdMin(thresholdMin), _norm(norm), _tempFolder(
					pTempFolder), _superBundleQueue(superBundleQueue), _tempFiles(
					tempFiles), _processThread(NULL), _kMersNumberCPU(0), _kMersNumberGPU(
					0), _uKMersNumberCPU(0), _uKMersNumberGPU(0), _btUKMersNumberCPU(
					0), _btUKMersNumberGPU(0), _maxKmcHashtableSize(
					maxKmcHashtableSize), _syncSplitterCounter(0), _kMerBundlesNumber(
					kMerBundlesNumber), _tempFilesOrder(tempFilesOrder)
	{

		// create array of threads
		_processSplitterThreads =
				new std::thread*[_processSplitterThreadsNumber];

		distributor = new KmerDistributer(_numCPUHasher, _numGPUHasher,
				tempFilesNumber);

		barrier = new Barrier(_processSplitterThreadsNumber);
	}

	/**
	 * Clean up.
	 */
	~KmerHasher() {
		delete[] _processSplitterThreads;
		delete distributor;
		delete barrier;
	}

	void saveHistogram() {
		FILE* file;
		file = fopen((_tempFolder + "histogram.csv").c_str(), "wb");
		fprintf(file, "counter; number of uk-mers\n");
		for(uint i = 1; i < HISTOGRAM_SIZE; ++i)
			fprintf(file, "  %3u; %9lu\n", i, _histogram[i]);
		fprintf(file, ">=%3u; %9lu\n", HISTOGRAM_SIZE, _histogram[0]);
		fclose(file);
	}

	/** Main working procedure for this class. */
	void process() {

		// decide which template specialization to use
		switch (_k) {
		case 8:
			process_template<8>();
			break;
		case 9:
			process_template<9>();
			break;
		case 10:
			process_template<10>();
			break;
		case 11:
			process_template<11>();
			break;
		case 12:
			process_template<12>();
			break;
		case 13:
			process_template<13>();
			break;
		case 14:
			process_template<14>();
			break;
		case 15:
			process_template<15>();
			break;
		case 16:
			process_template<16>();
			break;
		case 17:
			process_template<17>();
			break;
		case 18:
			process_template<18>();
			break;
		case 19:
			process_template<19>();
			break;
		case 20:
			process_template<20>();
			break;
		case 21:
			process_template<21>();
			break;
		case 22:
			process_template<22>();
			break;
		case 23:
			process_template<23>();
			break;
		case 24:
			process_template<24>();
			break;
		case 25:
			process_template<25>();
			break;
		case 26:
			process_template<26>();
			break;
		case 27:
			process_template<27>();
			break;
		case 28:
			process_template<28>();
			break;
		case 29:
			process_template<29>();
			break;
		case 30:
			process_template<30>();
			break;
		case 31:
			process_template<31>();
			break;
		case 32:
			process_template<32>();
			break;
		case 33:
			process_template<33>();
			break;
		case 34:
			process_template<34>();
			break;
		case 35:
			process_template<35>();
			break;
		case 36:
			process_template<36>();
			break;
		case 37:
			process_template<37>();
			break;
		case 38:
			process_template<38>();
			break;
		case 39:
			process_template<39>();
			break;
		case 40:
			process_template<40>();
			break;
		case 41:
			process_template<41>();
			break;
		case 42:
			process_template<42>();
			break;
		case 43:
			process_template<43>();
			break;
		case 44:
			process_template<44>();
			break;
		case 45:
			process_template<45>();
			break;
		case 46:
			process_template<46>();
			break;
		case 47:
			process_template<47>();
			break;
		case 48:
			process_template<48>();
			break;
		case 49:
			process_template<49>();
			break;
		case 50:
			process_template<50>();
			break;
		case 51:
			process_template<51>();
			break;
		case 52:
			process_template<52>();
			break;
		case 53:
			process_template<53>();
			break;
		case 54:
			process_template<54>();
			break;
		case 55:
			process_template<55>();
			break;
		case 56:
			process_template<56>();
			break;
		case 57:
			process_template<57>();
			break;
		case 58:
			process_template<58>();
			break;
		case 59:
			process_template<59>();
			break;
		case 60:
			process_template<60>();
			break;
		case 61:
			process_template<61>();
			break;
		case 62:
			process_template<62>();
			break;
		case 63:
			process_template<63>();
			break;
		case 64:
			process_template<64>();
			break;
		case 65:
			process_template<65>();
			break;
		case 66:
			process_template<66>();
			break;
		case 67:
			process_template<67>();
			break;
		case 68:
			process_template<68>();
			break;
		case 69:
			process_template<69>();
			break;
		case 70:
			process_template<70>();
			break;
		case 71:
			process_template<71>();
			break;
		case 72:
			process_template<72>();
			break;
		case 73:
			process_template<73>();
			break;
		case 74:
			process_template<74>();
			break;
		case 75:
			process_template<75>();
			break;
		case 76:
			process_template<76>();
			break;
		case 77:
			process_template<77>();
			break;
		case 78:
			process_template<78>();
			break;
		case 79:
			process_template<79>();
			break;
		case 80:
			process_template<80>();
			break;
		case 81:
			process_template<81>();
			break;
		case 82:
			process_template<82>();
			break;
		case 83:
			process_template<83>();
			break;
		case 84:
			process_template<84>();
			break;
		case 85:
			process_template<85>();
			break;
		case 86:
			process_template<86>();
			break;
		case 87:
			process_template<87>();
			break;
		case 88:
			process_template<88>();
			break;
		case 89:
			process_template<89>();
			break;
		case 90:
			process_template<90>();
			break;
		case 91:
			process_template<91>();
			break;
		case 92:
			process_template<92>();
			break;
		case 93:
			process_template<93>();
			break;
		case 94:
			process_template<94>();
			break;
		case 95:
			process_template<95>();
			break;
		case 96:
			process_template<96>();
			break;
		case 97:
			process_template<97>();
			break;
		case 98:
			process_template<98>();
			break;
		case 99:
			process_template<99>();
			break;
		case 100:
			process_template<100>();
			break;
		case 101:
			process_template<101>();
			break;
		case 102:
			process_template<102>();
			break;
		case 103:
			process_template<103>();
			break;
		case 104:
			process_template<104>();
			break;
		case 105:
			process_template<105>();
			break;
		case 106:
			process_template<106>();
			break;
		case 107:
			process_template<107>();
			break;
		case 108:
			process_template<108>();
			break;
		case 109:
			process_template<109>();
			break;
		case 110:
			process_template<110>();
			break;
		case 111:
			process_template<111>();
			break;
		case 112:
			process_template<112>();
			break;
		case 113:
			process_template<113>();
			break;
		case 114:
			process_template<114>();
			break;
		case 115:
			process_template<115>();
			break;
		case 116:
			process_template<116>();
			break;
		case 117:
			process_template<117>();
			break;
		case 118:
			process_template<118>();
			break;
		case 119:
			process_template<119>();
			break;
		case 120:
			process_template<120>();
			break;
		case 121:
			process_template<121>();
			break;
		case 122:
			process_template<122>();
			break;
		case 123:
			process_template<123>();
			break;
		case 124:
			process_template<124>();
			break;
		case 125:
			process_template<125>();
			break;
		case 126:
			process_template<126>();
			break;
		case 127:
			process_template<127>();
			break;
		case 128:
			process_template<128>();
			break;
		case 129:
			process_template<129>();
			break;
		case 130:
			process_template<130>();
			break;
		case 131:
			process_template<131>();
			break;
		case 132:
			process_template<132>();
			break;
		case 133:
			process_template<133>();
			break;
		case 134:
			process_template<134>();
			break;
		case 135:
			process_template<135>();
			break;
		case 136:
			process_template<136>();
			break;
		case 137:
			process_template<137>();
			break;
		case 138:
			process_template<138>();
			break;
		case 139:
			process_template<139>();
			break;
		case 140:
			process_template<140>();
			break;
		case 141:
			process_template<141>();
			break;
		case 142:
			process_template<142>();
			break;
		case 143:
			process_template<143>();
			break;
		case 144:
			process_template<144>();
			break;
		case 145:
			process_template<145>();
			break;
		case 146:
			process_template<146>();
			break;
		case 147:
			process_template<147>();
			break;
		case 148:
			process_template<148>();
			break;
		case 149:
			process_template<149>();
			break;
		case 150:
			process_template<150>();
			break;
		case 151:
			process_template<151>();
			break;
		case 152:
			process_template<152>();
			break;
		case 153:
			process_template<153>();
			break;
		case 154:
			process_template<154>();
			break;
		case 155:
			process_template<155>();
			break;
		case 156:
			process_template<156>();
			break;
		case 157:
			process_template<157>();
			break;
		case 158:
			process_template<158>();
			break;
		case 159:
			process_template<159>();
			break;
		case 160:
			process_template<160>();
			break;
		case 161:
			process_template<161>();
			break;
		case 162:
			process_template<162>();
			break;
		case 163:
			process_template<163>();
			break;
		case 164:
			process_template<164>();
			break;
		case 165:
			process_template<165>();
			break;
		case 166:
			process_template<166>();
			break;
		case 167:
			process_template<167>();
			break;
		case 168:
			process_template<168>();
			break;
		case 169:
			process_template<169>();
			break;
		case 170:
			process_template<170>();
			break;
		case 171:
			process_template<171>();
			break;
		case 172:
			process_template<172>();
			break;
		case 173:
			process_template<173>();
			break;
		case 174:
			process_template<174>();
			break;
		case 175:
			process_template<175>();
			break;
		case 176:
			process_template<176>();
			break;
		case 177:
			process_template<177>();
			break;
		case 178:
			process_template<178>();
			break;
		case 179:
			process_template<179>();
			break;
		case 180:
			process_template<180>();
			break;
		case 181:
			process_template<181>();
			break;
		case 182:
			process_template<182>();
			break;
		case 183:
			process_template<183>();
			break;
		case 184:
			process_template<184>();
			break;
		case 185:
			process_template<185>();
			break;
		case 186:
			process_template<186>();
			break;
		case 187:
			process_template<187>();
			break;
		case 188:
			process_template<188>();
			break;
		case 189:
			process_template<189>();
			break;
		case 190:
			process_template<190>();
			break;
		case 191:
			process_template<191>();
			break;
		case 192:
			process_template<192>();
			break;
		case 193:
			process_template<193>();
			break;
		case 194:
			process_template<194>();
			break;
		case 195:
			process_template<195>();
			break;
		case 196:
			process_template<196>();
			break;
		case 197:
			process_template<197>();
			break;
		case 198:
			process_template<198>();
			break;
		case 199:
			process_template<199>();
			break;
		case 200:
			process_template<200>();
			break;
		case 201:
			process_template<201>();
			break;
		case 202:
			process_template<202>();
			break;
		case 203:
			process_template<203>();
			break;
		case 204:
			process_template<204>();
			break;
		case 205:
			process_template<205>();
			break;
		case 206:
			process_template<206>();
			break;
		case 207:
			process_template<207>();
			break;
		case 208:
			process_template<208>();
			break;
		case 209:
			process_template<209>();
			break;
		case 210:
			process_template<210>();
			break;
		case 211:
			process_template<211>();
			break;
		case 212:
			process_template<212>();
			break;
		case 213:
			process_template<213>();
			break;
		case 214:
			process_template<214>();
			break;
		case 215:
			process_template<215>();
			break;
		case 216:
			process_template<216>();
			break;
		case 217:
			process_template<217>();
			break;
		case 218:
			process_template<218>();
			break;
		case 219:
			process_template<219>();
			break;
		case 220:
			process_template<220>();
			break;
		case 221:
			process_template<221>();
			break;
		case 222:
			process_template<222>();
			break;
		case 223:
			process_template<223>();
			break;
		case 224:
			process_template<224>();
			break;
		case 225:
			process_template<225>();
			break;
		case 226:
			process_template<226>();
			break;
		case 227:
			process_template<227>();
			break;
		case 228:
			process_template<228>();
			break;
		case 229:
			process_template<229>();
			break;
		case 230:
			process_template<230>();
			break;
		case 231:
			process_template<231>();
			break;
		case 232:
			process_template<232>();
			break;
		case 233:
			process_template<233>();
			break;
		case 234:
			process_template<234>();
			break;
		case 235:
			process_template<235>();
			break;
		case 236:
			process_template<236>();
			break;
		case 237:
			process_template<237>();
			break;
		case 238:
			process_template<238>();
			break;
		case 239:
			process_template<239>();
			break;
		case 240:
			process_template<240>();
			break;
		case 241:
			process_template<241>();
			break;
		case 242:
			process_template<242>();
			break;
		case 243:
			process_template<243>();
			break;
		case 244:
			process_template<244>();
			break;
		case 245:
			process_template<245>();
			break;
		case 246:
			process_template<246>();
			break;
		case 247:
			process_template<247>();
			break;
		case 248:
			process_template<248>();
			break;
		case 249:
			process_template<249>();
			break;
		case 250:
			process_template<250>();
			break;
		case 251:
			process_template<251>();
			break;
		case 252:
			process_template<252>();
			break;
		case 253:
			process_template<253>();
			break;
		case 254:
			process_template<254>();
			break;
		case 255:
			process_template<255>();
			break;
		case 256:
			process_template<256>();
			break;
		case 257:
			process_template<257>();
			break;
		case 258:
			process_template<258>();
			break;
		case 259:
			process_template<259>();
			break;
		case 260:
			process_template<260>();
			break;
		case 261:
			process_template<261>();
			break;
		case 262:
			process_template<262>();
			break;
		case 263:
			process_template<263>();
			break;
		case 264:
			process_template<264>();
			break;
		case 265:
			process_template<265>();
			break;
		case 266:
			process_template<266>();
			break;
		case 267:
			process_template<267>();
			break;
		case 268:
			process_template<268>();
			break;
		case 269:
			process_template<269>();
			break;
		case 270:
			process_template<270>();
			break;
		case 271:
			process_template<271>();
			break;
		case 272:
			process_template<272>();
			break;
		case 273:
			process_template<273>();
			break;
		case 274:
			process_template<274>();
			break;
		case 275:
			process_template<275>();
			break;
		case 276:
			process_template<276>();
			break;
		case 277:
			process_template<277>();
			break;
		case 278:
			process_template<278>();
			break;
		case 279:
			process_template<279>();
			break;
		case 280:
			process_template<280>();
			break;
		case 281:
			process_template<281>();
			break;
		case 282:
			process_template<282>();
			break;
		case 283:
			process_template<283>();
			break;
		case 284:
			process_template<284>();
			break;
		case 285:
			process_template<285>();
			break;
		case 286:
			process_template<286>();
			break;
		case 287:
			process_template<287>();
			break;
		case 288:
			process_template<288>();
			break;
		case 289:
			process_template<289>();
			break;
		case 290:
			process_template<290>();
			break;
		case 291:
			process_template<291>();
			break;
		case 292:
			process_template<292>();
			break;
		case 293:
			process_template<293>();
			break;
		case 294:
			process_template<294>();
			break;
		case 295:
			process_template<295>();
			break;
		case 296:
			process_template<296>();
			break;
		case 297:
			process_template<297>();
			break;
		case 298:
			process_template<298>();
			break;
		case 299:
			process_template<299>();
			break;
		case 300:
			process_template<300>();
			break;
		default:
			throw std::runtime_error(
					std::string("Gerbil Error: Unsupported k"));
		}
	}

	void join() {
		_processThread->join();
		delete _processThread;
		IF_DEB(printf("all KmerHashers are rdy...\n"));
	}

	void print() {
		printf("kmers (CPU)     : %12lu\n", _kMersNumberCPU);
		printf("kmers (GPU)     : %12lu\n", _kMersNumberGPU);
		printf("ukmers (CPU)    : %12lu\n", _uKMersNumberCPU);
		printf("ukmers (GPU)    : %12lu\n", _uKMersNumberGPU);
		printf("below th (CPU)  : %12lu\n", _btUKMersNumberCPU);
		printf("below th (GPU)  : %12lu\n", _btUKMersNumberGPU);
	}

	SyncSwapQueueMPSC<KmcBundle> *getKmcSyncSwapQueue() {
		return &_kmcSyncSwapQueue;
	}
}
;

}

#endif /* KMERHASHER_H_ */

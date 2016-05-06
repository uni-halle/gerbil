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

#ifndef KMERCOUNTINGHASHTABLE_H_
#define KMERCOUNTINGHASHTABLE_H_

#include "../cuda_ds/CountingHashTable.h"
#include "KMer.h"
#include "SyncQueue.h"
#include "Bundle.h"

namespace gerbil {
namespace gpu {

namespace cd = ::cuda_ds;
namespace cdi = ::cuda_ds::internal;

/**
 * Extend the cuda_ds::CountingHashTable to specialize the counting of KMers.
 * Just wrapped the constructor, the add method and added an additional
 * extracting method.
 *
 * @param K K
 * @param intsPerKMer The number of ints per KMer
 */
template<uint32_t K>
class KmerCountingHashTable: public cd::CountingHashTable<KMer<K>,
GPU_KMER_BUNDLE_DATA_SIZE> {

private:
	typedef cdi::KeyValuePair<cdi::intsPerKey<KMer<K>>()> _KeyValuePair;

	SyncSwapQueueMPSC<KmcBundle>* kmcQueue;	// output queue
	uint32_t threshold;						// minimal number of kmer occurence
	uint64_t kMersNumber;					// number of inserted kmers
	uint64_t uKMersNumber;					// number of inserted unique kmers
	uint64_t btUKMersNumber;// number of inserted unique kmers below threshold

	char* rawbuffer;						// buffer for copy jobs

	// buffers needed for emergency evacuation of the noSuccessArea of the Hash table
	FailureBuffer<K>* failureBuffer;
	FailureBuffer<K>* tempBuffer;
	KMerBundle<K>* _kmb;					// working kmer bundle
	KmcBundle* curKmcBundle;			// working kmc bundle

public:
	/**
	 * Main Constructor.
	 * @param devID the cuda device id on which the table should be allocated.
	 * @param kmcQueue The output queue where kmer counts are stored.
	 * @param threshold The minimal number of counts to be output.
	 */
	KmerCountingHashTable(const uint32_t devID,
			SyncSwapQueueMPSC<KmcBundle>* kmcQueue, uint32_t threshold,
			std::string tempPath) :
			cd::CountingHashTable<KMer<K>, GPU_KMER_BUNDLE_DATA_SIZE>(devID), kmcQueue(
					kmcQueue), threshold(threshold), kMersNumber(0), uKMersNumber(
					0), btUKMersNumber(0), rawbuffer(nullptr), failureBuffer(
					nullptr), tempBuffer(nullptr), _kmb(nullptr), curKmcBundle(
					nullptr) {

		// allocate host memory buffer
		rawbuffer = (char*) malloc(GPU_COPY_BUFFER_SIZE);
		failureBuffer = new FailureBuffer<K>(
		FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, tempPath, 1024 + devID,
				0);
		tempBuffer = new FailureBuffer<K>(
		FAILUREBUFFER_KMER_BUNDLES_NUMBER_PER_THREAD, tempPath, 1024 + devID,
				1);
		_kmb = new KMerBundle<K>();
		curKmcBundle = new KmcBundle();
	}

	~KmerCountingHashTable() {
		free(rawbuffer);
		delete failureBuffer;
		delete tempBuffer;
		delete _kmb;
		delete curKmcBundle;
	}

	/**
	 * Get number of kmers inserted until last call of extract().
	 */
	uint64_t getKMersNumber() const {
		return kMersNumber;
	}

	/**
	 * Get number of distinct kmers inserted until last call of extract().
	 */
	uint64_t getUKMersNumber() const {
		return uKMersNumber;
	}

	/**
	 * Get number of distinct kmers below the threshold until last call of extract().
	 */
	uint64_t getUKMersNumberBelowThreshold() const {
		return btUKMersNumber;
	}

	/**
	 * Add a bundle of kmers to the table.
	 * @param kmb Pointer to the kmer bundle to be inserted.
	 * @return True, if the kmers in this bundle could be stored successfully.
	 * False, if an emergency extract is neccessary.
	 */
	void addBundle(gpu::KMerBundle<K>* kmb) {

		// try to insert the bundle into the hash table
		bool success = this->insert((const KMer<K>*) kmb->getData(),
				kmb->count());

		if (!success) {
			printf("EMERGENCY!\n");
			emergencyExtract();
			// try to add again
			success = this->insert((const KMer<K>*) kmb->getData(),
					kmb->count());
		}

		kmb->clear();
	}

	/**
	 * Extract the table entries and push them into the kmc bundle queue.
	 */
	void extractTableEntries() {

		// set device id
		cudaSetDevice(this->devID);

		// entry-based view on the buffer
		_KeyValuePair *kmerCounts = (_KeyValuePair*) rawbuffer;

		// temporary working variables
		KMer<K>* curKey;
		_KeyValuePair* curKeyValuePair;
		uint32_t curCount;

		// Load table to host memory in small chunks of size GPU_COPY_BUFFER_SIZE
		const uint64_t sizeOfEntry = sizeof(_KeyValuePair);
		const uint32_t entriesPerBuffer = GPU_COPY_BUFFER_SIZE / sizeOfEntry;

		// wait for threads to success
		cudaError_t err = cudaDeviceSynchronize();
		if (err != cudaSuccess) {
			throw std::runtime_error(
					std::string("gerbil::exception: ")
							+ std::string(cudaGetErrorName(err)) + " "
							+ std::string(cudaGetErrorString(err)));
		}

		// compress table entries
		//const uint32_t numDifferentKeys = this->compressEntries();

		// take back all memory from table
		for (uint32_t i = 0; i * entriesPerBuffer < this->getNumEntries();
				i++) {

			const uint32_t itemsToFetch =
					(i + 1) * entriesPerBuffer <= this->getNumEntries() ?
							entriesPerBuffer :
							this->getNumEntries() % entriesPerBuffer;

			// copy chunk of data
			cudaMemcpy(rawbuffer, &this->table[i * entriesPerBuffer],
					itemsToFetch * sizeOfEntry, cudaMemcpyDeviceToHost);

			// copy kmer counts to result vector
			for (uint32_t j = 0; j < itemsToFetch; j++) {

				curKeyValuePair = (_KeyValuePair*) &kmerCounts[j];
				curKey = (KMer<K>*) &curKeyValuePair->getKey();
				curCount = curKeyValuePair->getCount();

				if (!curKey->isEmpty() && curCount > 0) {

					if (curCount >= threshold) {

						// copy kmer count to output queue
						if (!curKmcBundle->add<K>(*curKey, curCount)) {
							kmcQueue->swapPush(curKmcBundle);
							curKmcBundle->add<K>(*curKey, curCount);
						}
						curKey->clear();
					} else {
						btUKMersNumber++;
					}

					uKMersNumber++;
					kMersNumber += curCount;
				}
			}
		}
	}

	void emergencyExtract() {

		// set device id
		cudaSetDevice(this->devID);

		// entry-based view on the buffer
		KMer<K>* kmers = (KMer<K>*) rawbuffer;

		// Determine if there are any unsuccessfully inserted kmers
		uint32_t numNoSuccess;
		cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint32_t),
				cudaMemcpyDeviceToHost);

		printf("copy %u kmers back to host\n", numNoSuccess);

		// if there are some
		if (numNoSuccess > 0) {

			// Load kmers (in small chunks) to host memory
			const uint32_t kmersPerBuffer = GPU_COPY_BUFFER_SIZE
					/ sizeof(KMer<K> );

			// for all chunks
			for (uint32_t i = 0; i * kmersPerBuffer < numNoSuccess; i++) {

				const uint32_t itemsToFetch =
						(i + 1) * kmersPerBuffer <= numNoSuccess ?
								kmersPerBuffer : numNoSuccess % kmersPerBuffer;

				// copy back to host memory
				cudaMemcpy(kmers, &this->noSuccessArea[i * kmersPerBuffer],
						itemsToFetch * sizeof(KMer<K> ),
						cudaMemcpyDeviceToHost);

				// store kmers in failure buffer
				for (uint32_t j = 0; j < itemsToFetch; j++) {
					failureBuffer->addKMer(kmers[j]);
				}
			}

		}

		// reset the no Success Area
		cudaError_t err = cudaMemset(this->numNoSuccessPtr, 0,
				sizeof(uint32_t));
		if (err != cudaSuccess) {
			throw std::runtime_error(
					std::string("cuda_ds::exception: ")
							+ std::string(cudaGetErrorName(err)) + " "
							+ std::string(cudaGetErrorString(err)));
		}

		printf("done\n");
	}

	/**
	 * Extract, compress and output the entries stored at the noSuccessArea
	 * and push them as kmer counts into the kmc bundle queue.
	 */
	void extractNoSuccessArea() {

		// set device id
		cudaSetDevice(this->devID);

		// entry-based view on the buffer
		KMer<K>* kmers = (KMer<K>*) rawbuffer;

		// working variables
		uint32_t curCount;

		// Determine if there are any unsuccessfully inserted kmers
		uint32_t numNoSuccess;
		cudaMemcpy(&numNoSuccess, this->numNoSuccessPtr, sizeof(uint32_t),
				cudaMemcpyDeviceToHost);

		// better assure that maxNumNoSuccess has not been reached
		numNoSuccess = std::min(numNoSuccess, this->maxNumNoSuccess);

		// if there are some
		if (numNoSuccess > 0) {

			printf("numNoSuccess=%u\n", numNoSuccess);

			// sort area of shame
			cdi::sortKeys<cdi::intsPerKey<KMer<K>>()>(this->noSuccessArea,
					numNoSuccess);

			// Load kmers (in small chunks)
			const uint32_t kmersPerBuffer = GPU_COPY_BUFFER_SIZE
					/ sizeof(KMer<K> );

			// load first buffer
			uint32_t itemsToFetch =
					kmersPerBuffer <= numNoSuccess ?
							kmersPerBuffer : numNoSuccess % kmersPerBuffer;

			cudaMemcpy(rawbuffer, this->noSuccessArea,
					itemsToFetch * sizeof(KMer<K> ), cudaMemcpyDeviceToHost);

			// compare with last seen key
			KMer<K> lastKey = kmers[0];
			curCount = 1;

			// compress and add to result vector
			for (uint32_t j = 1; j < itemsToFetch; j++) {

				if (kmers[j] != lastKey) {

					if (!lastKey.isEmpty()) {

						if (curCount >= threshold) {

							// copy kmer count to output queue
							if (!curKmcBundle->add<K>(lastKey, curCount)) {
								kmcQueue->swapPush(curKmcBundle);
								curKmcBundle->add<K>(lastKey, curCount);
							}
						} else {
							btUKMersNumber++;
						}

						uKMersNumber++;
						kMersNumber += curCount;
					}

					lastKey = kmers[j];
					curCount = 1;
				} else
					curCount++;
			}

			// load other buffers (if any)
			for (uint32_t i = 1; i * kmersPerBuffer < numNoSuccess; i++) {

				itemsToFetch =
						(i + 1) * kmersPerBuffer <= numNoSuccess ?
								kmersPerBuffer : numNoSuccess % kmersPerBuffer;

				cudaMemcpy(rawbuffer, &this->noSuccessArea[i * kmersPerBuffer],
						itemsToFetch * sizeof(KMer<K> ),
						cudaMemcpyDeviceToHost);

				// compress and add to result vector
				for (uint32_t j = 0; j < itemsToFetch; j++) {
					if (kmers[j] != lastKey) {

						// put out kmer count
						if (!lastKey.isEmpty()) {

							if (curCount >= threshold) {

								// copy kmer count to output queue
								if (!curKmcBundle->add<K>(lastKey, curCount)) {
									kmcQueue->swapPush(curKmcBundle);
									curKmcBundle->add<K>(lastKey, curCount);
								}
							} else {
								btUKMersNumber++;
							}

							uKMersNumber++;
							kMersNumber += curCount;
						}

						lastKey = kmers[j];
						curCount = 1;
					} else
						curCount++;
				}
			}

			// insert last seen kmer
			if (!lastKey.isEmpty()) {

				if (curCount >= threshold) {
					// copy kmer count to output queue
					if (!curKmcBundle->add<K>(lastKey, curCount)) {
						kmcQueue->swapPush(curKmcBundle);
						curKmcBundle->add<K>(lastKey, curCount);
					}
				} else {
					btUKMersNumber++;
				}

				uKMersNumber++;
				kMersNumber += curCount;
			}
		}
	}

	/**
	 * Extract kmer counts from the hash table and insert them
	 * as kmc bundles into the output queue. In the end, the
	 * table is empty again and has same size.
	 */
	void extractAndClear() {

		extractTableEntries();
		extractNoSuccessArea();

		// clear the table
		this->clear();

		// handle failures as long as some exist
		while (!failureBuffer->isEmpty()) {

			// the original failureBuffer may be filled again
			std::swap(tempBuffer, failureBuffer);

			// for all kmer bundles in failure buffer
			while (tempBuffer->getNextKMerBundle(_kmb)) {
				// add bundle to the table
				addBundle(_kmb);
			}

			// extract kmer counts
			extractTableEntries();
			extractNoSuccessArea();

			this->clear();

			// all kmer are read from failure buffer
			tempBuffer->clear();
		}

		// submit last kmc bundle
		if (!curKmcBundle->isEmpty()) {
			kmcQueue->swapPush(curKmcBundle);
		}
	}

};

}
}

#endif /* KMERCOUNTINGHASHTABLE_H_ */

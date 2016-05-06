/*
 * KmerDistributer.cpp
 *
 *  Created on: Jan 22, 2016
 *      Author: rechner
 */

#include "../../include/gerbil/KmerDistributer.h"
#include <stdio.h>
#include <cmath>
#include <cstring>
#include <climits>
#include <map>
#include <algorithm>

namespace gerbil {

KmerDistributer::KmerDistributer(const uint32_t numCPUHasherThreads,
		const uint32_t numGPUHasherThreads, const uint32_t maxNumFiles) :
		numThreadsCPU(numCPUHasherThreads), numThreadsGPU(numGPUHasherThreads) {

	const uint32_t numThreads = numCPUHasherThreads + numGPUHasherThreads;

	prefixSum = new float[numThreads + 1];
	tmp = new float[numThreads];
	capacity = new uint64_t[numThreads];

	throughput = new float[numThreads];
	memset(throughput, 0, numThreads * sizeof(float)); // fill throughput array with initial zeros

	// init all capacities with infinity
	for (int i = 0; i < numThreads; i++)
		capacity[i] = ULONG_MAX;

	// initialize file specific variables
	ratio.resize(maxNumFiles);
	buckets.resize(maxNumFiles);
	lock.resize(maxNumFiles);
	for (int i = 0; i < maxNumFiles; i++) {
		ratio[i] = new float[numThreads];
		buckets[i] = new uint32_t[BUCKETSIZE];
		lock[i] = 1;				// set lock status of each file to locked
		memset(ratio[i], 0, numThreads * sizeof(float));
	}
}

KmerDistributer::~KmerDistributer() {
	delete[] prefixSum;
	delete[] tmp;
	delete[] capacity;
	delete[] throughput;
	for (int i = 0; i < buckets.size(); i++) {
		delete[] buckets[i];
		delete[] ratio[i];
	}
}

void KmerDistributer::updateThroughput(const bool gpu, const uint32_t tId,
		const float t) {

	if (gpu) {
		//printf("gpu hasher %u: update throughput to %f kmers/ms\n", tId, t);
		throughput[numThreadsCPU + tId] = t;
	} else {
		//printf("cpu hasher %u: update throughput to %f kmers/ms\n", tId, t);
		throughput[tId] = t;
	}
}

void KmerDistributer::updateCapacity(const bool gpu, const uint32_t tId,
		const uint64_t cap) {
	if (gpu) {
		//printf("gpu hasher %u: update capacity to %u kmer\n", tId, cap);
		capacity[numThreadsCPU + tId] = cap;
	} else {
		//printf("cpu hasher %u: update capacity to %u kmer\n", tId, cap);
		capacity[tId] = cap;
	}
}

float KmerDistributer::getSplitRatio(const bool gpu, const uint32_t tId,
		const uint_tfn curTempFileId) const {

	//return 1;

	// wait until updates are complete
	while (lock[curTempFileId])
		usleep(10);

	return gpu ?
			ratio[curTempFileId][numThreadsCPU + tId] :
			ratio[curTempFileId][tId];
}

void KmerDistributer::updateFileInformation(const uint32_t curFileId,
		const uint64_t kmersInFile) {

	// try to lock data structure
	if(CAS(&lock[curFileId], 1, 2) != 1)
		return;

	const float eps = 1e-3;
	const uint32_t numThreads = numThreadsCPU + numThreadsGPU;

	// compute sum of all throughput
	float sum = 0;
	for (int i = 0; i < numThreads; i++) {
		//printf("thread %u has throughput %f\n", i, throughput[i]);
		sum += throughput[i];

		// break if any thread has not yet provieded throughput information
		if (fabs(throughput[i]) < eps) {
			sum = 0;
			break;
		}
	}

	/**
	 * Update the ratio of kmers from latest throughput and capacity
	 */

	// if some thread has not provided throughput information
	if (fabs(sum) < eps) {
		// distribute kmers equally over all threads
		for (int i = 0; i < numThreads; i++) {
			ratio[curFileId][i] = 1.0f / numThreads;
			//printf("thread %u becomes ratio of %f\n", i, ratio[curFileId][i]);
		}
	}
	// if all threads have provided throughput information
	else {

		// define the new ratio of kmers for each hasher thread
		for (int i = 0; i < numThreads; i++) {
			/* Version 1: Each Ratio is defined as average between throughput ratio and uniform distribution */
			//cpuRatio[i] = 0.5 * cpu_throughput[i] / sum + 0.5 / (numCPUHasherThreads + numGPUHasherThreads);
			/* Version 2: Each Ratio is defined as average between troughput ratio and old ratio */
			ratio[curFileId][i] = 0.5 * ratio[lastFileId][i]
					+ 0.5 * throughput[i] / sum;

			/* Version 3: Each Ratio is defined as pure throughput ratio */
			//ratio[curFileId][i] = throughput[lastFileId][i] / sum;
			//printf("thread %u becomes ratio of %f\n", i, ratio[curFileId][i]);
		}
	}

	// Check whether any capacity constraint will be violated.

	int numPositiveCapacity = 0;
	for (int i = 0; i < numThreads; i++) {
		float maxRatio = (float) capacity[i] / (float) kmersInFile;
		tmp[i] = maxRatio - ratio[curFileId][i];
		if (tmp[i] > 0)
			numPositiveCapacity++;
		//printf("thread %u has diff of %f\n", i, tmp[i]);
	}

	for (int i = 0; i < numThreads; i++) {

		// while ratio constraint of thread i is violated and there are threads that could compensate
		while (tmp[i] < -eps && numPositiveCapacity > 0) {

			// try to redistribute an amount of kmers to other hasher threads
			/*printf(
					"capacity of thread %u is violated! redistribute a relative amount of %f kmers!\n",
					i, -tmp[i]);*/

			// try to redistribute uniformly to all threads with positive capacity
			float y = -tmp[i] / (float) numPositiveCapacity;

			for (int k = 0; k < numThreads; k++) {

				// if thread k has still open capcacity
				if (tmp[k] > eps) {

					float x = std::min(y, tmp[k]);
					ratio[curFileId][i] -= x;// ratio of thread i becomes smaller by x
					ratio[curFileId][k] += x;// ratio of thread k becomes larger by x
					tmp[i] += x;
					tmp[k] -= x;
					//printf("give %f to thread %u\n", x, k);
					if (tmp[k] < eps)
						numPositiveCapacity--;
				}
			}
		}

		// if difference is still smaller than zero: exception!
		if (numPositiveCapacity == 0 && tmp[i] > 0)
			throw std::runtime_error(
					"gerbil::exception: system has not enough memory for this number of files. Try to increase the number of temporary files.\n");
	}

	/**
	 *  Update prefix sums
	 */
	int j = 0;
	float x = 0;
	for (int i = 0; i < numThreads; i++) {
		prefixSum[j] = x;
		x += ratio[curFileId][i];
		j++;
	//	printf("thread %u becomes final ratio of %f\n", i, ratio[curFileId][i]);
	}
	prefixSum[numThreads] = 1.0f;

	/**
	 * Update Bucket table for the current file.
	 */
	int bucketId = 0;
	int i = 0;
	// for each hasher thread
	while (bucketId < numThreads) {

		// for each bucket entry with corresponds to current hasher thread
		while (i < BUCKETSIZE
				&& (float) i / (float) BUCKETSIZE <= prefixSum[bucketId + 1]) {
			buckets[curFileId][i] = bucketId;
			i++;
		}
		bucketId++;
	}

	// unlock current file
	lock[curFileId] = 0;

	lastFileId = curFileId;
}

} /* namespace gerbil */

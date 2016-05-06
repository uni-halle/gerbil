#include <iostream>
#include <map>
#include <algorithm>
#include <chrono>
#include <unordered_map>

#include "../../include/cuda_ds/CountingHashTable.h"
#include "Kmer.hpp"

using namespace std;
using namespace cuda_ds;

namespace std {

template<uint32_t k>
struct hash<Kmer<k>> {
	std::size_t operator()(const Kmer<k>& kmer) const {

		const int intsPerKmer = (k + 15) / 16;
		const int* kmerInts = (const int*) &kmer;
		size_t res = 0;
		for (int i = 0; i < intsPerKmer; i++) {
			res += 31 * kmerInts[i];
		}
		return res;
	}
};

template<uint32_t k>
struct key_comparer {
	bool operator()(const pair<Kmer<k>, uint32_t>& l, const pair<Kmer<k>, uint32_t>& r) const {
		return l.first < r.first;
	}
};

template<uint32_t k>
struct value_comparer {
	bool operator()(const pair<Kmer<k>, uint32_t>& l, const pair<Kmer<k>, uint32_t>& r) const {
		return l.second < r.second;
	}
};

}

int main2() {

	const uint32_t k = 28;
	const uint32_t numKmers = 10000000;
	const uint32_t bundleBufferSize = 512*1024;
	const uint32_t kmersPerBundle = bundleBufferSize / sizeof(Kmer<k> );
	const uint32_t numBundles = (numKmers + kmersPerBundle - 1)
			/ kmersPerBundle;
	srand(40);

	printf("k                = %i\n", k);
	printf("sizeof(Kmer<k>)  = %i\n", sizeof(Kmer<k>));
	printf("numKmers         = %i\n", numKmers);
	printf("bundleBufferSize = %i\n", bundleBufferSize);
	printf("kmersPerBundle   = %i\n", kmersPerBundle);
	printf("numBundles       = %i\n", numBundles);
	printf("\n");

	Kmer<k> **bundles = new Kmer<k>*[numBundles];
	uint32_t* bundleSizes = new uint32_t[numBundles];

	for (int i = 0; i < numBundles; i++) {
		bundles[i] = new Kmer<k> [kmersPerBundle];
	}

	const char proteins[4] = { 'A', 'C', 'G', 'T' };

	for (int l = 0; l < 100; l++) {

		cuda_ds::CountingHashTable<Kmer<k>, bundleBufferSize> table;

		unordered_map<Kmer<k>, int> mymap;
		//map<Kmer<k>, int> mymap;

		vector<pair<Kmer<k>, uint32_t>> hostCount;
		vector<pair<Kmer<k>, uint32_t>> gpuCount;
		hostCount.reserve(numKmers);
		gpuCount.reserve(numKmers);

		for (int ii = 0; ii < numBundles; ii++) {

			bundleSizes[ii] =
					(ii + 1) * kmersPerBundle < numKmers ?
							kmersPerBundle : numKmers % kmersPerBundle;

			//cout << "bundle " << ii << endl << "---------" << endl;
			for (int i = 0; i < bundleSizes[ii]; i++) {
				char dna[k];
				for (int j = 0; j < k; j++)
					dna[j] = proteins[rand() % 4];
				bundles[ii][i] = Kmer<k>(dna);
				int* ptr = (int*) &bundles[ii][i];
				//cout << i << " " << bundles[ii][i] << " " << *ptr << " "
				//		<< endl;
			}
			//cout << endl;
		}

		// cpu run
		auto t0 = std::chrono::high_resolution_clock::now();
		for (int ii = 0; ii < numBundles; ii++) {

			// insert kmer
			for (int i = 0; i < bundleSizes[ii]; i++) {
				mymap[bundles[ii][i]]++;
			}

		}
		auto t1 = std::chrono::high_resolution_clock::now();

		// get results
		for (auto it = mymap.begin(); it != mymap.end(); ++it) {
			hostCount.push_back(*it);
		}

		auto t2 = std::chrono::high_resolution_clock::now();
		auto hostTime1 = ::chrono::duration_cast < std::chrono::milliseconds
				> (t1 - t0);
		auto hostTime2 = ::chrono::duration_cast < std::chrono::milliseconds
				> (t2 - t1);
		auto hostTime = ::chrono::duration_cast < std::chrono::milliseconds
				> (t2 - t0);

		printf("gpu run\n");

		// gpu run
		t0 = std::chrono::high_resolution_clock::now();
		table.init(numKmers*2);

		for (int ii = 0; ii < numBundles; ii++) {
			// insert kmer
			table.insert(bundles[ii], bundleSizes[ii]);
		}

	//	std::cout << table << std::endl;

		t1 = std::chrono::high_resolution_clock::now();
		// get results
		table.extract(gpuCount);

		//std::cout << table << std::endl;

		t2 = std::chrono::high_resolution_clock::now();
		auto devTime1 = ::chrono::duration_cast < std::chrono::milliseconds
				> (t1 - t0);
		auto devTime2 = ::chrono::duration_cast < std::chrono::milliseconds
				> (t2 - t1);
		auto devTime = ::chrono::duration_cast < std::chrono::milliseconds
				> (t2 - t0);

		printf("              cpu      gpu\n");
		printf("------------------------------------\n");
		printf("phase 1: %6ims %6ims\n", (int) hostTime1.count(),
				(int) devTime1.count());
		printf("phase 2: %6ims %6ims\n", (int) hostTime2.count(),
				(int) devTime2.count());
		printf("------------------------------------\n");
		printf("         %6ims %6ims\tSpeedup=%.2f\n", (int) (hostTime.count()),
				(int) (devTime.count()),
				(hostTime.count() / (float) devTime.count()));

		// compare cpu with gpu results

		//sort by value using std::sort
		std::sort(hostCount.begin(), hostCount.end(), value_comparer<k>());
		std::sort(gpuCount.begin(), gpuCount.end(), value_comparer<k>());

		//sort by key using std::stable_sort
		std::stable_sort(hostCount.begin(), hostCount.end(), key_comparer<k>());
		std::stable_sort(gpuCount.begin(), gpuCount.end(), key_comparer<k>());

		//cerr << table << endl;

		if (hostCount.size() != gpuCount.size()) {
			cerr << "Error ! unequal list lengths: " << hostCount.size() << " "
					<< gpuCount.size() << endl;

			//cerr << table << endl;
			cerr << "gpuCount:" << endl;
			for (int i = 0; i < gpuCount.size(); i++) {
				cerr << i << ": " << gpuCount[i].first << ": "
						<< gpuCount[i].second << endl;
			}

			cerr << "hostCount:" << endl;
			for (int i = 0; i < hostCount.size(); i++) {
				cerr << i << ": " << hostCount[i].first << ": "
						<< hostCount[i].second << endl;
			}

			if (gpuCount.size() > hostCount.size()) {
				for (int i = 0; i < gpuCount.size(); i++) {

					if (i < gpuCount.size() - 1
							&& gpuCount[i].first == gpuCount[i + 1].first) {
						cerr << gpuCount[i].first << " counted twice!" << endl;
					}
					bool found = false;
					for (int j = 0; j < hostCount.size(); j++) {
						if (gpuCount[i].first == hostCount[j].first) {
							found = true;
							break;
						}
					}
					if (!found)
						cerr << gpuCount[i].first << " is missing in hostCount"
								<< endl;
				}
			} else {

				for (int i = 0; i < hostCount.size(); i++) {

					if (i < hostCount.size() - 1
							&& hostCount[i].first == hostCount[i + 1].first) {
						cerr << hostCount[i].first << " counted twice!" << endl;
					}
					bool found = false;
					for (int j = 0; j < gpuCount.size(); j++) {
						if (hostCount[i].first == gpuCount[j].first) {
							found = true;
							break;
						}
					}
					if (!found)
						cerr << hostCount[i].first << " is missing in gpuCount"
								<< endl;
				}

			}

		} else {

			cout << "unique kmers=" << hostCount.size() << endl;

			// output pairwise
			for (int i = 0; i < hostCount.size(); i++) {

				if (hostCount[i].first != gpuCount[i].first
						|| hostCount[i].second != gpuCount[i].second) {
					cerr << "Error: ";
					cerr << hostCount[i].first << ": " << hostCount[i].second
							<< " " << gpuCount[i].second << endl;
					return 0;
				}
			}
		}
	}

	for (int i = 0; i < numBundles; i++)
		delete bundles[i];
	delete[] bundles;
	delete[] bundleSizes;

	cudaDeviceReset();

	cout << "Success!" << endl;
	return 0;
}

int main() {

	return main2();
}

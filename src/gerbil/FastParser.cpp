/*
 * FastParser.cpp
 *
 *  Created on: 21.05.2015
 *      Author: marius
 */

#include "../../include/gerbil/FastParser.h"

void gerbil::Decompressor::setCompressedFastBundle(FastBundle* pCompressedFastBundle) {
	if (_compressedFastBundle)
		std::cerr << "ERROR" << __FILE__ << " " << __LINE__ << std::endl;
	else {
		_compressedFastBundle = pCompressedFastBundle;
		strm.avail_in = _compressedFastBundle->size;
		strm.next_in = (byte *) _compressedFastBundle->data;
		//printf("set_compressedFastBundle->size = %u\n", _compressedFastBundle->size);

		//printf("fastBundle.size = %u\n", fastBundle.size);
	}
}

uint xxxx = 0;

bool gerbil::Decompressor::decompress(uint tId) {
	if(!strm.avail_in)
		return _compressedFastBundle->size != 0;
	strm.avail_out = FAST_BUNDLE_DATA_SIZE_B - fastBundle.size;
	strm.next_out = ((byte *) fastBundle.data) + fastBundle.size;
	//printf("decompress_in: %u\n", strm.avail_in);
	//std::cout << strm.avail_in << " --> ";
	ret = inflate(&strm, Z_NO_FLUSH);
	//std::cout << "strm.avail_in = " << strm.avail_in << std::endl;
	//ret = inflateInit2(&strm, 16+MAX_WBITS);
	// ignore
	assert(ret != Z_STREAM_ERROR);  /* state not clobbered */
	switch (ret) {
		case Z_NEED_DICT:
			ret = Z_DATA_ERROR;     /* and fall through */
		case Z_DATA_ERROR:
		case Z_MEM_ERROR:
			(void) inflateEnd(&strm);
			//std::cerr << "ret[" << tId << "] = " << ret << std::endl;
			printf("ERROR: ret[%u] = %i\n", tId, ret);
			exit(66);
			//case Z_OK:
			//printf("OK\n");
	}

	fastBundle.size = FAST_BUNDLE_DATA_SIZE_B - strm.avail_out;
	//std::cout << "decompress_out: " << fastBundle.size << std::endl;
#if false
	xxxx++;
	if (xxxx >= 1179) {
		char s1 = fastBundle.data[fastBundle.size - 1];
		char s2 = fastBundle.data[500];
		fastBundle.data[fastBundle.size - 1] = '\0';
		fastBundle.data[500] = '\0';
		std::cout << fastBundle.data << "\n...\n" << (fastBundle.data + fastBundle.size - 500) << std::endl;
		fastBundle.data[fastBundle.size - 1] = s1;
		fastBundle.data[500] = s2;
		std::cout << "@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@" << std::endl;
		//exit(3);
	}
#endif
	return fastBundle.size < FAST_BUNDLE_DATA_SIZE_B;
}

void gerbil::FastParser::nextPart(char *&bp, char *&bp_end, const size_t& tId) {
	//printf("nextPart[%u]\n", tId);
	//IF_MESS_FASTPARSER(_sw[tId].hold();)
	FastBundle *curFastBundle = _curFastBundles[tId];
	bool finish = false;
	if(curFastBundle->compressType == fc_none || curFastBundle->size == 0) {
		curFastBundle->clear();
		finish = !_fastSyncSwapQueues[tId]->swapPop(curFastBundle);
		_curFastBundles[tId] = curFastBundle;
		if(!finish && curFastBundle->compressType != fc_none) {
			//printf("new bundle[%u]\n", tId);
			//printf("_decompressors[%u]->setCompressedFastBundle(%u)\n", tId, curFastBundle->size);
			_decompressors[tId]->reset();
			_decompressors[tId]->setCompressedFastBundle(curFastBundle);
		}
	}

	////////////////////
	if(finish) {
		//_curFastBundles[tId]->finalize(_curFastBundles[tId]->compressType);
		//printf("finish[%u]\n", tId);
		return;
	}

	if(curFastBundle->compressType == fc_none) {
		bp = curFastBundle->data;
		bp_end = curFastBundle->data + curFastBundle->size;
		//printf("curFastBundle->compressType == fc_none[%u]\n", tId);
	}
	else {
		if (curFastBundle->compressType == fc_DECOMPRESSOR) {
			_decompressors[tId]->fastBundle.clear();
			//printf("curFastBundle->compressType == fc_DECOMPRESSOR\n");
		}
		//printf("-->>>decompress[%u]\n", tId);
		bool isLast = false;
		while(_decompressors[tId]->decompress(tId)) {
			// aufraeumen
			curFastBundle = _decompressors[tId]->getCompressedFastBundle();
			curFastBundle->clear();
			if(_fastSyncSwapQueues[tId]->swapPop(curFastBundle)) {
				_decompressors[tId]->setCompressedFastBundle(curFastBundle);
			}
			else {
				//printf("last bundle[%u]\n", tId);
				//isLast = true;
				std::cerr << "ERROR " << __FILE__ << "  " << __LINE__ << std::endl;
				break;
			}
		}

		//printf("<<<--decompress[%u]\n", tId);
		if(curFastBundle->size == 0) {
			//printf("last bundle[%u]\n", tId);
			_curFastBundles[tId] = _decompressors[tId]->getCompressedFastBundle();
			_curFastBundles[tId]->compressType = fc_none;
			_decompressors[tId]->fastBundle.finalize(fc_DECOMPRESSOR);
		}
		else
			_curFastBundles[tId] = &(_decompressors[tId]->fastBundle);
		//printf("_curFastBundles[%u].size = %u\n", tId, _curFastBundles[tId]->size);
		bp = _decompressors[tId]->fastBundle.data;
		bp_end = bp + _decompressors[tId]->fastBundle.size;

		//_curFastBundles[tId]->finalize();
	}
	//IF_MESS_FASTPARSER(_sw[tId].proceed();)
}



inline void gerbil::FastParser::skipLineBreak(char *&bp, char *&bp_end, const size_t &tId) {
	if (!*(++bp))
		nextPart(bp, bp_end, tId);
	if (*bp == '\n' && !*(++bp))
		nextPart(bp, bp_end, tId);
}

inline void gerbil::FastParser::skipLine(char *&bp, char *&bp_end, const size_t &tId) {
	while (*bp != '\n' && *bp != '\r')
		if (!*(++bp))
			nextPart(bp, bp_end, tId);
	skipLineBreak(bp, bp_end, tId);
}

inline void gerbil::FastParser::skipLine(char *&bp, char *&bp_end, const size_t &l, const size_t &tId) {
	if (bp + l < bp_end)
		bp += l; //skip +\n
	else {
		size_t skip = l - (bp_end - bp);
		nextPart(bp, bp_end, tId);
		bp += skip;
	}
	skipLineBreak(bp, bp_end, tId);
}

inline void gerbil::FastParser::storeLine(
		char *&bp, char *&bp_end, size_t &l,
		ReadBundle *&readBundle, ReadBundle *&rbs, const size_t &tId, const char &skip
) {
	//std::cout << "storeLine" << std::endl;

	bool first = true;

	do {

		char *sl(bp);
		while (*bp != '\n' && *bp != '\r' && *bp) {
			++bp;
		}
		l = bp - sl;
		if (!*bp) {
			// store first part
			if (first) {
				if (!readBundle->add(l, sl)) {
					//readBundle->print();
					IF_MESS_FASTPARSER(_sw[tId].hold();)
					_syncQueue.swapPush(readBundle);
					IF_MESS_FASTPARSER(_sw[tId].proceed();)
					readBundle->add(l, sl);
				}
				}
				else if (!readBundle->expand(l, sl)) {
					//readBundle->print();
					readBundle->transfer(rbs);

					IF_MESS_FASTPARSER(_sw[tId].hold();)
					_syncQueue.swapPush(readBundle);
					IF_MESS_FASTPARSER(_sw[tId].proceed();)

					ReadBundle *swap = rbs;
					rbs = readBundle;
					readBundle = swap;
					readBundle->expand(l, sl);
				}
			nextPart(bp, bp_end, tId);
			// store second part
			sl = bp;
			while (*bp != '\n' && *bp != '\r') ++bp;
			l += bp - sl;
			if (!readBundle->expand(bp - sl, sl)) {
				//readBundle->print();
				readBundle->transfer(rbs);

				IF_MESS_FASTPARSER(_sw[tId].hold();)
				_syncQueue.swapPush(readBundle);
				IF_MESS_FASTPARSER(_sw[tId].proceed();)

				ReadBundle *swap = rbs;
				rbs = readBundle;
				readBundle = swap;
				readBundle->expand(bp - sl, sl);
			}
		} else {
			// store read
			if (first) {
				if (!readBundle->add(l, sl)) {
					//readBundle->print();
					IF_MESS_FASTPARSER(_sw[tId].hold();)
					_syncQueue.swapPush(readBundle);
					IF_MESS_FASTPARSER(_sw[tId].proceed();)
					readBundle->add(l, sl);
				}
			} else {
				if (!readBundle->expand(l, sl)) {
					//readBundle->print();
					readBundle->transfer(rbs);

					IF_MESS_FASTPARSER(_sw[tId].hold();)
					_syncQueue.swapPush(readBundle);
					IF_MESS_FASTPARSER(_sw[tId].proceed();)

					ReadBundle *swap = rbs;
					rbs = readBundle;
					readBundle = swap;
					readBundle->expand(l, sl);
				}
			}
		}
		skipLineBreak(bp, bp_end, tId);
		first = false;
	} while (*bp && skip && *bp != skip);
}


gerbil::FastParser::FastParser(
		uint32 &readBundlesNumber, TFileType fileType,
		TSeqType seqType, SyncSwapQueueSPSC<FastBundle> **fastSyncSwapQueues,
		const uint32_t &_readerParserThreadsNumber
) : _syncQueue(readBundlesNumber), _threadsNumber(_readerParserThreadsNumber), _processThreads(NULL) {
	_fileType = fileType;
	_seqType = seqType;
	_fastSyncSwapQueues = fastSyncSwapQueues;

	_readsNumber = 0;

	_curFastBundles = new FastBundle *[_threadsNumber];
	_decompressors  = new Decompressor *[_threadsNumber];
	_processThreads = new std::thread *[_threadsNumber];
	IF_MESS_FASTPARSER(
			_sw = new StopWatch[_threadsNumber];
			for (uint32_t i(0); i < _threadsNumber; ++i)
				_sw[i].setMode(CLOCK_THREAD_CPUTIME_ID);
	)
}

void gerbil::FastParser::processFastq(const size_t &tId) {
	//std::cout << "processFastq" << std::endl;
	char *bp;
	char *bp_end;
	size_t l;

	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();

	while (true) {
		nextPart(bp, bp_end, tId);
		//printf("see: _curFastBundles[%u]->size = %u\n", tId, _curFastBundles[tId]->size);
		if (!_curFastBundles[tId]->size)
			break;
		//printf("see 2\n");
		while (*bp) {
			//printf("see 3\n");
			// skip description
			skipLine(bp, bp_end, tId);

			// store read
			storeLine(bp, bp_end, l, readBundle, rbs, tId, '+');

			// skip + [description]
			skipLine(bp, bp_end, tId);

			if(xxxx == 1180)
				std::cout << bp << std::endl;
			// skip qualifiers
			do
				skipLine(bp, bp_end, l, tId);
			while (*bp && *bp != '@');

			++_readsNumber;
		}
	}
	if (!readBundle->isEmpty())
		_syncQueue.swapPush(readBundle);
	delete readBundle;
	delete rbs;
}

void gerbil::FastParser::processFasta(const size_t &tId) {
	char *bp;
	char *bp_end;
	size_t l;

	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();

	while (true) {
		nextPart(bp, bp_end, tId);
		if (!_curFastBundles[tId]->size)
			break;
		while (*bp) {
			// skip description
			skipLine(bp, bp_end, tId);

			storeLine(bp, bp_end, l, readBundle, rbs, tId, '>');

			++_readsNumber;
		}
	}
	if (!readBundle->isEmpty())
		_syncQueue.swapPush(readBundle);
	delete readBundle;
	delete rbs;
}

void gerbil::FastParser::processMultiline(const size_t &tId) {
	char *bp;
	char *bp_end;
	size_t l;

	ReadBundle *readBundle = new ReadBundle();
	ReadBundle *rbs = new ReadBundle();

	while (true) {
		nextPart(bp, bp_end, tId);
		if (!_curFastBundles[tId]->size)
			break;
		while (*bp) {

			storeLine(bp, bp_end, l, readBundle, rbs, tId, 0);

			++_readsNumber;
		}
	}
	if (!readBundle->isEmpty())
		_syncQueue.swapPush(readBundle);
	delete readBundle;
	delete rbs;
}


void gerbil::FastParser::process() {
	if (_seqType != st_reads) {
		std::cerr << "unsupported sequence type for fastq" << std::endl;
		exit(1);
	}
	for (uint32_t i = 0; i < _threadsNumber; i++) {
		_processThreads[i] = new std::thread([this](uint64_t tId) {
			IF_MESS_FASTPARSER(_sw[tId].start();)

			_curFastBundles[tId] = new FastBundle();
			_decompressors[tId] = new Decompressor();

			switch (_fileType) {
				case ft_fastq:
					processFastq(tId);
					break;
				case ft_fasta:
					processFasta(tId);
					break;
				case ft_multiline:
					processMultiline(tId);
					break;
				default:
					std::cerr << "unknown Filetype" << std::endl;
					exit(1);
			}

			delete _curFastBundles[tId];
			delete _decompressors[tId];

			IF_MESS_FASTPARSER(
					_sw[tId].stop();
					printf("time parser[%2lu]: %.3f s\n", tId, _sw[tId].get_s());
			)
		}, i);
	}
}

void gerbil::FastParser::join() {
	for (uint32_t i = 0; i < _threadsNumber; ++i) {
		_processThreads[i]->join();
		delete _processThreads[i];
	}
	_syncQueue.finalize();
	//printf("fastParser is rdy...\n");
}

gerbil::SyncSwapQueueMPMC<gerbil::ReadBundle> *gerbil::FastParser::getSyncQueue() {
	return &_syncQueue;
}

void gerbil::FastParser::print() {
	printf("number of reads        : %12lu\n", _readsNumber);
}

gerbil::FastParser::~FastParser() {
	delete[] _curFastBundles;
	delete[] _decompressors;
	delete[] _processThreads;
	IF_MESS_FASTPARSER(delete[] _sw;)
}


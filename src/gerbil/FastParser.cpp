/*
 * FastParser.cpp
 *
 *  Created on: 21.05.2015
 *      Author: marius
 */

#include "../../include/gerbil/FastParser.h"


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
	_processThreads = new std::thread *[_threadsNumber];
	IF_MESS_FASTPARSER(
			_sw = new StopWatch[_threadsNumber];
			for (uint32_t i(0); i < _threadsNumber; ++i)
				_sw[i].setMode(CLOCK_THREAD_CPUTIME_ID);
	)
}

void gerbil::FastParser::nextPart(char *&bp, char *&bp_end, const size_t &tId) {
	IF_MESS_FASTPARSER(_sw[tId].hold();)
	FastBundle *curFastBundle = _curFastBundles[tId];
	curFastBundle->clear();
	if (_fastSyncSwapQueues[tId]->swapPop(curFastBundle)) {
		bp = curFastBundle->data;
		bp_end = curFastBundle->data + curFastBundle->size;
	}
	_curFastBundles[tId] = curFastBundle;
	IF_MESS_FASTPARSER(_sw[tId].proceed();)
}

void gerbil::FastParser::processFastq(const size_t &tId) {
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

			// store read
			storeLine(bp, bp_end, l, readBundle, rbs, tId, '+');

			// skip + [description]
			skipLine(bp, bp_end, tId);

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
	delete[] _processThreads;
	IF_MESS_FASTPARSER(delete[] _sw;)
}


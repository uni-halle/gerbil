/*
 * SuperReader.cpp
 *
 *  Created on: 07.06.2015
 *      Author: marius
 */

#include "../../include/gerbil/SuperReader.h"

gerbil::SuperReader::SuperReader(const uint32 &superBundlesNumber,
		TempFile* tempFiles, const uint_tfn& tempFilesNumber) :
		_syncSwapQueue(superBundlesNumber) {
	_tempFiles = tempFiles;
	_tempFilesNumber = tempFilesNumber;
	_processThread = NULL;
	_tempFilesOrder = new uint_tfn[_tempFilesNumber];
}

gerbil::SuperReader::~SuperReader() {
	delete[] _tempFiles;
	delete[] _tempFilesOrder;
}

gerbil::uint_tfn* gerbil::SuperReader::getTempFilesOrder() {
	return _tempFilesOrder;
}

void gerbil::SuperReader::processThread() {
	IF_MESS_SUPERREADER(
			StopWatch sw;
			sw.start();
	)
	SuperBundle* sb = new SuperBundle;
	std::vector<std::pair<uint32_t, uint32_t>> binFileOrder;
	for (uint_tfn tempFileId(0); tempFileId < _tempFilesNumber; ++tempFileId)
		binFileOrder.push_back(
				std::make_pair(_tempFiles[tempFileId].getKMersNumber(),
						tempFileId));

	std::sort(binFileOrder.begin(), binFileOrder.end());

	std::sort(binFileOrder.begin(), binFileOrder.begin() + _tempFilesNumber / 4,
			std::greater<std::pair<uint32_t, uint32_t>>());

	uint i = 0;
	for (std::vector<std::pair<uint32_t, uint32_t>>::iterator it =
			binFileOrder.begin(); it != binFileOrder.end(); ++it)
		_tempFilesOrder[i++] = it->second;

	for (std::vector<std::pair<uint32_t, uint32_t>>::iterator it =
			binFileOrder.begin(); it != binFileOrder.end(); ++it) {
		uint_tfn tempFileId = it->second;

		if (!_tempFiles[tempFileId].isEmpty()) {
			_tempFiles[tempFileId].openR();
			while (_tempFiles[tempFileId].read(sb)) {
				sb->tempFileId = tempFileId;
				IF_MESS_SUPERREADER(sw.hold();)
				_syncSwapQueue.swapPush(sb);
			IF_MESS_SUPERREADER(sw.proceed();)
		}
		_tempFiles[tempFileId].close();
	}
	_tempFiles[tempFileId].remove();
}
_syncSwapQueue.finalize();

delete sb;

IF_MESS_SUPERREADER(
	sw.stop();
	printf("time SuperReader: %.3f s\n", sw.get_s());
)
}

void gerbil::SuperReader::process() {
_processThread = new std::thread([this] {
IF_DEB(printf("SuperReader start...\n"));
IF_DEB(DEB_startThreadClock());

processThread();

IF_DEB(DEB_stopThreadClock("SuperReader"));
});
}

void gerbil::SuperReader::join() {
_processThread->join();
delete _processThread;
}

gerbil::SyncSwapQueueSPMC<gerbil::SuperBundle>* gerbil::SuperReader::getSuperBundleQueue() {
return &_syncSwapQueue;
}

/*
 * KmcWriter.cpp
 *
 *  Created on: 07.06.2015
 *      Author: marius
 */

#include "../../include/gerbil/KmcWriter.h"


gerbil::KmcWriter::KmcWriter(std::string fileName, SyncSwapQueueMPSC<KmcBundle>* kmcSyncSwapQueue, const uint32_t &k) {
	_processThread = NULL;
	_fileName = fileName;
	_k = k;
	_kmcSyncSwapQueue = kmcSyncSwapQueue;
	std::remove(_fileName.c_str());
	_file = fopen (_fileName.c_str() , "wb" );
	if(!_file) {
		std::cerr << "unable to create output-file" << std::endl;
		exit(21);
	}
	setbuf(_file, NULL);

	_fileSize = 0;
}

gerbil::KmcWriter::~KmcWriter() {

}


void gerbil::KmcWriter::process() {
	if(_processThread)
		return;
	_processThread = new std::thread([this]{
		IF_MESS_KMCWRITER(
			StopWatch sw;
			sw.start();
		)
		IF_DEB(
			printf("kmcWriter start...\n");
		)
		KmcBundle* kb = new KmcBundle;

		uint64 s = 0;
		IF_MESS_KMCWRITER(sw.hold();)
		while(_kmcSyncSwapQueue->swapPop(kb)) {
			IF_MESS_KMCWRITER(sw.proceed();)
			//printf("++++++++++++++++++++++++++++++++++++++\n");
			if(!kb->isEmpty()) {
				_fileSize += kb->getSize();
				if(true) {
					//printf(">>>>write: %u\n", kb->getSize());
					fwrite ((char*) kb->getData() , 1 , kb->getSize() , _file );
				}
				else {
					const byte* p = kb->getData();
					const byte* end = p + kb->getSize();
					const byte* lp;
					char* a;
					uint32 l;
					while(p < end) {
						//fprintf(_file, "%3.0f|", (double)*(p++));
						//continue;
						l = (uint32)*(p++);
						if(l >= 255) {
							l = *((uint32*)p);
							p += 4;
						}
						a = getByteCodedSeq(p, _k);
						lp = p;
						p += getKMerCompactByteNumbers(_k);
						fprintf(_file, ">%9u\n%s\n", l, a);
						//fprintf(_file, "%u\n",l);
						delete[] a;
						//fprintf(_file, "[0x %02x %02x %02x %02x %02x %02x %02x] %s: %9d\n", lp[0], lp[1], lp[2], lp[3], lp[4], lp[5], lp[6], a, l);
					}
				}
			}
			kb->clear();
			//putchar('!');
			IF_MESS_KMCWRITER(sw.hold();)
		}
		IF_MESS_KMCWRITER(sw.proceed();)

		delete kb;
		if(_file)
			fclose(_file);
		IF_MESS_KMCWRITER(
			sw.stop();
			printf("kmcWriter: %7.3f s\n", sw.get_s());
		)
	});
}

void gerbil::KmcWriter::join() {
	_processThread->join();
	delete _processThread;
}

void gerbil::KmcWriter::print() {
	printf("size of output  : %12.3f MB\n", (double)_fileSize / 1024 / 1024);
}

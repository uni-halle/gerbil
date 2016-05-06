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

#ifndef FASTPARSER_H_
#define FASTPARSER_H_

#include "SyncQueue.h"
#include "Bundle.h"

namespace gerbil {

/*
 * extracts reads from each FastBundle and stores them in ReadBundles
 */
class FastParser {
private:
  TSeqType _seqType;                          // sequence type
  TFileType _fileType;                        // file type

  uint64 _readsNumber;                        // number of reads
  SyncSwapQueueMPMC<ReadBundle> _syncQueue;   // SyncSwapQueue for ReadBundles

  uint8 _threadsNumber;                       // number of threads
  std::thread** _processThreads;              // list of threads

  SyncSwapQueueSPSC<FastBundle>** _fastSyncSwapQueues;  // SyncSwapQueue for FastBundles

  FastBundle** _curFastBundles;                         // current FastBundles

  inline void skipLineBreak(char* &bp, char* &bp_end, const size_t &tId);
  inline void skipLine(char* &bp, char* &bp_end, const size_t &tId);
  inline void skipLine(char* &bp, char* &bp_end, const size_t &l, const size_t &tId);
  inline void storeLine(
    char* &bp, char* &bp_end, size_t &l,
    ReadBundle* &readBundle, ReadBundle* &rbs, const size_t &tId
  );

  void nextPart(char* &bp, char* &bp_end, const size_t &tId);

  void processFastq(const size_t &tId);
  void processFasta(const size_t &tId);
  void processMultiline(const size_t &tId);

  StopWatch* _sw;
public:
  SyncSwapQueueMPMC<ReadBundle>* getSyncQueue();          // returns SyncSwapQueue of ReadBundles

  inline uint64 getReadsNumber() { return _readsNumber; } // returns total number of reads

  /*
   * constructor
   */
  FastParser(
      uint32 &readBundlesNumber, TFileType fileType, TSeqType seqType,
      SyncSwapQueueSPSC<FastBundle>** _fastSyncSwapQueues,
      const uint8 &_readerParserThreadsNumber
  );

  /*
   * starts the entire working process
   */
  void process();

  /*
   * joins all threads
   */
  void join();

  /*
   * prints some statistical outputs
   */
  void print();

  /*
   * destructor
   */
  ~FastParser();
};

}


#endif /* FASTPARSER_H_ */

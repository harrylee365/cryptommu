#ifndef DMACONTROLLER_H
#define DMACONTROLLER_H

#include <stdint.h>
#include <map>
#include <vector>
#include <list>
#include <string>
#include "SPMInterface.hh"
#include "NetworkInterface.hh"
#include "modules/Common/BaseCallbacks.hh"
#include "modules/linked-prefetch-tile/prefetcher-tile.hh"
#include "modules/MsgLogger/MsgLogger.hh"

namespace LCAcc
{

class TLBEntry
{
public:
  uint64_t vpBase;
  uint64_t ppBase;
  bool free;
  uint64_t mruTick;
  TLBEntry() : vpBase(0), ppBase(0), free(true), mruTick(0) {}
  void setMRU()
  {
    mruTick = GetSystemTime();
  };
};

class BaseTLBMemory
{
public:
  virtual bool lookup(uint64_t vp_base, uint64_t& pp_base, bool set_mru = true) = 0;
  virtual void insert(uint64_t vp_base, uint64_t pp_base) = 0;
  virtual void flushAll() = 0;
};
// assoc mean associativity
// sets mean set
class TLBMemory : public BaseTLBMemory
{
  int numEntries;
  int assoc;
  int sets;
  TLBEntry **entries;

protected:
  TLBMemory() {}

public:
  TLBMemory(int _numEntries, int associativity) :
    numEntries(_numEntries), assoc(associativity)
  {
    if (assoc == 0) {
      assoc = numEntries;
    }

    assert(numEntries % assoc == 0);
    sets = numEntries / assoc;
    //std::cout << "class: sets" << sets << std::endl;
    //std::cout << "class: assoc" << assoc << std::endl;
    entries = new TLBEntry*[sets];

    for (int i = 0; i < sets; i++) {
      entries[i] = new TLBEntry[assoc];
    }
  }
  virtual ~TLBMemory()
  {
    for (int i = 0; i < assoc; i++) {
      delete [] entries[i];
    }

    delete [] entries;
  }

  virtual bool lookup(uint64_t vp_base, uint64_t& pp_base, bool set_mru = true);
  virtual void insert(uint64_t vp_base, uint64_t pp_base);
  virtual void flushAll();
};

class InfiniteTLBMemory : public BaseTLBMemory
{
  std::map<uint64_t, uint64_t> entries;
public:
  InfiniteTLBMemory() {}
  ~InfiniteTLBMemory() {}

  bool lookup(uint64_t vp_base, uint64_t& pp_base, bool set_mru = true)
  {
    auto it = entries.find(vp_base);
    if (it != entries.end()) {
      pp_base = it->second;
      return true;
    } else {
      pp_base = 0;
      return false;
    }
  }
  void insert(uint64_t vp_base, uint64_t pp_base)
  {
    entries[vp_base] = pp_base;
  }
  void flushAll() {}
};

class DMAController
{
public:
  class AccessType
  {
  public:
    static const int Read;
    static const int Write;
    static const int ReadLock;
    static const int WriteLock;
    static const int WriteUnlock;
    static const int Unlock;
  };
private:
  class SignalEntry
  {
  public:
    int ID;
    uint64_t dstAddr;
    uint64_t requestedAddr;
    unsigned int size;
    CallbackBase* onFinish;
    bool Matches(const SignalEntry& s)
    {
      return s.ID == ID && s.dstAddr == dstAddr && s.requestedAddr == requestedAddr && s.size == size;
    }
  };
  uint64_t GetPageAddr(uint64_t addr) const;
  std::vector<SignalEntry> remoteSignals;
  std::vector<SignalEntry> localSignals;
  std::map<uint64_t, std::list<TransferData*> > MSHRs;
  std::map<uint64_t, std::list<TransferData*> > MAC_queue;
  DMAEngineHandle *dmaDevice;
  prftch_direct_interface_t* dmaInterface;
  SPMInterface* spm;
  NetworkInterface* network;
  Arg1CallbackBase<uint64_t>* onAccViolation;
  Arg3CallbackBase<uint64_t, uint64_t, uint64_t>* onTLBMiss;
  Arg1CallbackBase<uint64_t>* onMACver;
  // std::map<uint64_t, uint64_t> tlbMap;
  int buffer;
  void OnAccessError(uint64_t logicalAddr);
  void OnNetworkMsg(int src, const void* buffer, unsigned int bufSize);
  typedef Arg1MemberCallback<DMAController, uint64_t, &DMAController::OnAccessError> OnAccessErrorCB;
  typedef Arg3MemberCallback<DMAController, int, const void*, unsigned int, &DMAController::OnNetworkMsg> OnNetworkMsgCB;
  bool isHookedToMemory;
public:
  DMAController(NetworkInterface* ni, SPMInterface* spmInterface, Arg3CallbackBase<uint64_t, uint64_t, uint64_t>* TLBMiss, Arg1CallbackBase<uint64_t>* accessViolation, Arg1CallbackBase<uint64_t>* MACver);
  ~DMAController();
  void BeginTransfer(int srcSpm, uint64_t srcAddr, const std::vector<unsigned int>& srcSize, const std::vector<int>& srcStride, int dstSpm, uint64_t dstAddr, const std::vector<unsigned int>& dstSize, const std::vector<int>& dstStride, size_t elementSize, CallbackBase* finishedCB);
  void BeginTransfer(int srcSpm, uint64_t srcAddr, const std::vector<unsigned int>& srcSize, const std::vector<int>& srcStride, int dstSpm, uint64_t dstAddr, const std::vector<unsigned int>& dstSize, const std::vector<int>& dstStride, size_t elementSize, int priority, CallbackBase* finishedCB);
  void PrefetchMemory(uint64_t baseAddr, const std::vector<unsigned int>& size, const std::vector<int>& stride, size_t elementSize);
  void BeginSingleElementTransfer(int mySPM, uint64_t src, uint64_t dst, uint32_t size, int type, CallbackBase* finishedCB);
  void SetBuffer(int buf);
  void FlushTLB();
  void AddTLBEntry(uint64_t vAddr, uint64_t pAddr);
  void HookToMemoryController(const std::string& deviceName);

  enum TransferStatus {
    Running,
    TlbWait
  };

  typedef MemberCallback2<DMAController, uint64_t, uint64_t, &DMAController::AddTLBEntry> AddTLBEntryCB;

protected:

  inline std::string GetDeviceName()
  {
    char s[100];
    sprintf(s, "dma%02d", network->GetNetworkPort());
    return s;
  }

  // Access Stats
  uint64_t hits;
  uint64_t mshrhits;
  uint64_t misses;
  uint64_t flushTlb;
  uint64_t tlbCycles;
  uint64_t BCCMshrhits;
  uint64_t* protection_table_Memory;

  TransferStatus transferStatus;
  uint64_t timeStamp;

public:
  // private TLB entries
  int numEntries;

  int associativity;

  int hitLatency;

  int AESLatency = 40;

  int acc_count=0;

  int MacLatency = AESLatency + hitLatency;

  BaseTLBMemory *tlbMemory;

  void beginTranslateTiming(TransferData* td);

  void translateTiming(TransferData* td);

  void finishTranslation(uint64_t vp_base, uint64_t pp_base, uint64_t MAC);
   void finishVerification(uint64_t vp_base);
  void flushAll();

  typedef Arg1MemberCallback<DMAController, TransferData*, &DMAController::beginTranslateTiming> beginTranslateTimingCB;
  typedef MemberCallback1<DMAController, TransferData*, &DMAController::translateTiming> translateTimingCB;

  uint64_t getTlbHits()
  {
    return hits + mshrhits;
  }
  uint64_t getTlbMisses()
  {
    return misses;
  }
  uint64_t getTlbFlush()
  {
    return flushTlb;
  }
  uint64_t getTlbCycles()
  {
    return tlbCycles;
  }
   uint64_t getBCCMshrhits()
  {
    return BCCMshrhits;
  }
};

}

#endif

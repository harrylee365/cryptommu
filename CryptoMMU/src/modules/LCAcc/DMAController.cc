#include <string>
#include <cassert>
#include <list>

#include "DMAController.hh"
#include "LCAccCommandListing.hh"
#include "SimicsInterface.hh"
#include "../MsgLogger/MsgLogger.hh"
#include "arch/isa_traits.hh"
#include "../Common/BitConverter.hh"
#include "../MsgLogger/MsgLogger.hh"
#include "mem/ruby/common/Global.hh"
#include "../TLBHack/TLBHack.hh"
#include "arch/vtophys.hh"

using namespace LCAcc;

const int opmode_read = 1;
const int opmode_write = 2;
const int opmode_send = 3;
const int opmode_recv = 4;

namespace LCAcc
{

const int DMAController::AccessType::Read = 0;
const int DMAController::AccessType::Write = 1;
const int DMAController::AccessType::ReadLock = 3;
const int DMAController::AccessType::WriteLock = 4;
const int DMAController::AccessType::WriteUnlock = 5;
const int DMAController::AccessType::Unlock = 6;

DMAController::DMAController(NetworkInterface* ni, SPMInterface* spm,
                             Arg3CallbackBase<uint64_t, uint64_t, uint64_t>* TLBMiss,
                             Arg1CallbackBase<uint64_t>* accessViolation,Arg1CallbackBase<uint64_t>* MACver)
{
  //values on rhs will be different based on new objects which will be assigned to class variable
  dmaDevice = g_dmaDevice[RubySystem::deviceIDtoAccID(spm->GetID())];
  //std::cout << "DEVICE ID" << RubySystem::deviceIDtoAccID(spm->GetID()) << std::endl;
  //std::cout << "ACC ID" << spm->GetID() << std::endl;
  dmaInterface = g_dmaInterface;
  network = ni;
  onAccViolation = accessViolation;
  onTLBMiss = TLBMiss;
  onMACver = MACver;
  buffer = -1;
  this->spm = spm;
  //std::cout <<"This DMA controller" << this << std::endl;
  dmaInterface->Configure(dmaDevice, ni->GetNodeID(), spm->GetID(),
                          beginTranslateTimingCB::Create(this),
                          OnAccessErrorCB::Create(this));

  ni->RegisterRecvHandler(OnNetworkMsgCB::Create(this));

  isHookedToMemory = false;

  numEntries = RubySystem::getLCAccTLBSize();
  hitLatency = RubySystem::getLCAccTLBLatency();
  associativity = RubySystem::getLCAccTLBAssoc();
  //std::cout <<"input TLB Size" << numEntries << std::endl;
  //std::cout <<"input TLB associativity" << associativity << std::endl;
  if (numEntries > 0) {
    tlbMemory = new TLBMemory(numEntries, associativity);
    //std::cout<<"TLB memory" << tlbMemory << std::endl;
  } else {
    tlbMemory = new InfiniteTLBMemory();
  }

  hits = 0;
  mshrhits = 0;
  misses = 0;
  flushTlb = 0;
  transferStatus = Running;
  timeStamp = 0;
  BCCMshrhits = 0;
  //protection_table_Memory =new int[1024*1024];
  
}

DMAController::~DMAController()
{
  if (isHookedToMemory) {
    assert(dmaDevice);
    assert(dmaInterface);
    dmaInterface->UnhookMemoryPort(dmaDevice);
  }
}

uint64_t
DMAController::GetPageAddr(uint64_t addr) const
{
  return addr - (addr % TheISA::PageBytes);
}

void
DMAController::finishTranslation(uint64_t vp_base, uint64_t pp_base, uint64_t MAC_return)
{
 
  std::list<TransferData*> &tds = MSHRs[vp_base];
  std::list<TransferData*>::iterator it;
  for (it = tds.begin(); it != tds.end(); it++) {
    TransferData* td = (*it);
    if (MAC_return==1 )
    {
         acc_count++;
         uint64_t offset = td->getVaddr() % TheISA::PageBytes;
         td->setPaddr(pp_base+offset);
         dmaInterface->finishTranslation(dmaDevice, td); 
    }

    else if (MAC_return==0)
    {
    uint64_t offset = td->getVaddr() % TheISA::PageBytes;
    td->setPaddr(pp_base + offset);
    dmaInterface->finishTranslation(dmaDevice, td);
    }
    
  }
     
    tlbMemory->insert(vp_base, pp_base);
    MSHRs.erase(vp_base);
   
    if (MSHRs.empty()) {
    transferStatus = Running;
    tlbCycles += GetSystemTime() - timeStamp;

  }

}

void
DMAController::beginTranslateTiming(TransferData* td)
{
  if (RubySystem::idealMMU()) {
    // The below code implements a perfect TLB with instant access to the
    // host page table.
    uint64_t vaddr = td->getVaddr();
    uint64_t vp_base = GetPageAddr(vaddr);
    uint64_t offset = vaddr - vp_base;

    if (g_TLBHack_interface->PageKnown(1, vp_base)) {
      uint64_t pp_base = g_TLBHack_interface->Lookup(1, vp_base);
      td->setPaddr(pp_base + offset);
      dmaInterface->finishTranslation(dmaDevice, td);
    } else {
      System *m5_system = *(System::systemList.begin());
      ThreadContext* cpu = m5_system->getThreadContext(0);
      uint64_t pp_base = TheISA::vtophys(cpu, vp_base);
      td->setPaddr(pp_base + offset);
      dmaInterface->finishTranslation(dmaDevice, td);
    }
  } else {

    SimicsInterface::RegisterCallback(
      translateTimingCB::Create(this, td), hitLatency);
  }
}

void
DMAController::translateTiming(TransferData* td)
{
 
  uint64_t vaddr = td->getVaddr();
  uint64_t vp_base = GetPageAddr(vaddr);
  uint64_t offset = vaddr - vp_base;
  uint64_t pp_base;
  uint64_t MAC_dma;

  
  if (tlbMemory->lookup(vp_base, pp_base)) {
    hits++;
    //BCCMshrhits++;
    //std::cout <<"No Read Acceleration is running" << std::endl;
    td->setPaddr(pp_base + offset);
    td->MAC_ver =1; //if it is a read request immediately finish translation
    MSHRs[vp_base].push_back(td);
    MAC_dma =1;
    MAC_verfication++;
    //std::cout <<"MAC_verifcation" << hits << std::endl;
    onTLBMiss->Call(vp_base, MAC_dma,pp_base);
  }
  else
  {
    misses++;
    //std::cout <<"Miss in LcAcc TLB" << std::endl;
    if (MSHRs.empty()) {
      transferStatus = TlbWait;
      timeStamp = GetSystemTime();
    }
    td->MAC_ver =0;
    MAC_dma =0;
    pp_base=0;
    MSHRs[vp_base].push_back(td);
    onTLBMiss->Call(vp_base, MAC_dma,pp_base);

  }
   
}

void
DMAController::OnAccessError(uint64_t logicalAddr)
{
  onAccViolation->Call(logicalAddr);
}

void
DMAController::OnNetworkMsg(int src, const void* buf, unsigned int bufSize)
{
  assert(bufSize >= 4);
  assert(buf);
  uint32_t* b = (uint32_t*)buf;
  uint32_t opcode = b[0];

  if (opcode == DMA_MEMORY_REQUEST) {
    assert(bufSize == sizeof(uint32_t) * 6);
    SignalEntry target;
    target.ID = src;
    BitConverter bc_dst;
    bc_dst.u32[0] = b[1];
    bc_dst.u32[1] = b[2];
    target.dstAddr = bc_dst.u64[0];
    BitConverter bc_request;
    bc_request.u32[0] = b[3];
    bc_request.u32[1] = b[4];
    target.requestedAddr = bc_request.u64[0];
    target.size = b[5];
    target.onFinish = NULL;

    for (std::vector<SignalEntry>::iterator i = localSignals.begin(); i != localSignals.end(); i++) {
      if (i->Matches(target)) {
        std::vector<uint32_t> msg;
        msg.push_back(DMA_MEMORY_RESPONSE);
        bc_request.u64[0] = i->requestedAddr;
        msg.push_back(bc_request.u32[0]);
        msg.push_back(bc_request.u32[1]);
        bc_dst.u64[0] = i->dstAddr;
        msg.push_back(bc_dst.u32[0]);
        msg.push_back(bc_dst.u32[1]);
        msg.push_back(i->size);

        for (unsigned int j = 0; j < i->size / sizeof(uint32_t) + ((i->size % sizeof(uint32_t)) ? 1 : 0); j++) {
          msg.push_back(0);
        }

        spm->Read(i->requestedAddr, i->size, &(msg[6]));
        network->SendMessage(i->ID, &(msg[0]), msg.size() * sizeof(int32_t));
        assert(i->onFinish);
        SimicsInterface::RegisterCallback(i->onFinish, 0);
        localSignals.erase(i);
        return;
      }
    }

    remoteSignals.push_back(target);
  } else if (opcode == DMA_MEMORY_RESPONSE) {
    assert(bufSize >= (int)sizeof(int32_t) * 6);
    uint32_t* args = (uint32_t*) & (b[1]);
    SignalEntry r;
    r.ID = src;
    BitConverter bc_dst;
    bc_dst.u32[0] = args[2];
    bc_dst.u32[1] = args[3];
    r.dstAddr = bc_dst.u64[0];
    BitConverter bc_request;
    bc_request.u32[0] = args[0];
    bc_request.u32[1] = args[1];
    r.requestedAddr = bc_request.u64[0];
    r.size = args[4];
    assert(bufSize == sizeof(uint32_t) * 6 + r.size);

    // bool resolved = false;
    for (std::vector<SignalEntry>::iterator i = localSignals.begin(); i != localSignals.end(); i++) {
      if (i->Matches(r)) {
        // resolved = true;
        assert(i->onFinish);
        SimicsInterface::RegisterCallback(i->onFinish, 0);
        spm->Write(r.dstAddr, r.size, &(args[5]));
        localSignals.erase(i);
        return;
      }
    }

    //std::cerr << "DMA Memory Response came unsolicited" << std::endl;
    assert(0);
  }
}

void
DMAController::BeginTransfer(int srcSpm, uint64_t srcAddr,
                             const std::vector<unsigned int>& srcSize,
                             const std::vector<int>& srcStride, int dstSpm, uint64_t dstAddr,
                             const std::vector<unsigned int>& dstSize,
                             const std::vector<int>& dstStride, size_t elementSize,
                             CallbackBase* finishCB)
{
  BeginTransfer(srcSpm, srcAddr, srcSize, srcStride,
                dstSpm, dstAddr, dstSize, dstStride, elementSize, 0, finishCB);
}

void
DMAController::BeginTransfer(int srcSpm, uint64_t srcAddr,
                             const std::vector<unsigned int>& srcSize,
                             const std::vector<int>& srcStride, int dstSpm, uint64_t dstAddr,
                             const std::vector<unsigned int>& dstSize,
                             const std::vector<int>& dstStride, size_t elementSize, int priority,
                             CallbackBase* finishCB)
{
  assert(finishCB);
  assert(srcSize.size() == srcStride.size());
  assert(dstSize.size() == dstStride.size());
  assert(srcSize.size() > 0);
  assert(dstSize.size() > 0);
  assert(elementSize <= 8);
  assert(elementSize > 0);
  assert(srcSpm != -1 || dstSpm != -1);
  assert(srcSpm == spm->GetID() || dstSpm == spm->GetID());

  if (srcSpm == -1 || dstSpm == -1) {
    dmaInterface->StartTransferPrio(dmaDevice,
                                    srcSpm, srcAddr, srcSize.size(), &(srcSize[0]), &(srcStride[0]),
                                    dstSpm, dstAddr, dstSize.size(), &(dstSize[0]), &(dstStride[0]),
                                    elementSize, priority, buffer, finishCB);
  } else if (dstSpm == spm->GetID()) {
    uint32_t msg[6];
    msg[0] = DMA_MEMORY_REQUEST;
    BitConverter bc_dst;
    bc_dst.u64[0] = dstAddr;
    msg[1] = bc_dst.u32[0];
    msg[2] = bc_dst.u32[1];
    BitConverter bc_src;
    bc_src.u64[0] = srcAddr;
    msg[3] = bc_src.u32[0];
    msg[4] = bc_src.u32[1];
    uint32_t size = elementSize;

    for (size_t i = 0; i < dstSize.size(); i++) {
      size *= dstSize[i];
    }

    msg[5] = size;
    network->SendMessage(srcSpm, msg, 6 * sizeof(uint32_t));
    SignalEntry target;
    target.ID = srcSpm;
    target.dstAddr = dstAddr;
    target.requestedAddr = srcAddr;
    target.size = elementSize;

    for (size_t i = 0; i < srcSize.size(); i++) {
      target.size *= srcSize[i];
    }

    target.onFinish = finishCB;
    localSignals.push_back(target);
  } else {
    assert(srcSpm == spm->GetID());
    SignalEntry target;
    target.ID = dstSpm;
    target.dstAddr = dstAddr;
    target.requestedAddr = srcAddr;
    target.size = elementSize;

    for (size_t i = 0; i < srcSize.size(); i++) {
      target.size *= srcSize[i];
    }

    target.onFinish = finishCB;

    for (std::vector<SignalEntry>::iterator i = remoteSignals.begin(); i != remoteSignals.end(); i++) {
      if (i->Matches(target)) {
        std::vector<uint32_t> msg;
        msg.push_back(DMA_MEMORY_RESPONSE);
        BitConverter bc_request;
        bc_request.u64[0] = i->requestedAddr;
        msg.push_back(bc_request.u32[0]);
        msg.push_back(bc_request.u32[1]);
        BitConverter bc_dst;
        bc_dst.u64[0] = i->dstAddr;
        msg.push_back(bc_dst.u32[0]);
        msg.push_back(bc_dst.u32[1]);
        msg.push_back(i->size);

        for (unsigned int j = 0; j < i->size / sizeof(uint32_t) + ((i->size % sizeof(uint32_t)) ? 1 : 0); j++) {
          msg.push_back(0);
        }

        spm->Read(i->requestedAddr, i->size, &(msg[6]));
        network->SendMessage(i->ID, &(msg[0]), msg.size() * sizeof(int32_t));
        SimicsInterface::RegisterCallback(target.onFinish, 0);
        assert(i->onFinish == NULL);
        remoteSignals.erase(i);
        return;
      }
    }

    localSignals.push_back(target);
  }
}

void
DMAController::PrefetchMemory(uint64_t baseAddr,
                              const std::vector<unsigned int>& size,
                              const std::vector<int>& stride, size_t elementSize)
{
  dmaInterface->Prefetch(dmaDevice, baseAddr, size.size(),
                         &(size[0]), &(stride[0]), elementSize);
}

void
DMAController::SetBuffer(int buf)
{
  buffer = buf;
}

void
DMAController::FlushTLB()
{
  flushTlb++;
  tlbMemory->flushAll();
}

void
DMAController::AddTLBEntry(uint64_t vAddr, uint64_t pAddr)
{
  uint64_t vp_base = GetPageAddr(vAddr);
  uint64_t pp_base = GetPageAddr(pAddr);

  assert(vAddr - vp_base == pAddr - pp_base);

  tlbMemory->insert(vp_base, pp_base);
}

void
DMAController::HookToMemoryController(const std::string& deviceName)
{
  assert(dmaInterface);
  assert(dmaDevice);

  if (isHookedToMemory) {
    dmaInterface->UnhookMemoryPort(dmaDevice);
  }

  dmaInterface->HookToMemoryPort(dmaDevice, deviceName.c_str());
}

void
DMAController::BeginSingleElementTransfer(int mySPM, uint64_t spmAddr,
    uint64_t memAddr, uint32_t size, int type, CallbackBase* finishedCB)
{
  //std::cout << "Accessing single element " << spmAddr << std::endl;
  bool isRead = (type == AccessType::Read || type == AccessType::ReadLock);
  bool isLock = (type == AccessType::ReadLock || type == AccessType::WriteLock);
  bool isUnlock = (type == AccessType::WriteUnlock || type == AccessType::Unlock);
  assert(!isLock || !isUnlock);
  assert(!isLock);//for now, since locking is disabled
  assert(!isUnlock);//for now, since locking is disabled
  //std::cout << "isRead" << isRead << std::endl;
  //std::cout << "memAddr" << pAddr << std::endl;
  dmaInterface->StartSingleTransferPrio(dmaDevice, isRead ? -1 : mySPM,
                                        isRead ? memAddr : spmAddr, isRead ? mySPM : -1,
                                        isRead ? spmAddr : memAddr, size, 0, buffer, finishedCB);
}

void
TLBMemory::flushAll()
{
  for (int set = 0; set < sets; set++) {
    for (int assc = 0; assc < assoc; assc++) {
      //std::cout<< "I am breaking in Flush" << std::endl;
      //std::cout << "sets" << set   << "assc" << assc << std::endl;
      //std::cout << "sets" << sets   << "assoc" << assoc << std::endl;
        if(entries[set][assc].free == false)
              {
                  entries[set][assc].free = true;
                  entries[set][assc].vpBase = 0;
                  entries[set][assc].ppBase = 0;
              }
               
      
    }
  }
}

bool
TLBMemory::lookup(uint64_t vp_base, uint64_t& pp_base, bool set_mru)
{
  int set = (vp_base / TheISA::PageBytes) % sets;
  //std::cout << "set" << set << std::endl;
  for (int i = 0; i < assoc; i++) {
    if (entries[set][i].vpBase == vp_base && !entries[set][i].free) {
      pp_base = entries[set][i].ppBase;
      assert(entries[set][i].mruTick > 0);

      if (set_mru) {
        entries[set][i].setMRU();
      }

      return true;
    }
  }

  pp_base = 0;
  return false;
}

void
TLBMemory::insert(uint64_t vp_base, uint64_t pp_base)
{
  uint64_t a;

  if (lookup(vp_base, a)) {
    return;
  }

  int set = (vp_base / TheISA::PageBytes) % sets;
  TLBEntry* entry = NULL;
  uint64_t minTick = GetSystemTime();

  for (int i = 0; i < assoc; i++) {
    if (entries[set][i].free) {
      entry = &entries[set][i];
      break;
    } else if (entries[set][i].mruTick < minTick) {
      minTick = entries[set][i].mruTick;
      entry = &entries[set][i];
    }
  }

  assert(entry);

  entry->vpBase = vp_base;
  entry->ppBase = pp_base;
  entry->free = false;
  entry->setMRU();
}

}

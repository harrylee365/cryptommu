#include "LCAccDevice.hh"
#include "SimicsInterface.hh"
#include "DMAController.hh"
#include "SPMInterface.hh"
#include "NetworkInterface.hh"
#include "memInterface.hh"
#include <iostream>
#include <fstream>
#include <cassert>
#include <list>
#include "LCAccCommandListing.hh"
#include "LCAccOperatingModeInclude.hh"
#include "../MsgLogger/MsgLogger.hh"
#include "../Common/TransferDescription.hh"
#include "../Common/ComputeDescription.hh"
#include "../Common/BitConverter.hh"
#include "../../mem/ruby/system/SimicsNetworkPortInterface.hh"
#include "../NetworkInterrupt/NetworkInterrupts.hh"
#include "memInterface.hh"
using namespace LCAcc;

#define NO_SPM_ID -1
#define NO_DEVICE_ID ""

namespace LCAcc
{
bool LCAccDevice::initializedModeFactory = false;
ObjFactory<std::string, LCAccOperatingMode> LCAccDevice::modeFactory;
std::map<std::string, int> LCAccDevice::nameToMode;
std::map<int, std::string> LCAccDevice::modeToName;
std::map<int, LCAccOperatingMode*> LCAccDevice::modeSampleSet;
int PrefetcherID[4] = {1000, 1001, 1002, 1003};

void
LCAccDevice::PrefetchOrder(const TransferOrder& to)
{
  if (to.src.spm != NO_SPM_ID && to.dst.spm != NO_SPM_ID) {
    return;
  }

  const EndNodeTransferData& entd =
    (to.src.spm == NO_SPM_ID) ? to.src : to.dst;
  dma->PrefetchMemory(entd.addr, entd.sizeSet, entd.strideSet,
                      to.elementSize);
}

void
LCAccDevice::EmitActionReads(Action* a)
{
  assert(a);
  assert(a->computeReady == false);
  assert(a->retireReady == false);
  assert(a->readsInFlight == 0 && a->writesInFlight == 0);

  for (size_t i = 0; i < a->prefetchSet.size(); i++) {
    PrefetchOrder(a->prefetchSet[i]);
  }

  if (a->readSet.size() > 0) {
    // dma read required
    // ML_LOG(GetDeviceName(), "Begin DMA read task " << a->taskID);
    pendingReadSet.push(a);

    a->readsInFlight = a->readSet.size();

    for (size_t i = 0; i < a->readSet.size(); i++) {
      TransferOrder& t = a->readSet[i];
      assert(t.src.IsSPM() || t.dst.IsSPM());

      if (t.src.IsSPM()) {
        //ML_LOG(GetDeviceName(), "BEGIN Transfer Chain Read "
        //       << & (a->readSet[i]) << ' ' << t.src.spm << " -> " << t.dst.spm);
      } else {
        //ML_LOG(GetDeviceName(), "BEGIN Transfer Read " << i
        //       << " Task " << a->taskID);
      }

      dma->BeginTransfer(t.src.spm, (uint64_t)t.src.addr,
                         t.src.sizeSet, t.src.strideSet, t.dst.spm,
                         (uint64_t)t.dst.addr, t.dst.sizeSet, t.dst.strideSet,
                         t.elementSize, Priority::Medium,
                         FinishedActionReadCB::Create(this, a, i));
    }
  } else {
    // no read required, begin compute
    EmitActionCompute(a);
  }
}

void
LCAccDevice::FinishedActionRead(Action* a, size_t index)
{
  assert(a);
  assert(a->readsInFlight > 0);
  assert(a->readSet.size() > index);
  assert(a->computeReady == false);

  a->readsInFlight--;
  TransferOrder& t = a->readSet[index];

  if (t.src.IsSPM()) {
    //ML_LOG(GetDeviceName(), "END Transfer Chain Read "
    //       << & (a->readSet[index]) << " "  << t.src.spm << " -> " << t.dst.spm);
  } else {
    //ML_LOG(GetDeviceName(), "END Transfer Read " << index
    //       << " Task " << a->taskID);
  }

  if (a->readsInFlight == 0) {
    // done reading, begin compute
    // ML_LOG(GetDeviceName(), "End DMA read task " << a->taskID);
    pendingReadSet.pop();
    StartTask();
    EmitActionCompute(a);
  }
}

void
LCAccDevice::EmitActionCompute(Action* a)
{
  assert(a);
  assert(a->readsInFlight == 0 && a->writesInFlight == 0);
  assert(a->computeReady == false);

  a->computeReady = true;
  assert(!awaitingComputeActionSet.empty());

  while (!awaitingComputeActionSet.empty() &&
         awaitingComputeActionSet.front()->computeReady) {

    Action* ca = awaitingComputeActionSet.front();
    assert(ca);
    awaitingComputeActionSet.pop();
    std::vector<uint64_t> readBufferAddrs;
    std::vector<uint64_t> seedAddrs;

    assert(ca->compute);
    ca->compute->pendingComputeElements = ca->compute->computeCount;

    for (size_t i = 0; i < ca->compute->computeCount; i++) {
      seedAddrs.clear();
      readBufferAddrs.clear();

      for (size_t x = 0; x < ca->compute->argumentAddressGen.size(); x++) {
        seedAddrs.push_back(ca->compute->CalcArgAddress(x, i));
      }

      ca->compute->mode->GetSPMReadIndexSet(i,
                                            ca->compute->maxCompute, ca->taskID, seedAddrs,
                                            ca->compute->argumentActive, readBufferAddrs);
      ComputeElement* ce = ca->compute->NextComputeElement();
      ce->readsRemaining = readBufferAddrs.size();

      if (readBufferAddrs.size() > 0) {
        // spm read required
        for (size_t x = 0; x < readBufferAddrs.size(); x++) {
          spm->TimedRead(readBufferAddrs[x], 1,
                         SPMReadCompleteCB::Create(this, ce));
        }
      } else {
        // no read required
        SimicsInterface::RegisterCallback(
          IssueComputeElementCB::Create(this, ce), 0);
      }
    }
  }
}

void
LCAccDevice::SPMReadComplete(ComputeElement* ce)
{
  assert(ce);
  assert(ce->readsRemaining > 0);
  ce->readsRemaining--;

  if (ce->readsRemaining == 0) {
    SimicsInterface::RegisterCallback(
      IssueComputeElementCB::Create(this, ce), 0);
  }
}

void
LCAccDevice::IssueComputeElement(ComputeElement* ce)
{
  assert(ce);
  assert(ce->host);
  assert(ce->host->mode);
  assert(ce->host->computesToEmit > 0);
  assert(ce->host->computesToEmit <= ce->host->computesToFinish);

  uint64_t startOffset =
    (nextComputeAvailable > SimicsInterface::GetSystemTime()) ?
    (nextComputeAvailable - SimicsInterface::GetSystemTime()) : 0;
  uint64_t cycleTime = (mCycleTime <= 0) ?
                       ce->host->mode->CycleTime() : mCycleTime;
  uint64_t depth = (mPipelineDepth <= 0) ?
                   ce->host->mode->PipelineDepth() : mPipelineDepth;
  uint64_t ii = (mInitiationInterval <= 0) ?
                ce->host->mode->InitiationInterval() : mInitiationInterval;

  nextComputeAvailable =
    SimicsInterface::GetSystemTime() + startOffset + ii * cycleTime;
  SimicsInterface::RegisterCallback(
    PerformComputeElementCB::Create(this, ce),
    startOffset + depth * cycleTime);

  if (ce->host->computesToEmit == ce->host->computeCount) {
    SimicsInterface::RegisterCallback(
      AnnounceComputeBeginCB::Create(this, ce), startOffset);
  }

  ce->host->computesToEmit--;
}

void
LCAccDevice::AnnounceComputeBegin(ComputeElement* ce)
{
  assert(ce);
  //ML_LOG(GetDeviceName(), "BEGIN Compute task " << ce->taskID);
}

void
LCAccDevice::PerformComputeElement(ComputeElement* ce)
{
  assert(ce);
  assert(ce->readsRemaining == 0);
  assert(ce->host);
  assert(ce->host->mode);
  assert(ce->host->computesToFinish > 0);

  std::vector<uint64_t> seedAddrs;

  for (size_t i = 0; i < ce->host->argumentAddressGen.size(); i++) {
    seedAddrs.push_back(ce->host->CalcArgAddress(i, ce->index));
  }

  ce->host->mode->Compute(ce->index, ce->maxCompute, ce->taskID,
                          seedAddrs, ce->host->argumentActive);
  std::vector<uint64_t> writeAddrs;
  ce->host->mode->GetSPMWriteIndexSet(ce->index, ce->maxCompute,
                                      ce->taskID, seedAddrs, ce->host->argumentActive, writeAddrs);
  ce->writesRemaining = writeAddrs.size();
  ce->randomAccessesRemaining = ce->host->mode->MemoryAccessCount();
  ce->host->computesToFinish--;

  if (writeAddrs.size() > 0 || ce->randomAccessesRemaining > 0) {
    // spm write required
    for (size_t i = 0; i < writeAddrs.size(); i++) {
      spm->TimedWrite(writeAddrs[i], 1,
                      SPMWriteCompleteCB::Create(this, ce));
    }

    // random dma access required
    for (size_t i = 0; i < ce->randomAccessesRemaining; i++) {
      assert(ce->host->mode->MemoryAccessCount() > 0);
      uint64_t spmAddr;
      uint64_t memAddr;
      uint32_t size;
      int type;
      ce->host->mode->GetNextMemoryAccess(spmAddr, memAddr, size, type);
      spmAddr = spmAddr + ce->host->indexOffset;
      //std::cout << "mySPM" << spm->GetID() << std::endl;
      //std::cout << "spmAddr" << spmAddr << std::endl;
      //std::cout << "memAddr" << memAddr << std::endl;
      //std::cout << "size" << size << std::endl;
      //std::cout << "type" << type << std::endl;
      // ML_LOG(GetDeviceName(), "random access on 0x" << std::hex << memAddr);

      dma->BeginSingleElementTransfer(spm->GetID(),
                                      spmAddr, memAddr, size, type,
                                      ComputeRandomAccessCompleteCB::Create(this, ce, ce->index,
                                          ce->maxCompute, ce->taskID, spmAddr, memAddr, type));
    }
  } else {
    // no write required
    TryRetireComputeElement(ce);
  }
}

void
LCAccDevice::ComputeRandomAccessComplete(ComputeElement* ce,
    int iteration, int maxIteration, int taskID, uint64_t spmAddr,
    uint64_t memAddr, int accessType)
{
  assert(ce);
  assert(ce->randomAccessesRemaining > 0);
  assert(ce->host);
  assert(ce->host->mode);

  std::vector<uint64_t> seedAddrs;

  for (size_t i = 0; i < ce->host->argumentAddressGen.size(); i++) {
    seedAddrs.push_back(ce->host->CalcArgAddress(i, ce->index));
  }

  ce->host->mode->MemoryAccessComplete(iteration, maxIteration,
                                       taskID, spmAddr, memAddr, accessType, seedAddrs);
  size_t newAccesses = ce->host->mode->MemoryAccessCount();
  ce->randomAccessesRemaining += newAccesses;

  for (size_t i = 0; i < newAccesses; i++) {
    assert(ce->host->mode->MemoryAccessCount() > 0);
    uint64_t spmAddr;
    uint64_t memAddr;
    uint32_t size;
    int type;
    ce->host->mode->GetNextMemoryAccess(spmAddr, memAddr, size, type);
    dma->BeginSingleElementTransfer(spm->GetID(),
                                    spmAddr, memAddr, size, type,
                                    ComputeRandomAccessCompleteCB::Create(this, ce, ce->index,
                                        ce->maxCompute, ce->taskID, spmAddr, memAddr, type));
  }

  ce->randomAccessesRemaining--;
  TryRetireComputeElement(ce);
}

void
LCAccDevice::SPMWriteComplete(ComputeElement* ce)
{
  assert(ce);
  assert(ce->writesRemaining > 0);
  ce->writesRemaining--;
  TryRetireComputeElement(ce);
}

void
LCAccDevice::TryRetireComputeElement(ComputeElement* ce)
{
  assert(ce);

  if (ce->writesRemaining == 0 && ce->randomAccessesRemaining == 0) {
    assert(ce->host);
    assert(ce->host->pendingComputeElements > 0);
    ce->host->pendingComputeElements--;

    if (ce->host->pendingComputeElements == 0) {
      //ML_LOG(GetDeviceName(), "END Compute task "  << ce->taskID);
      assert(ce->host->host);
      SimicsInterface::RegisterCallback(
        EmitActionWritesCB::Create(this, ce->host->host), 0);
    }

    delete ce;
  }
}

void
LCAccDevice::EmitActionWrites(Action* a)
{
  assert(a);
  assert(a->readsInFlight == 0 && a->writesInFlight == 0);

  if (a->writeSet.size() > 0) {
    // dma write required
    // ML_LOG(GetDeviceName(), "Begin DMA write task " << a->taskID);
    a->writesInFlight = a->writeSet.size();

    for (size_t i = 0; i < a->writeSet.size(); i++) {
      TransferOrder& t = a->writeSet[i];
      assert(t.src.IsSPM() || t.dst.IsSPM());

      if (t.dst.IsSPM()) {
        //ML_LOG(GetDeviceName(), "BEGIN Transfer Chain Write "
        //       << & (a->writeSet[i]) << " " << t.src.spm << " -> " << t.dst.spm);
      } else {
        //ML_LOG(GetDeviceName(), "BEGIN Transfer Write " << i
        //       << " Task " << a->taskID);
      }

      dma->BeginTransfer(t.src.spm, (uint64_t)t.src.addr,
                         t.src.sizeSet, t.src.strideSet, t.dst.spm,
                         (uint64_t)t.dst.addr, t.dst.sizeSet, t.dst.strideSet,
                         t.elementSize, Priority::Medium,
                         FinishedActionWriteCB::Create(this, a, i));
    }
  } else {
    // no write required
    RetireAction(a);
  }
}

void
LCAccDevice::FinishedActionWrite(Action* a, size_t index)
{
  assert(a);
  assert(a->writesInFlight > 0);
  a->writesInFlight--;
  assert(a->writeSet.size() > index);
  TransferOrder& t = a->writeSet[index];

  if (t.dst.IsSPM()) {
    //ML_LOG(GetDeviceName(), "END Transfer Chain Write "
    //       << & (a->writeSet[index]) << " "  << t.src.spm << " -> " << t.dst.spm);
  } else {
    //ML_LOG(GetDeviceName(), "END Transfer Write " << index
    //       << " Task " << a->taskID);
  }

  if (a->writesInFlight == 0) {
    // done dma write
    // ML_LOG(GetDeviceName(), "End DMA write task " << a->taskID);
    RetireAction(a);
  }
}

void
LCAccDevice::RetireAction(Action* a)
{
  assert(a);
  assert(a->readsInFlight == 0 && a->writesInFlight == 0);
  assert(a->computeReady);
  assert(a->retireReady == false);
  a->retireReady = true;
  assert(!awaitingRetireActionSet.empty());

  while (!awaitingRetireActionSet.empty() &&
         awaitingRetireActionSet.front()->retireReady) {
    Action* curAc = awaitingRetireActionSet.front();
    awaitingRetireActionSet.pop();

    if (curAc->notifyOnComplete) {
      // the last task in the job
      //ML_LOG(GetDeviceName(), "Finish job "
      //       << a->taskID / a->numberOfTasks);

      uint32_t msg[2];
      msg[0] = LCACC_CMD_TASK_COMPLETED;
      msg[1] = currentUser;
      netPort->SendMessage(currentUserProc, &msg, sizeof(msg));
      // for (size_t i = 0; i < prefetcherSet.size(); i++) {
      //     uint8_t stopMsg = 2;//stop prefetching
      //     netPort->SendMessage(prefetcherSet[i], &stopMsg,
      //         sizeof(stopMsg));
      // }
    }

    delete curAc;
    StartTask();
  }
}

void
LCAccDevice::StartTask()
{
  if (awaitingRetireActionSet.size() < spmWindowCount &&
      !pendingActionSet.empty() && pendingReadSet.empty()) {
    // issue width equals to spmwindowcount
    Action* a = pendingActionSet.front();
    pendingActionSet.pop();
    assert(a->compute != NULL || a->readSet.size() > 0 || a->writeSet.size() > 0);
    awaitingComputeActionSet.push(a);
    awaitingRetireActionSet.push(a);
    EmitActionReads(a);
  }

  if (awaitingRetireActionSet.empty() && pendingActionSet.empty()) {
    assert(spm);
    assert(lcaActive);
    spm->Clear();
    lcaActive = false;
  }
}

void 
LCAccDevice::ProtectionTableAccess()
{
  //std::cout << "Protection Table" << protection_table_request.size() << std::endl;
  assert(protection_table_request.size());
   
    int32_t* msg = protection_table_request.front();
    protection_table_request.pop_front();
    int command = msg[0];
    int process = msg[1];
    BitConverter bc_vAddr;
    bc_vAddr.u32[0] = msg[2];
    bc_vAddr.u32[1] = msg[3];
    uint64_t vAddr = bc_vAddr.u64[0];
    BitConverter bc;
    bc.s32[0] = msg[4];
    bc.s32[1] = msg[5];
    uint64_t pAddr = bc.u64[0];
    BitConverter bmc;
    bmc.u32[0] = msg[6];
    bmc.u32[1] = msg[7];
    uint64_t MAC_return = bmc.u64[0];
    BitConverter bno;
    bno.u32[0] = msg[8];
    bno.u32[1] = msg[9];
    uint64_t node_return = bno.u64[0];
    //std::cout << "data" << data << std::endl;
    //std::cout << "spmAddr" << pAddr << std::endl;
    //std::cout << "memAddr" << pAddr << std::endl;
    //uint8_t addr=4294967295;
    int size=1;
    uint8_t* data = new uint8_t [size];
    int cpu=2;
    int buffer;
    //TimedReadComplete(addr,buffer);
    uint64_t memAddr = 2133079400;
    //*data =0;
    //MemoryInterface::Instance()->functionalRead(memAddr, data, size);
    //int buffer;
    //std::cout <<std::hex<< "memAddr" << pAddr << std::endl;
    //std::cout << "node_return" << node_return << std::endl;
    //std::cout << "Failure 1"  << std::endl;
    SimicsInterface::TimedBufferRead(node_return, memAddr, size, buffer,TimedReadCompleteCB::Create(this));
    if(data)
    {
      cout_protection++;
      //std::cout << "Protection verfication Count" << cout_protection << std::endl;
          
          dma->finishTranslation(vAddr, pAddr,MAC_return);
    }
    
    if(protection_table_request.size())
    {
      ProtectionTableAccess();
    }
    /*
    std::cout<<"data1" << (int*)data << std::endl;
    std::cout<<"data2" << (int*)(data+1) << std::endl;
    delete [] data;
    */
    

  
  
}
void 
LCAccDevice::TimedReadComplete()
{    
    
    
    /*
    if(data)
    {
    std::cout << "data" << data << std::endl;
    std::cout << "address" << addr << std::endl;
    dma->finishTranslation(vAddr, pAddr,MAC_return);
    }
    */
    
    if(protection_table_request.size())
    {
      ProtectionTableAccess();
    }
    
    
}


void 
LCAccDevice::ProtectionTableRead()
{
   //std::cout << "Protection Table Read" << std::endl;
}


void LCAccDevice::MsgHandler(int src, const void* packet, unsigned int packetSize)
{
  //std::cout << "Size of packet " << packetSize << std::endl;
  assert(packetSize >= 9);
  int32_t* msg = (int32_t*)packet;
  int command = msg[0];
  int process = msg[1];
  // std::cout << "Accelerator " << netPort->GetNetworkPort() << " recvd from " << src << " for " << process << std::endl;

  switch (command) {
  case (GAM_CMD_ACCNOTIFY): {
    assert(packetSize == sizeof(int32_t) * 4);
    assert(pendingActionSet.size() == 0);
    assert(pendingTasksToRead.size() == 0);
    //std::cout << "GAM notify recved " << netPort->GetNetworkPort() << " for " << process << " from " << src << std::endl;
    currentUser = process;
    bufferID = msg[2];
    bufferSize = msg[3];//unused for the time being
    assert(bufferID == -1 || bufferSize > 0);
    dma->FlushTLB();

    spm->SetBuffer(bufferID);
    dma->SetBuffer(bufferID);
    int32_t returnMessage[2];
    returnMessage[0] = GAM_CMD_ACCCONFIRM;
    returnMessage[1] = process;
    netPort->SendMessage(src, returnMessage, sizeof(int32_t) * 2);
  }
  break;

  case (LCACC_CMD_BEGIN_TASK_SIGNATURE): {
    assert(!lcaActive);
    lcaActive = true;
    assert(packetSize == 10 * sizeof(int32_t));
    // std::cout << "process:" << process << std::endl;
    assert(process == currentUser);
    TaskBlock t;
    currentUserProc = src;
    assert(currentUserProc != 0);
    BitConverter bc_vAddr;
    bc_vAddr.u32[0] = msg[2];
    bc_vAddr.u32[1] = msg[3];
    t.startingAddress = bc_vAddr.u64[0];
    BitConverter bc;
    bc.u32[0] = msg[4];
    bc.u32[1] = msg[5];
    uint64_t pAddr = bc.u64[0];

    if (t.startingAddress) {
      assert(pAddr);
      dma->AddTLBEntry(t.startingAddress, pAddr);
    }

    t.bufferSize = msg[6];
    //std::cout << "Begin task signature recved on " << netPort->GetNetworkPort() << " " << t.startingAddress << " -- " << pAddr << " " << t.bufferSize << std::endl;
    t.lowMark = 1000000000;//msg[4];
    t.taskType = tt_signature;
    assert(pendingTasksToRead.size() == 0);
    pendingTasksToRead.push(t);
    StartTaskRead();
  }
  break;

  case (LCACC_CMD_BEGIN_EMBEDDED_TASK_SIGNATURE): {

    // std::ofstream of("Msg", std::ios::out | std::ios::app);
    // of << "******* Task Recieved" << std::endl;
    // of << src << std::endl;
    // of << netPort->GetNetworkPort() << std::endl;
    // of << packetSize << std::endl;
    // for(size_t i = 0; i < packetSize; i++)
    // {
    //     of << " " << (int)(((uint8_t*)packet)[i]) << std::endl;
    // }
    // of << std::endl;
    // of.close();

    assert(!lcaActive);
    lcaActive = true;
    currentUserProc = src;
    currentUser = process;
    bufferID = msg[2];
    bufferSize = msg[3];//unused for the time being
    assert(bufferID == -1 || bufferSize > 0);

    // dma->FlushTLB();
    // ML_LOG(GetDeviceName(), "Flushing DMA TLB");

    spm->SetBuffer(bufferID);
    dma->SetBuffer(bufferID);
    PacketReader pr(&(msg[4]), packetSize - sizeof(int32_t) * 4);
    ParseTaskSignature(pr);
    assert(pr.SizeRemaining() < 4);
  }
  break;

  case (LCACC_CMD_TLB_SERVICE): {
 
    currentUserProc = src;
    BitConverter bc_vAddr;
    bc_vAddr.u32[0] = msg[2];
    bc_vAddr.u32[1] = msg[3];
    uint64_t vAddr = bc_vAddr.u64[0];
    BitConverter bc;
    bc.s32[0] = msg[4];
    bc.s32[1] = msg[5];
    uint64_t pAddr = bc.u64[0];
    BitConverter bmc;
    bmc.u32[0] = msg[6];
    bmc.u32[1] = msg[7];
    uint64_t MAC_return = bmc.u64[0];
    BitConverter bno;
    bno.u32[0] = msg[8];
    bno.u32[1] = msg[9];
    uint64_t node_return = bno.u64[0];

    //ML_LOG(GetDeviceName(), "END TLB miss 0x" << std::hex << vAddr);
    //ML_LOG(GetDeviceName(), "TLB miss serviced 0x" << std::hex
    //    << vAddr << " -> 0x" << pAddr);
    dma->finishTranslation(vAddr, pAddr,MAC_return);
  
    
  }
  break;

  case (LCACC_CMD_CANCEL_TASK): {
    assert(0);//not yet supported
  }
  break;

  case (DMA_MEMORY_REQUEST):
  case (DMA_MEMORY_RESPONSE):
    break;

  default:
    std::cout << "ERROR :" << command << std::endl;
    assert(0);
  }//*/
}


void LCAccDevice::ReadTaskSignatureDone()
{
  ML_LOG(GetDeviceName(), "End task list read");
  TaskBlock task = pendingTasksToRead.front();
  pendingTasksToRead.pop();
  uint8_t* buf = new uint8_t[task.bufferSize];
  spm->Read(spmTaskLoadAddr, task.bufferSize, buf);
  PacketReader pr(buf, task.bufferSize);
  delete [] buf;
  pr.SetEndianSwap(false);
  //pr.SetEndianSwap(true);

  ParseTaskSignature(pr);
  assert(pr.SizeRemaining() < 4);//packet size is word-rounded upward, and parsed packet may not be.  Thus the only thing that can be said the remaining un-parsed segment should not exceed the rounding amount.
}

void
LCAccDevice::ParseTaskSignature(PacketReader& pr)
{
  std::queue<std::vector<TransferOrder> > prefetchQueue;
  std::queue<std::vector<TransferOrder> > readQueue;
  std::queue<std::vector<TransferOrder> > writeQueue;
  std::queue<std::vector<ComputeOrder> > computeQueue;
  unsigned int transfersPerTask;
  unsigned int computesPerTask;
  unsigned int skipTasks;
  unsigned int numberOfTasks;
  spmWindowSize = 0;
  spmWindowCount = 0;
  std::vector<PolyhedralAddresser> transferSrcStartAddr;
  std::vector<PolyhedralAddresser> transferDstStartAddr;
  std::vector<std::vector<uint32_t> > transferSrcSize;
  std::vector<std::vector<int32_t> > transferSrcStride;
  std::vector<std::vector<uint32_t> > transferDstSize;
  std::vector<std::vector<int32_t> > transferDstStride;
  std::vector<size_t> transferElementSize;
  std::vector<int> transferSrcSPM;
  std::vector<int> transferDstSPM;
  std::vector<size_t> transferMoveSize;
  std::vector<ComputeOrder> computeSet;

  //decode signature
  transfersPerTask = pr.Read<uint8_t>();
  computesPerTask = pr.Read<uint8_t>();
  skipTasks = pr.Read<uint32_t>();
  numberOfTasks = pr.Read<uint32_t>();
  //ML_LOG(GetDeviceName(), "Receiving job " << skipTasks / 8);
  //ML_LOG(GetDeviceName(), "transfers per task: " << transfersPerTask);
  //ML_LOG(GetDeviceName(), "computes per task: " << computesPerTask);
  //ML_LOG(GetDeviceName(), "tasks in job: " << skipTasks << " - "
  //       << numberOfTasks + skipTasks - 1);

  assert(transfersPerTask > 0);
  assert(computesPerTask > 0);
  assert(numberOfTasks > 0);
  std::vector<TransferDescription> tdSet;
  std::vector<ComputeDescription> cdSet;

  for (unsigned int i = 0; i < computesPerTask; i++) {
    ComputeDescription cd;
    cd.ReadIn(pr);
    // FIXME
    // cd.spmWindowCount = 2;
    assert(modes.find(cd.opcode) != modes.end());
    std::map<unsigned int, unsigned int> addrMap;
    CalculateSPMMapping(spmBanks, spmBankMod, spmReadPortCount,
                        spmWritePortCount, cd, addrMap);
    spm->SetAddressMap(addrMap);
    cdSet.push_back(cd);
  }

  for (unsigned int i = 0; i < transfersPerTask; i++) {
    TransferDescription td;
    td.ReadIn(pr);
    tdSet.push_back(td);
  }

  if (pr.SizeRemaining() > 2 * sizeof(uint32_t) + sizeof(uint64_t)) {
    // TLB preload entries
    uint32_t tlbCounter = pr.Read<uint32_t>();
    //ML_LOG(GetDeviceName(), "Preload TLB with "
    //       << tlbCounter << "entries");

    for (uint32_t i = 0; i < tlbCounter; i++) {
      uint32_t logicalPage = pr.Read<uint32_t>();
      uint64_t physicalPage = pr.Read<uint64_t>();
      //ML_LOG(GetDeviceName(), "0x" << std::hex << logicalPage
      //       << " -> 0x" << std::hex << physicalPage);
      dma->AddTLBEntry(logicalPage, physicalPage);
    }
  }

  for (unsigned int i = 0; i < transfersPerTask; i++) {
    TransferDescription& td = tdSet[i];
    std::vector<uint32_t> srcSizePrefix, srcSizeSuffix;
    std::vector<int32_t> srcStridePrefix, srcStrideSuffix;
    std::vector<uint32_t> dstSizePrefix, dstSizeSuffix;
    std::vector<int32_t> dstStridePrefix, dstStrideSuffix;
    int src, dst;
    uint64_t srcBaseAddr, dstBaseAddr;
    unsigned int srcDimensionality, dstDimensionality;
    unsigned int srcDimensionalityLower, dstDimensionalityLower;
    unsigned int srcDimensionalityUpper, dstDimensionalityUpper;
    size_t elementSize;

    src = td.srcDevice;
    srcBaseAddr = td.srcBaseAddress;
    srcDimensionalityLower = td.srcSplit;
    srcDimensionalityUpper = td.srcSize.size() - td.srcSplit;;
    dst = td.dstDevice;
    dstBaseAddr = td.dstBaseAddress;
    dstDimensionalityLower = td.dstSplit;
    dstDimensionalityUpper = td.dstSize.size() - td.dstSplit;
    srcDimensionality = srcDimensionalityLower + srcDimensionalityUpper;
    dstDimensionality = dstDimensionalityLower + dstDimensionalityUpper;

    // std::cout << "Src: " << src << " from " << srcBaseAddr << " dim " << srcDimensionality << " split: " << srcDimensionalityLower << "/" << srcDimensionalityUpper << std::endl;
    // std::cout << "Dst: " << dst << " from " << dstBaseAddr << " dim " << dstDimensionality << " cutoff " << dstDimensionalityLower << "/" << dstDimensionalityUpper << std::endl;
    for (unsigned int x = 0; x < srcDimensionalityLower; x++) {
      srcSizePrefix.push_back(td.srcSize[x]);
      srcStridePrefix.push_back(td.srcStride[x]);
    }

    for (unsigned int x = 0; x < srcDimensionalityUpper; x++) {
      srcSizeSuffix.push_back(td.srcSize[x + srcDimensionalityLower]);
      srcStrideSuffix.push_back(td.srcStride[x + srcDimensionalityLower]);
    }

    for (unsigned int x = 0; x < dstDimensionalityLower; x++) {
      dstSizePrefix.push_back(td.dstSize[x]);
      dstStridePrefix.push_back(td.dstStride[x]);
    }

    for (unsigned int x = 0; x < dstDimensionalityUpper; x++) {
      dstSizeSuffix.push_back(td.dstSize[x + dstDimensionalityLower]);
      dstStrideSuffix.push_back(td.dstStride[x + dstDimensionalityLower]);
    }

    elementSize = td.elementSize;
    // std::cout << "Src manifest ";
    // for(size_t x = 0; x < srcSizePrefix.size(); x++) std::cout << "[" << srcSizePrefix[x] << "," << srcStridePrefix[x] << "] ";
    // for(size_t x = 0; x < srcSizeSuffix.size(); x++) std::cout << "[" << srcSizeSuffix[x] << "," << srcStrideSuffix[x] << "] ";
    // std::cout << std::endl;
    // std::cout << "Dst manifest ";
    // for(size_t x = 0; x < dstSizePrefix.size(); x++) std::cout << "[" << dstSizePrefix[x] << "," << dstStridePrefix[x] << "] ";
    // for(size_t x = 0; x < dstSizeSuffix.size(); x++) std::cout << "[" << dstSizeSuffix[x] << "," << dstStrideSuffix[x] << "] ";
    // std::cout << std::endl;
    // std::cout << "Element size " << elementSize << std::endl;
    unsigned int srcMoveSize = elementSize, dstMoveSize = elementSize;

    for (size_t x = 0; x < srcSizeSuffix.size(); x++) {
      srcMoveSize *= srcSizeSuffix[x];
    }

    for (size_t x = 0; x < dstSizeSuffix.size(); x++) {
      dstMoveSize *= dstSizeSuffix[x];
    }

    assert(srcMoveSize == dstMoveSize);
    transferSrcStartAddr.push_back(PolyhedralAddresser(srcBaseAddr, srcSizePrefix, srcStridePrefix));
    transferDstStartAddr.push_back(PolyhedralAddresser(dstBaseAddr, dstSizePrefix, dstStridePrefix));
    transferSrcSize.push_back(srcSizeSuffix);
    transferSrcStride.push_back(srcStrideSuffix);
    transferDstSize.push_back(dstSizeSuffix);
    transferDstStride.push_back(dstStrideSuffix);
    transferSrcSPM.push_back(src);
    transferDstSPM.push_back(dst);
    transferElementSize.push_back(elementSize);
    transferMoveSize.push_back(srcMoveSize);
    // if( (src == NO_SPM_ID) || ((dst == NO_SPM_ID) && (prefetcherSet.size() > 0)) )
    // {//configure prefetchers for this stream
    //     PacketBuilder pb;
    //     uint64_t baseAddr;
    //     uint8_t direction;
    //     std::vector<int32_t>* stridePref;
    //     std::vector<int32_t>* strideSuf;
    //     std::vector<uint32_t>* sizePref;
    //     std::vector<uint32_t>* sizeSuf;
    //     if(src == NO_SPM_ID)
    //     {
    //         direction = 0;//a read
    //         baseAddr = srcBaseAddr;
    //         stridePref = &srcStridePrefix;
    //         strideSuf = &srcStrideSuffix;
    //         sizePref = &srcSizePrefix;
    //         sizeSuf = &srcSizeSuffix;
    //     }
    //     else
    //     {
    //         direction = 1;//a write
    //         baseAddr = dstBaseAddr;
    //         stridePref = &dstStridePrefix;
    //         strideSuf = &dstStrideSuffix;
    //         sizePref = &dstSizePrefix;
    //         sizeSuf = &dstSizeSuffix;
    //     }
    //     pb.Write(direction);
    //     pb.Write(baseAddr);
    //     pb.Write((uint8_t)(sizePref->size()));
    //     pb.Write((uint8_t)(sizeSuf->size()));
    //     for(size_t i = 0; i < sizePref->size(); i++)
    //     {
    //         pb.Write((*sizePref)[i]);
    //         pb.Write((*stridePref)[i]);
    //     }
    //     for(size_t i = 0; i < sizeSuf->size(); i++)
    //     {
    //         pb.Write((*sizeSuf)[i]);
    //         pb.Write((*strideSuf)[i]);
    //     }
    //     for(size_t i = 0; i < prefetcherSet.size(); i++)
    //     {
    //         netPort->SendMessage(prefetcherSet[i], pb.GetBuffer(), pb.GetBufferSize());
    //     }
    // }
  }

  for (unsigned int i = 0; i < computesPerTask; i++) {
    ComputeDescription& cd = cdSet[i];
    assert(spmWindowSize == 0);
    spmWindowSize = cd.spmWindowSize;
    spmWindowCount = cd.spmWindowCount;
    unsigned int argCount;
    unsigned int registerCount;
    std::vector<std::vector<uint32_t> > indexSize;
    std::vector<std::vector<int32_t> > indexStride;
    std::vector<uint64_t> baseIndex;
    std::vector<int64_t> controlRegisterValue;
    argCount = cd.argBaseAddr.size();

    //ML_LOG(GetDeviceName(), "spmwindowcount: " << spmWindowCount);
    //ML_LOG(GetDeviceName(), "compute opcode: " << cd.opcode
    //       << " f(" << argCount << ")");

    for (unsigned int x = 0; x < argCount; x++) {
      baseIndex.push_back(cd.argBaseAddr[x]);
      unsigned int dimensionality = cd.argSize[x].size();
      std::vector<uint32_t> size;
      std::vector<int32_t> stride;

      for (unsigned int y = 0; y < dimensionality; y++) {
        size.push_back(cd.argSize[x][y]);
        stride.push_back(cd.argStride[x][y]);
      }

      indexSize.push_back(size);
      indexStride.push_back(stride);

      //std::cout << "Arg " << x << " base:" << baseIndex[x] << " ";

      for (size_t i = 0; i < size.size(); i++) {
       // std::cout << "size: " << size[i]
        //          << " stride: " << stride[i] << " ";
      }

      //std::cout << "element size: " << cd.argElementSize[x] << std::endl;
    }

    registerCount = cd.controlRegister.size();
    //std::cout << "Register count: " << registerCount << " ";

    for (unsigned int x = 0; x < registerCount; x++) {
      controlRegisterValue.push_back(cd.controlRegister[x]);
      //std::cout << "[" << (int32_t)controlRegisterValue[x] << "] ";
    }

    std::cout << std::endl;
    unsigned int computeSize = 0;
    std::vector<PolyhedralAddresser> argAddr;
    std::vector<uint32_t> argSize;
    std::vector<uint64_t> registers;

    for (size_t x = 0; x < indexSize.size(); x++) {
      PolyhedralAddresser pac(baseIndex[x], indexSize[x],
                              indexStride[x]);
      computeSize = (pac.TotalSize() > computeSize) ?
                    pac.TotalSize() : computeSize;
      argAddr.push_back(pac);
    }

    assert(computeSize > 0);

    for (size_t x = 0; x < controlRegisterValue.size(); x++) {
      registers.push_back(controlRegisterValue[x]);
    }

    modes[cd.opcode]->SetSPM(spm);
    modes[cd.opcode]->BeginComputation();
    modes[cd.opcode]->SetArgumentWidth(cd.argElementSize);
    modes[cd.opcode]->SetRegisterValues(registers);
    computeSet.push_back(ComputeOrder(modes[cd.opcode],
                                      computeSize, computeSet.size(),
                                      (computeSet.size() % spmWindowCount) * spmWindowSize,
                                      argAddr, cd.argActive, registers));
  }

  // fill queues
  for (unsigned int i = skipTasks; i < skipTasks + numberOfTasks; i++) {
    std::vector<TransferOrder> reads;
    std::vector<TransferOrder> writes;
    std::vector<ComputeOrder> computes;
    std::vector<TransferOrder> prefetches;

    for (unsigned int x = 0; x < transfersPerTask; x++) {
      TransferOrder t;
      t.src.spm = transferSrcSPM[x];
      t.dst.spm = transferDstSPM[x];
      t.src.sizeSet = transferSrcSize[x];
      t.src.strideSet = transferSrcStride[x];
      t.dst.sizeSet = transferDstSize[x];
      t.dst.strideSet = transferDstStride[x];
      t.src.addr = transferSrcStartAddr[x].GetAddr(i);
      t.dst.addr = transferDstStartAddr[x].GetAddr(i);
      t.elementSize = transferElementSize[x];

      if (t.dst.spm == spm->GetID()) {
        reads.push_back(t);
      } else {
        writes.push_back(t);
      }

      if (prefetchDistance > 0) {
        prefetches.push_back(t);
      }
    }

    assert(spmWindowSize > 0);
    assert(spmWindowCount > 1);

    for (size_t x = 0; x < computeSet.size(); x++) {
      ComputeOrder co = computeSet[x];
      co.indexOffset += (i % spmWindowCount) * spmWindowSize;
      co.taskID = i;
      computes.push_back(co);
    }

    readQueue.push(reads);
    writeQueue.push(writes);
    computeQueue.push(computes);

    if (prefetchDistance > 0) {
      prefetchQueue.push(prefetches);
    }
  }

  // translate queue data into actions
  int taskID = skipTasks;
  assert(readQueue.size() == writeQueue.size() &&
         readQueue.size() == computeQueue.size());

  for (int i = 0; i < prefetchDistance && prefetchQueue.size() > 0; i++) {
    prefetchQueue.pop();
  }

  while (!readQueue.empty()) {
    Action* a = new Action(taskID);
    a->numberOfTasks = numberOfTasks;
    // set the flag for the last task in the job
    a->notifyOnComplete = (readQueue.size() > 1) ? false : true;
    a->user = currentUser;
    std::vector<TransferOrder>& read = readQueue.front();
    std::vector<TransferOrder>& write = writeQueue.front();
    std::vector<ComputeOrder>& compute = computeQueue.front();

    for (size_t i = 0; i < read.size(); i++) {
      a->readSet.push_back(read[i]);
    }

    for (size_t i = 0; i < write.size(); i++) {
      a->writeSet.push_back(write[i]);
    }

    if (prefetchQueue.size() > 0) {
      std::vector<TransferOrder>& prefetch = prefetchQueue.front();

      for (size_t i = 0; i < prefetch.size(); i++) {
        a->prefetchSet.push_back(prefetch[i]);
      }

      prefetchQueue.pop();
    }

    assert(compute.size() == 1);

    for (size_t i = 0; i < compute.size(); i++) {
      // These are hackish changes because I'm not sure how to deal
      // with multiple computes in a single effort under the new
      // system. It comes with some problems, such as when to set
      // register state if I allow near-by computes to blend due to
      // SPM MSHR's.  To circumvent this problem, I'm restricting to
      // a single compute per LCA per task.  Thus, the register state
      // stays the same for the entire duration of the task.
      a->compute = new ComputeOrder(compute[i]);
      a->compute->host = a;
    }

    pendingActionSet.push(a);
    readQueue.pop();
    writeQueue.pop();
    computeQueue.pop();
    taskID++;
  }

  if (!pendingTasksToRead.empty()) {
    StartTaskRead();
  } else {
    //ML_LOG(GetDeviceName(), "Starting job "
    //       << skipTasks / numberOfTasks);
    StartTask();
  }
}

void LCAccDevice::StartTaskRead()
{
  //ML_LOG(GetDeviceName(), "Begin task list read");
  TaskBlock t = pendingTasksToRead.front();
  assert(t.taskType == tt_signature);
  std::vector<unsigned int> size;
  std::vector<int> stride;
  size.push_back(t.bufferSize);
  stride.push_back(1);
  //std::cout << "StartTaskRead: TaskDescription" << std::endl;
  /*uint32_t* data = (uint32_t *)t.startingAddress;
  for(int i = 0; i < t.bufferSize / 4; i++)
    std::cout << data[i] << ", ";
    std::cout << std::endl;//*/
  dma->BeginTransfer(NO_SPM_ID, t.startingAddress, size, stride,
                     spm->GetID(), spmTaskLoadAddr, size, stride, 1, Priority::Low,
                     ReadTaskSignatureDoneCB::Create(this));
}
void LCAccDevice::HandleTLBMiss(uint64_t addr, uint64_t MAC_code, uint64_t paddr)
{
  uint32_t msg[10];
  msg[0] = LCACC_CMD_TLB_MISS;
  msg[1] = currentUser;
  BitConverter bc;
  bc.u64[0] = addr;
  msg[2] = bc.u32[0];
  msg[3] = bc.u32[1];
  BitConverter bcp;
  bcp.u64[0] = paddr;
  msg[4] = bcp.u32[0];
  msg[5] = bcp.u32[1];
  BitConverter bmc;
  bmc.u64[0] = MAC_code;
  msg[6] = bmc.u32[0];
  msg[7] = bmc.u32[1];
  BitConverter bno;
  bno.u64[0] = netPort->GetNodeID();
  msg[8] = bno.u32[0];
  msg[9] = bno.u32[1];
  //std::cout << "Current User Proc" << currentUserProc << std::endl;
  //std::cout << "Current User " << currentUser << std::endl;
  //std::cout << "Node Id Handle TLB miss " << netPort->GetNodeID() << std::endl;
    //assert(msg[2] == addr);
  //ML_LOG(GetDeviceName(), "meta vaddr: paddr : MAC_code " << std::hex << addr << paddr << MAC_code);
  //ML_LOG(GetDeviceName(), "BEGIN TLB miss 0x" << std::hex << addr);
  netPort->SendMessage(currentUserProc, msg, sizeof(msg));
  
  
}

void LCAccDevice::HandleMACver(uint64_t addr)
{
  uint32_t msg[4];
  msg[0] = IOMMU_MACver_REQUEST;
  msg[1] = currentUser;
  BitConverter bc;
  bc.u64[0] = addr;
  msg[2] = bc.u32[0];
  msg[3] = bc.u32[1];
  //assert(msg[2] == addr);

  //ML_LOG(GetDeviceName(), "BEGIN TLB miss 0x" << std::hex << addr);
  netPort->SendMessage(currentUserProc, msg, sizeof(msg));
}
void LCAccDevice::HandleAccessViolation(uint64_t)
{
  //std::cerr << "[ERROR] Access violation picked up in LCAccDevice." << std::endl;
  assert(0);
}
void LCAccDevice::AddOperatingMode(const std::string& name)
{
  assert(nameToMode.find(name) != nameToMode.end());
  int modeID = nameToMode[name];
  LCAccOperatingMode* mode = modeFactory.Create(name);
  assert(mode);
  assert(modes.find(modeID) == modes.end());
  modes[modeID] = mode;
  lcaccDeviceName = name;
  //ML_LOG(GetDeviceName(), "LCAcc registered at port "
  //       << netPort->GetNodeID());
}
void LCAccDevice::AddNetworkInterface(NetworkInterface* port)
{
  assert(netPort == NULL);
  netPort = port;
}
void LCAccDevice::AddProgrammablePrefetcher(int networkID)
{
  prefetcherSet.push_back(networkID);
}

void
LCAccDevice::CalculateSPMMapping(int bankCount, int bankMod,
                                 int readPorts, int writePorts, const ComputeDescription& compute,
                                 std::map<unsigned int, unsigned int>& addrMap)
{
  addrMap.clear();
  std::vector<PolyhedralAddresser> addrStream;
  assert(modeSampleSet.find(compute.opcode) != modeSampleSet.end());
  LCAccOperatingMode* m = modeSampleSet[compute.opcode];
  assert(m);
  m->SetRegisterValues(compute.controlRegister);
  m->SetArgumentWidth(compute.argElementSize);
  size_t elmCount = 0;
  std::vector<uint64_t> addrSeed;

  for (size_t i = 0; i < compute.argSize.size(); i++) {
    addrStream.push_back(PolyhedralAddresser(compute.argBaseAddr[i],
                         compute.argSize[i], compute.argStride[i]));
    elmCount = (elmCount > addrStream[i].TotalSize()) ?
               elmCount : addrStream[i].TotalSize();
    addrSeed.push_back(0);
  }

  std::vector<int> bankPopulation;

  for (int i = 0; i < bankCount; i++) {
    bankPopulation.push_back(0);
  }

  std::vector<uint64_t> addrResult;

  for (size_t i = 0; i < elmCount; i++) {
    for (size_t x = 0; x < bankPopulation.size(); x++) {
      bankPopulation[x] = 0;
    }

    for (size_t x = 0; x < addrStream.size(); x++) {
      addrSeed[x] = addrStream[x].GetAddr(i);
    }

    addrResult.clear();
    m->GetSPMReadIndexSet(i, elmCount, 0, addrSeed,
                          compute.argActive, addrResult);
    m->GetSPMWriteIndexSet(i, elmCount, 0, addrSeed,
                           compute.argActive, addrResult);

    for (size_t x = 0; x < addrResult.size(); x++) {
      unsigned int addr = addrResult[x];
      assert((uint64_t)addr == addrResult[x]);

      if (addrMap.find(addr) != addrMap.end()) {
        bankPopulation[(addrMap[addr] / bankMod) % bankCount]++;
      } else {
        size_t bankID = 0;

        for (size_t z = 1; z < bankPopulation.size(); z++) {
          if (bankPopulation[bankID] > bankPopulation[z]) {
            bankID = z;
          }
        }

        bankPopulation[bankID]++;

        //std::cout << "Adding address map " << i << "/" << elmCount << "-" << x << "/" << addrResult.size() << ": addr=" << addr << std::endl;
        for (uint32_t spmIndex = 0; spmIndex < compute.spmWindowCount; spmIndex++) {
          addrMap[addr + spmIndex * compute.spmWindowSize] =
            bankMod * bankID;
        }
      }
    }
  }
}

void
LCAccDevice::Initialize(uint32_t threadID)
{
  assert(netPort);

  spm = new SPMInterface(GetDeviceName(), netPort->GetNodeID(),
                         netPort->GetNetworkPort(), 1, spmBanks, spmBankMod,
                         spmReadPortCount, spmReadPortLatency, spmWritePortCount,
                         spmWritePortLatency);
  dma = new DMAController(netPort, spm, HandleTLBMissCB::Create(this),
                          HandleAccessViolationCB::Create(this),HandleMACverCB::Create(this));
  std::vector<int> m;

  for (std::map<int, LCAccOperatingMode*>::iterator it = modes.begin(); it != modes.end(); it++) {
    m.push_back(it->first);
  }

  SimicsInterface::RegisterLCAcc(netPort->GetNetworkPort(), netPort->GetNodeID(), m, threadID);
  netPort->RegisterRecvHandler(MsgHandlerCB::Create(this));
}

bool
LCAccDevice::IsIdle()
{
  return !lcaActive;
}

void
LCAccDevice::SetPrefetchDistance(int dist)
{
  prefetchDistance = dist;
}

void
LCAccDevice::SetSPMConfig(int banks, int bankMod, int readPorts,
                          int readLatency, int writePorts, int writeLatency)
{
  assert(banks > 0);
  assert(bankMod > 0);
  assert(readPorts > 0);
  assert(writePorts > 0);
  assert(readLatency >= 0);
  assert(writeLatency >= 0);
  spmBanks = banks;
  spmBankMod = bankMod;
  spmReadPortCount = readPorts;
  spmWritePortCount = writePorts;
  spmReadPortLatency = readLatency;
  spmWritePortLatency = writeLatency;

  //ML_LOG(GetDeviceName(),
  //       "set SPM banks: " << banks << " bankmod: " << bankMod <<
  //       " readport: " << readPorts << " writeports: " << writePorts <<
  //       " readLatency: " << readLatency << " writeLatency: " << writeLatency);
}

void
LCAccDevice::SetTiming(unsigned int ii, unsigned int pipeline,
                       unsigned int cycle)
{
  mCycleTime = cycle;
  mPipelineDepth = pipeline;
  mInitiationInterval = ii;
  //ML_LOG(GetDeviceName(), "set cycle: " << cycle << " pipeline: "
  //       << pipeline << " ii: " << ii);
}

void
LCAccDevice::HookToMemoryDevice(const std::string& name)
{
  assert(memoryDeviceHookName == NO_DEVICE_ID);
  assert(name != NO_DEVICE_ID);
  memoryDeviceHookName = name;
}

LCAccDevice::LCAccDevice()
{
  mCycleTime = -1;
  mPipelineDepth = -1;
  mInitiationInterval = -1;
  prefetchDistance = 0;
  spmBanks = 1;
  spmBankMod = 1;
  spmWritePortCount = spmReadPortCount = 9999999;
  spmReadPortLatency = spmWritePortLatency = 0;
  memoryDeviceHookName = NO_DEVICE_ID;

  if (!initializedModeFactory) {
#define REGISTER_OPMODE(x) \
            {\
                std::cout << "Adding operating mode " << x::GetModeName() << " with mode # " << x::GetOpCode() << std::endl;\
                modeFactory.RegisterGenerator(x::GetModeName(), new ObjFactory<std::string, LCAccOperatingMode>::TemplateGenerator<x>);\
                assert(nameToMode.find(x::GetModeName()) == nameToMode.end()); \
                nameToMode[x::GetModeName()] = x::GetOpCode();\
                assert(modeToName.find(x::GetOpCode()) == modeToName.end()); \
                assert(modeSampleSet.find(x::GetOpCode()) == modeSampleSet.end()); \
                modeToName[x::GetOpCode()] = x::GetModeName(); \
                modeSampleSet[x::GetOpCode()] = modeFactory.Create(x::GetModeName()); \
            }
    std::cout << "Initializing operating mode factory" << std::endl;
#include "LCAccOperatingModeListing.hh"
    initializedModeFactory = true;
  }

  nextComputeAvailable = 0;
  spm = NULL;
  netPort = NULL;
  dma = NULL;
  spmTaskLoadAddr = 0;
  lcaActive = false;
}
LCAccDevice::~LCAccDevice()
{
  assert(!lcaActive);
  SimicsInterface::UnregisterLCAcc(netPort->GetNetworkPort());
  assert(spm);
  assert(dma);
  delete spm;
  delete dma;
  assert(netPort);
  delete netPort;
}
}

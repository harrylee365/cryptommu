/*
 * Copyright (c) 1999-2011 Mark D. Hill and David A. Wood
 * All rights reserved.
 *
 * Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are
 * met: redistributions of source code must retain the above copyright
 * notice, this list of conditions and the following disclaimer;
 * redistributions in binary form must reproduce the above copyright
 * notice, this list of conditions and the following disclaimer in the
 * documentation and/or other materials provided with the distribution;
 * neither the name of the copyright holders nor the names of its
 * contributors may be used to endorse or promote products derived from
 * this software without specific prior written permission.
 *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 */

#include <fcntl.h>
#include <zlib.h>

#include <cstdio>

#include "base/intmath.hh"
#include "base/statistics.hh"
#include "debug/RubyCacheTrace.hh"
#include "debug/RubySystem.hh"
#include "mem/ruby/common/Address.hh"
#include "mem/ruby/network/Network.hh"
#include "mem/ruby/system/System.hh"
#include "sim/eventq.hh"
#include "sim/simulate.hh"
#include "modules/Synchronize/Synchronize.hh"

using namespace std;
//0. SW version
// #define SIM_SW

//1. ARC version
//#define SIM_ARC

//2. ARC PLUS version
#define SIM_ARC_PLUS

//3. CHARM NON-ISLAND version
//#define SIM_CHARM

//4. CHARM ISLAND version
//#define SIM_CHARM_ISLAND

//for compiling purpose
#if defined(SIM_SW)
#define SIM_ARC
#endif

#if defined(SIM_ARC_PLUS) || defined(SIM_CHARM) || defined(SIM_CHARM_ISLAND)
#define SIM_TD
#endif

#define SIM_NET_PORTS
#ifdef SIM_NET_PORTS
#include <iostream>
#include "SimicsNetworkPortInterface.hh"
#include "SimicsNetworkPort.hh"
#include "sim/system.hh"
#endif

int RubySystem::m_random_seed;
bool RubySystem::m_randomization;
uint32_t RubySystem::m_block_size_bytes;
uint32_t RubySystem::m_block_size_bits;
uint32_t RubySystem::m_memory_size_bits;
uint32_t RubySystem::m_td_tlb_size;
uint32_t RubySystem::m_td_tlb_latency;
uint32_t RubySystem::m_td_tlb_assoc;
uint32_t RubySystem::m_lcacc_tlb_size;
uint32_t RubySystem::m_lcacc_tlb_mshr;
uint32_t RubySystem::m_lcacc_tlb_latency;
uint32_t RubySystem::m_lcacc_tlb_assoc;
uint32_t RubySystem::m_dma_issue_width;
#ifdef SIM_NET_PORTS
std::vector<std::string> RubySystem::accTypes;
int RubySystem::m_num_simics_net_ports;
int RubySystem::m_num_accelerators;
int RubySystem::m_num_TDs;
int RubySystem::m_num_acc_instances;
int RubySystem::m_host_ptw_latency;
bool RubySystem::m_tlb_hack;
bool RubySystem::m_ideal_mmu;
bool RubySystem::m_iommu;
/* flags to simulate three varaints of CryptoMMU */
bool RubySystem::m_cryptommu_serial;
bool RubySystem::m_cryptommu_parallel;
bool RubySystem::m_cryptommu_parallel_ra;

#endif

RubySystem::RubySystem(const Params *p)
    : ClockedObject(p)
{
    if (g_system_ptr != NULL)
        fatal("Only one RubySystem object currently allowed.\n");

    m_random_seed = p->random_seed;
    srandom(m_random_seed);
    m_randomization = p->randomization;

    m_block_size_bytes = p->block_size_bytes;
    assert(isPowerOf2(m_block_size_bytes));
    m_block_size_bits = floorLog2(m_block_size_bytes);
    m_memory_size_bits = p->memory_size_bits;

    m_td_tlb_size       = p->td_tlb_size;
    m_td_tlb_latency    = p->td_tlb_latency;
    m_lcacc_tlb_size    = p->lcacc_tlb_size;
    m_lcacc_tlb_latency = p->lcacc_tlb_latency;
    m_lcacc_tlb_assoc   = p->lcacc_tlb_assoc;
    m_lcacc_tlb_mshr    = p->lcacc_tlb_mshr;

    m_dma_issue_width   = p->dma_issue_width;


#ifdef SIM_NET_PORTS
    parseAccTypes(p->acc_types);

    m_iommu     = p->iommu;
    m_tlb_hack  = p->tlb_hack;
    m_ideal_mmu = p->ideal_mmu;

    m_cryptommu_serial = p->cryptommu_serial;
    m_cryptommu_parallel = p->cryptommu_parallel;
    m_cryptommu_parallel_ra = p->cryptommu_parallel_ra;
    m_num_simics_net_ports = p->num_simics_net_ports;
    // m_num_accelerators = p->num_accelerators;
    m_num_TDs = p->num_TDs;
    m_num_acc_instances = p->num_acc_instances;
    m_host_ptw_latency = p->host_ptw_latency;
#endif

    m_warmup_enabled = false;
    m_cooldown_enabled = false;

    // Setup the global variables used in Ruby
    g_system_ptr = this;

    // Resize to the size of different machine types
    g_abs_controls.resize(MachineType_NUM);

    // Collate the statistics before they are printed.
    Stats::registerDumpCallback(new RubyStatsCallback(this));
    // Create the profiler
    m_profiler = new Profiler(p);
    m_phys_mem = p->phys_mem;
}

void
RubySystem::parseAccTypes(std::string acc_types)
{
    std::string delimiter = ",";
    size_t pos = 0;
    std::string token;
    std::vector<std::string> typeDeviceNames;
    while ((pos = acc_types.find(delimiter)) != std::string::npos) {
        token = acc_types.substr(0, pos);
        accTypes.push_back(token);
        typeDeviceNames = m_acc_dict[token];
        m_num_accelerators += typeDeviceNames.size();
        acc_types.erase(0, pos + delimiter.length());
    }

    accTypes.push_back(acc_types);
    typeDeviceNames = m_acc_dict[acc_types];
    m_num_accelerators += typeDeviceNames.size();
}

void
RubySystem::registerNetwork(Network* network_ptr)
{
  m_network = network_ptr;
}

void
RubySystem::registerAbstractController(AbstractController* cntrl)
{
  m_abs_cntrl_vec.push_back(cntrl);

  MachineID id = cntrl->getMachineID();
  g_abs_controls[id.getType()][id.getNum()] = cntrl;
}

RubySystem::~RubySystem()
{
    delete m_network;
    delete m_profiler;
}

void
RubySystem::writeCompressedTrace(uint8_t *raw_data, string filename,
                                 uint64 uncompressed_trace_size)
{
    // Create the checkpoint file for the memory
    string thefile = Checkpoint::dir() + "/" + filename.c_str();

    int fd = creat(thefile.c_str(), 0664);
    if (fd < 0) {
        perror("creat");
        fatal("Can't open memory trace file '%s'\n", filename);
    }

    gzFile compressedMemory = gzdopen(fd, "wb");
    if (compressedMemory == NULL)
        fatal("Insufficient memory to allocate compression state for %s\n",
              filename);

    if (gzwrite(compressedMemory, raw_data, uncompressed_trace_size) !=
        uncompressed_trace_size) {
        fatal("Write failed on memory trace file '%s'\n", filename);
    }

    if (gzclose(compressedMemory)) {
        fatal("Close failed on memory trace file '%s'\n", filename);
    }
    delete[] raw_data;
}

void
RubySystem::serialize(std::ostream &os)
{
    m_cooldown_enabled = true;
    vector<Sequencer*> sequencer_map;
    Sequencer* sequencer_ptr = NULL;

    for (int cntrl = 0; cntrl < m_abs_cntrl_vec.size(); cntrl++) {
        sequencer_map.push_back(m_abs_cntrl_vec[cntrl]->getSequencer());
        if (sequencer_ptr == NULL) {
            sequencer_ptr = sequencer_map[cntrl];
        }
    }

    assert(sequencer_ptr != NULL);

    for (int cntrl = 0; cntrl < m_abs_cntrl_vec.size(); cntrl++) {
        if (sequencer_map[cntrl] == NULL) {
            sequencer_map[cntrl] = sequencer_ptr;
        }
    }

    // Store the cache-block size, so we are able to restore on systems with a
    // different cache-block size. CacheRecorder depends on the correct
    // cache-block size upon unserializing.
    uint64 block_size_bytes = getBlockSizeBytes();
    SERIALIZE_SCALAR(block_size_bytes);

    DPRINTF(RubyCacheTrace, "Recording Cache Trace\n");
    // Create the CacheRecorder and record the cache trace
    m_cache_recorder = new CacheRecorder(NULL, 0, sequencer_map,
                                         block_size_bytes);

    for (int cntrl = 0; cntrl < m_abs_cntrl_vec.size(); cntrl++) {
        m_abs_cntrl_vec[cntrl]->recordCacheTrace(cntrl, m_cache_recorder);
    }

    DPRINTF(RubyCacheTrace, "Cache Trace Complete\n");
    // save the current tick value
    Tick curtick_original = curTick();
    // save the event queue head
    Event* eventq_head = eventq->replaceHead(NULL);
    DPRINTF(RubyCacheTrace, "Recording current tick %ld and event queue\n",
            curtick_original);

    // Schedule an event to start cache cooldown
    DPRINTF(RubyCacheTrace, "Starting cache flush\n");
    enqueueRubyEvent(curTick());
    simulate();
    DPRINTF(RubyCacheTrace, "Cache flush complete\n");

    // Restore eventq head
    eventq_head = eventq->replaceHead(eventq_head);
    // Restore curTick
    setCurTick(curtick_original);

    // Aggergate the trace entries together into a single array
    uint8_t *raw_data = new uint8_t[4096];
    uint64 cache_trace_size = m_cache_recorder->aggregateRecords(&raw_data,
                                                                 4096);
    string cache_trace_file = name() + ".cache.gz";
    writeCompressedTrace(raw_data, cache_trace_file, cache_trace_size);

    SERIALIZE_SCALAR(cache_trace_file);
    SERIALIZE_SCALAR(cache_trace_size);

    m_cooldown_enabled = false;
}

void
RubySystem::readCompressedTrace(string filename, uint8_t *&raw_data,
                                uint64& uncompressed_trace_size)
{
    // Read the trace file
    gzFile compressedTrace;

    // trace file
    int fd = open(filename.c_str(), O_RDONLY);
    if (fd < 0) {
        perror("open");
        fatal("Unable to open trace file %s", filename);
    }

    compressedTrace = gzdopen(fd, "rb");
    if (compressedTrace == NULL) {
        fatal("Insufficient memory to allocate compression state for %s\n",
              filename);
    }

    raw_data = new uint8_t[uncompressed_trace_size];
    if (gzread(compressedTrace, raw_data, uncompressed_trace_size) <
            uncompressed_trace_size) {
        fatal("Unable to read complete trace from file %s\n", filename);
    }

    if (gzclose(compressedTrace)) {
        fatal("Failed to close cache trace file '%s'\n", filename);
    }
}

void
RubySystem::unserialize(Checkpoint *cp, const string &section)
{
    uint8_t *uncompressed_trace = NULL;

    // This value should be set to the checkpoint-system's block-size.
    // Optional, as checkpoints without it can be run if the
    // checkpoint-system's block-size == current block-size.
    uint64 block_size_bytes = getBlockSizeBytes();
    UNSERIALIZE_OPT_SCALAR(block_size_bytes);

    string cache_trace_file;
    uint64 cache_trace_size = 0;

    UNSERIALIZE_SCALAR(cache_trace_file);
    UNSERIALIZE_SCALAR(cache_trace_size);
    cache_trace_file = cp->cptDir + "/" + cache_trace_file;

    readCompressedTrace(cache_trace_file, uncompressed_trace,
                        cache_trace_size);
    m_warmup_enabled = true;

    vector<Sequencer*> sequencer_map;
    Sequencer* t = NULL;
    for (int cntrl = 0; cntrl < m_abs_cntrl_vec.size(); cntrl++) {
        sequencer_map.push_back(m_abs_cntrl_vec[cntrl]->getSequencer());
        if (t == NULL) t = sequencer_map[cntrl];
    }

    assert(t != NULL);

    for (int cntrl = 0; cntrl < m_abs_cntrl_vec.size(); cntrl++) {
        if (sequencer_map[cntrl] == NULL) {
            sequencer_map[cntrl] = t;
        }
    }

    m_cache_recorder = new CacheRecorder(uncompressed_trace, cache_trace_size,
                                         sequencer_map, block_size_bytes);
}

#ifdef SIM_NET_PORTS
int RubySystem::threadIDtoDeviceID(int threadID)
{
#ifdef SIM_ARC
  return threadID + 1;
#endif
#ifdef SIM_TD
  return threadID + m_num_TDs;
#endif
}
int RubySystem::deviceIDtoThreadID(int deviceID)
{
#ifdef SIM_ARC
  return deviceID - 1;
#endif
#ifdef SIM_TD
  return deviceID - m_num_TDs;
#endif
}
int RubySystem::accIDtoDeviceID(int accID)
{
  std::vector<System *>::iterator system_iterator = System::systemList.begin();
  System *m5_system = *system_iterator;
  assert(m5_system);
  int num_thread_contexts = m5_system->numContexts();
#ifdef SIM_ARC
  return accID + num_thread_contexts + 1;
#endif
#ifdef SIM_TD
  return accID + num_thread_contexts + m_num_TDs;
#endif
}
int RubySystem::deviceIDtoAccID(int deviceID)
{
  std::vector<System *>::iterator system_iterator = System::systemList.begin();
  System *m5_system = *system_iterator;
  assert(m5_system);
  int num_thread_contexts = m5_system->numContexts();
#ifdef SIM_ARC
  return deviceID - num_thread_contexts - 1;
#endif
#ifdef SIM_TD
  return deviceID - num_thread_contexts - m_num_TDs;
#endif
}
int RubySystem::tdIDtoL1CacheID(int tdID)
{
  std::vector<System *>::iterator system_iterator = System::systemList.begin();
  System *m5_system = *system_iterator;
  assert(m5_system);
  int num_thread_contexts = m5_system->numContexts();
  return tdID + num_thread_contexts;
}
int RubySystem::accIDtoL1CacheID(int accID)
{
  std::vector<System *>::iterator system_iterator = System::systemList.begin();
  System *m5_system = *system_iterator;
  assert(m5_system);
  int num_thread_contexts = m5_system->numContexts();
#ifdef SIM_ARC
  return accID + num_thread_contexts + 1;
#endif
#ifdef SIM_TD
  return accID + num_thread_contexts + m_num_TDs;
#endif
}
#define CENTER_ID 12
int RubySystem::mapPortID(int port)
{
  if(port >= CENTER_ID)
    port = (port + 1) % m_num_simics_net_ports;
  return port;
}
#endif

void
RubySystem::startup()
{

    // Ruby restores state from a checkpoint by resetting the clock to 0 and
    // playing the requests that can possibly re-generate the cache state.
    // The clock value is set to the actual checkpointed value once all the
    // requests have been executed.
    //
    // This way of restoring state is pretty finicky. For example, if a
    // Ruby component reads time before the state has been restored, it would
    // cache this value and hence its clock would not be reset to 0, when
    // Ruby resets the global clock. This can potentially result in a
    // deadlock.
    //
    // The solution is that no Ruby component should read time before the
    // simulation starts. And then one also needs to hope that the time
    // Ruby finishes restoring the state is less than the time when the
    // state was checkpointed.

#ifdef SIM_NET_PORTS
    /* malloc and initialize light weight interrupt interface */
    g_lwInt_interface = (lwInt_ifc_t *)malloc(sizeof(lwInt_ifc_t));
    memset(g_lwInt_interface, 0, sizeof(lwInt_ifc_t));
    warn("REGISTERING INTERRUPT INTERFACE %p\n", g_lwInt_interface);
    g_lwInt_interface->raiseLightWeightInt = &raiseLightWeightInt;
    g_lwInt_interface->isReady = &isReady;

    g_networkPort_interface = (SimicsNetworkPortInterface *)malloc(sizeof(SimicsNetworkPortInterface));
    memset(g_networkPort_interface, 0, sizeof(SimicsNetworkPortInterface));
    InitializeSimicsNetworkPortInterface(g_networkPort_interface);
    warn("Initialize SimicsNetworkPortInterface\n");//*/

    //Zhenman: Currently we initialize one networkInterrupts handler per ThreadContext
    int i;

    //suppose we just have one M5 system running
    //System::printSystems();
    std::vector<System *>::iterator system_iterator = System::systemList.begin();
    System *m5_system = *system_iterator;
    //assert(m5_system);
    int num_thread_contexts = m5_system->numContexts();
    g_network_interrupt_handle.resize(num_thread_contexts);
    for(i = 0; i < g_network_interrupt_handle.size(); i++) {
      g_network_interrupt_handle[i] = createNetworkInterruptHandle(i, RubySystem::threadIDtoDeviceID(i), i);
    }
    warn("Initialize network_interrupt_handle\n");//*/
     if(m_num_acc_instances==8*4) /*ports not included*/
    {
      m_num_acc_instances--;
      m_num_acc_instances--;
    }
    g_network_interrupts.resize(num_thread_contexts);
    for(i = 0; i < g_network_interrupts.size(); i++) {
      g_network_interrupts[i] = new NetworkInterrupts(g_network_interrupt_handle[i]);
    }
    warn("Initialize network_interrupts\n");//*/

    g_LCAccDeviceHandle.resize(m_num_accelerators * m_num_acc_instances);
    for(i = 0; i < g_LCAccDeviceHandle.size(); i++) {
      g_LCAccDeviceHandle[i] = LCAcc::CreateNewLCAccDeviceHandle();
    }//*/

    //Let's first test non AIM
    g_memObject = NULL;
    g_memInterface = NULL;

    g_dmaDevice.resize(m_num_accelerators * m_num_acc_instances);
    for(i = 0; i < g_dmaDevice.size(); i++) {
      g_dmaDevice[i] = CreateDMAEngineHandle();
    }
    warn("Initialize dmaDevice vector\n");

    g_dmaInterface = CreateDMAEngineInterface();
    warn("Initialize dmaInterface\n");

    g_spmObject.resize(m_num_accelerators * m_num_acc_instances);
    for(i = 0; i < g_spmObject.size(); i++) {
      g_spmObject[i] = CreateNewScratchPad();
    }
    warn("Initialize spmObject vector\n");

    g_spmInterface = CreateScratchPadInterface();
    warn("Initialize spmInterface\n");

    g_LCAccInterface = LCAcc::CreateLCAccInterface();
    warn("Initialize LCAccInterface\n");//*/
    int portID = RubySystem::accIDtoDeviceID(0);

    int TD_portID = num_thread_contexts;

#ifdef SIM_ARC
    //ARC version, uses GAM
    int GAM_portID = num_thread_contexts;
    g_gamObject = CreateGAMInterface();
    warn("Initialize GAMInterface\n");
    g_gamObject->SetNetPort(GAM_portID, 0);
    g_gamObject->SetActivityDelay(0);
    g_gamObject->SetIntervalLength(10000);
    g_gamObject->SetMinAllocCount(1);
#endif

#ifdef SIM_TD
    //uses TD
    g_TdDmaDevice.resize(m_num_TDs);
    for(i = 0; i < g_TdDmaDevice.size(); i++) {
      g_TdDmaDevice[i] = CreateDMAEngineHandle();
    }
    warn("Initialize TdDmaDevice vector\n");

    g_TdSpmObject.resize(m_num_TDs);
    for(i = 0; i < g_TdSpmObject.size(); i++) {
      g_TdSpmObject[i] = CreateNewScratchPad();
    }
    warn("Initialize TdSpmObject vector\n");

    g_TDInterface = CreateTDInterface();
#endif

#ifdef SIM_ARC_PLUS
    //ARC PLUS version
    g_TDHandle.resize(m_num_TDs);
    for(i = 0; i < g_TDHandle.size(); i++) {
      g_TDHandle[i] = CreateNewTDHandle();
      g_TDInterface->SetNetworkPort(g_TDHandle[i],TD_portID+i,i);
      g_TDInterface->SetSelector(g_TDHandle[i],"MinChain");
      g_TDInterface->SetTaskGrain(g_TDHandle[i],8);
      g_TDInterface->SetCFUAllocationPerIsland(g_TDHandle[i],999999);
      g_TDInterface->Initialize(g_TDHandle[i]);
    }
    // warn("Initialize TDInterface\n");
#endif

#ifdef SIM_CHARM
    //CHARM version
    g_TDHandle.resize(m_num_TDs);
    for(i = 0; i < g_TDHandle.size(); i++) {
      g_TDHandle[i] = CreateNewTDHandle();
      g_TDInterface->SetNetworkPort(g_TDHandle[i],TD_portID+i,i);
      //g_TDInterface->SetSelector(g_TDHandle[i],"SmoothMemory");
      g_TDInterface->SetSelector(g_TDHandle[i],"MinChain");
      g_TDInterface->SetTaskGrain(g_TDHandle[i],8);
      g_TDInterface->SetCFUAllocationPerIsland(g_TDHandle[i],999999);
      g_TDInterface->SetFpgaArea(g_TDHandle[i],0);
      g_TDInterface->Initialize(g_TDHandle[i]);
    }
    warn("Initialize TDInterface\n");
#endif

    //This is the number of applications to run, to synchronize
    //SetBarrierWidth(1);

#if defined(SIM_ARC) || defined(SIM_ARC_PLUS) && !defined(SIM_SW)

    int id = 0;
    for (i = 0; i < accTypes.size(); i++) {
        std::string bench = accTypes[i];
        std::vector<std::string> subAccNames = m_acc_dict[bench];
        int numSubAccs = subAccNames.size();

        for (int j = 0; j < m_num_acc_instances; j++) {
            for (int k = 0; k < numSubAccs; k++) {
                g_LCAccInterface->SetNetPort(
                    g_LCAccDeviceHandle[id],
                    portID + j + i * m_num_acc_instances,
                    RubySystem::accIDtoDeviceID(id));

                g_LCAccInterface->AddOperatingMode(
                    g_LCAccDeviceHandle[id],
                    subAccNames[k].c_str());

                id++;
            }
        }
    }

    //All benchmarks except for Texture_Synthesis, actually bind to one port
    //comment this out for Texture_Synthesis case
    // for (int k = 0; k < m_num_acc_instances; k++) {
    //     for (i = 0; i < m_num_accelerators; i++) {
    //         g_LCAccInterface->SetNetPort(
    //             g_LCAccDeviceHandle[i + k * m_num_accelerators], portID + k,
    //             RubySystem::accIDtoDeviceID(i + k * m_num_accelerators));
    //     }

        //vector add
        //g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0], "VectorAddSample");

        /*###Medical Imaging(MI) benchmarks##############*/
        // //deblur
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "denoise1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "blur1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "deblur1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[3+k*m_num_accelerators], "deblur2Mega");

        //denoise
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "denoise1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "denoise2Mega");

        // Denoise
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[index++],
        //     "denoise1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[index++],
        //     "denoise2Mega");

        //registration
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "registration1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "blur1Mega");

        //segmentation
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "segmentation1Mega");

        //CS_MR: ignore this benchmark
        /*g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "CSMR_stage1");
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "CSMR_stage2");
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "CSMR_stage3_1");
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[3+k*m_num_accelerators], "CSMR_stage3_2");
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[4+k*m_num_accelerators], "CSMR_stage_fft");
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[5+k*m_num_accelerators], "CSMR_stage_fft");//*/
        /*###Medical Imaging(MI) benchmarks##############*/

        /*###CoMmercial(CM) benchmarks###################*/
        //blackScholes
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[k], "blackScholes");

        // BlackScholes
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[index++],
        //     "blackScholes");

        //streamCluster
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "streamCluster1");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "streamCluster3");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "streamCluster4");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[3+k*m_num_accelerators], "streamCluster5");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[4+k*m_num_accelerators], "streamCluster6");//

        // //swaptions
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "swaptions1");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "swaptions2");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "swaptions3");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[3+k*m_num_accelerators], "swaptions4");//
        /*###CoMmercial(CM) benchmarks#######   ############*/

        /*###computer VISion(VIS) benchmarks#############*/
        //LPCIP
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[k*m_num_accelerators], "LPCIP");

        // //SURF
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "SURF1Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "SURF2Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "SURF3Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[3+k*m_num_accelerators], "SURF4Mega");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[4+k*m_num_accelerators], "SURF5Mega");//

        //Texture_Synthesis: still suffer from page fault
        /*g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[0], portID, RubySystem::accIDtoDeviceID(0));
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0], "TexSynth1");
        for(i = 1; i < 10; i++) {
          g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[i], portID+i-1, RubySystem::accIDtoDeviceID(i));
          g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "TexSynth2");
        }
        for(i = 10; i < 14; i++) {
          g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[i], portID+i-10, RubySystem::accIDtoDeviceID(i));
          g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "TexSynth3");
        }
        for(i = 14; i < 18; i++) {
          g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[i], portID+i-14, RubySystem::accIDtoDeviceID(i));
          g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "TexSynth4");
        }
        g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[18], portID, RubySystem::accIDtoDeviceID(18));
        g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[18], "TexSynth5");//*/
        /*###computer VISion(VIS) benchmarks#############*/

        /*###NAVigation(NAV) benchmarks##################*/
        //robot_localization
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[k*m_num_accelerators], "RobLoc");

        //disparity_map
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "DispMapCompSAD");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "DispMapFindDisp");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[2+k*m_num_accelerators], "DispMapIntegSum");//

        //ekf_slam
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[0+k*m_num_accelerators], "Jacobians");
        // g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[1+k*m_num_accelerators], "SphericalCoords");//
        /*###NAVigation(NAV) benchmarks##################*/
    // }

    for(i = 0; i < m_num_accelerators * m_num_acc_instances; i++) {
      g_LCAccInterface->SetPrefetchDistance(g_LCAccDeviceHandle[i], 0);
      g_LCAccInterface->SetSPMConfig(g_LCAccDeviceHandle[i], 1, 2048, 2, 1, 1, 1);
      //g_LCAccInterface->SetSPMConfig(g_LCAccDeviceHandle[i], 1, 64, 100000, 0, 100000, 0); //ideal SPM
      g_LCAccInterface->Initialize(g_LCAccDeviceHandle[i], 0);
    }//*/
#endif

#if defined(SIM_CHARM)
    for(i = 0; i < m_num_accelerators; i++) {
      g_LCAccInterface->SetNetPort(g_LCAccDeviceHandle[i], portID+i%24, RubySystem::accIDtoDeviceID(i));
    }

    /*###Medical Imaging(MI) benchmarks##############*/
    //deblur, denoise, registration segmentation
    assert(m_num_accelerators == 176);
    for(i = 0; i < 124; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "poly");
    }
    for(; i < 156; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "divide");
    }
    for(; i < 168; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sqrtf");
    }
    for(; i < 176; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "powf");
    } //*/
    /*###Medical Imaging(MI) benchmarks##############*/

    /*###CoMmercial(CM) benchmarks###################*/
    //blackScholes, streamCluster, swaptions
    /*assert(m_num_accelerators == 166);
    for(i = 0; i < 102; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "poly");
    }
    for(; i < 114; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "divide");
    }
    for(; i < 117; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sqrtf");
    }
    for(; i < 121; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "logf");
    }
    for(; i < 130; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "powf");
    }
    for(; i < 136; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectfi");
    }
    for(; i < 145; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectif");
    }
    for(; i < 151; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectff");
    }
    for(; i < 154; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectii");
    }
    for(; i < 160; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sample1");
    }
    for(; i < 163; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sampleWrite1");
    }
    for(; i < 166; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "Sum");
    } //*/
    /*###CoMmercial(CM) benchmarks###################*/

    /*###computer VISion(VIS) benchmarks#############*/
    //LPCIP, SURF, TextureSynthesis
    /*assert(m_num_accelerators == 747);
    for(i = 0; i < 343; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "poly");
    }
    for(; i < 373; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "divide");
    }
    for(; i < 378; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sqrtf");
    }
    for(; i < 387; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "powf");
    }
    for(; i < 394; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectfi");
    }
    for(; i < 408; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectif");
    }
    for(; i < 476; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "selectff");
    }
    for(; i < 490; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "Mod");
    }
    for(; i < 495; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "castif");
    }
    for(; i < 557; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "castfi");
    }
    for(; i < 566; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sample2");
    }
    for(; i < 692; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sample3");
    }
    for(; i < 732; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "Sum");
    }
    for(; i < 747; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "Min");
    }//*/
    /*###computer VISion(VIS) benchmarks#############*/

    /*###NAVigation(NAV) benchmarks##################*/
    //RobotLocalization, DisparityMap, EKF_SLAM
    /*assert(m_num_accelerators == 149);
    for(i = 0; i < 115; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "poly");
    }
    for(; i < 122; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "divide");
    }
    for(; i < 125; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "sqrtf");
    }
    for(; i < 149; i++) {
      g_LCAccInterface->AddOperatingMode(g_LCAccDeviceHandle[i], "Sum");
    } //*/
    /*###NAVigation(NAV) benchmarks##################*/

    for(i = 0; i < m_num_accelerators; i++) {
      g_LCAccInterface->SetPrefetchDistance(g_LCAccDeviceHandle[i], 0);
      g_LCAccInterface->SetSPMConfig(g_LCAccDeviceHandle[i], 1, 2048, 2, 1, 1, 1);
      g_LCAccInterface->Initialize(g_LCAccDeviceHandle[i], 0);
    }
#endif

    for(int i = 0; i < 30000; i++) {
      cycleCBRing.push_back(std::queue<CBContainer>());
    }

#endif

    if (m_warmup_enabled) {
        // save the current tick value
        Tick curtick_original = curTick();
        // save the event queue head
        Event* eventq_head = eventq->replaceHead(NULL);
        // set curTick to 0 and reset Ruby System's clock
        setCurTick(0);
        resetClock();

        // Schedule an event to start cache warmup
        enqueueRubyEvent(curTick());
        simulate();

        delete m_cache_recorder;
        m_cache_recorder = NULL;
        m_warmup_enabled = false;

        // Restore eventq head
        eventq_head = eventq->replaceHead(eventq_head);
        // Restore curTick and Ruby System's clock
        setCurTick(curtick_original);
        resetClock();
    }

    resetStats();
}

void
RubySystem::RubyEvent::process()
{
    if (ruby_system->m_warmup_enabled) {
        ruby_system->m_cache_recorder->enqueueNextFetchRequest();
    }  else if (ruby_system->m_cooldown_enabled) {
        ruby_system->m_cache_recorder->enqueueNextFlushRequest();
    }
}

void
RubySystem::resetStats()
{
    g_ruby_start = curCycle();
}

bool
RubySystem::functionalRead(PacketPtr pkt)
{
    Address address(pkt->getAddr());
    Address line_address(address);
    line_address.makeLineAddress();

    AccessPermission access_perm = AccessPermission_NotPresent;
    int num_controllers = m_abs_cntrl_vec.size();

    DPRINTF(RubySystem, "Functional Read request for %s\n",address);

    unsigned int num_ro = 0;
    unsigned int num_rw = 0;
    unsigned int num_busy = 0;
    unsigned int num_backing_store = 0;
    unsigned int num_invalid = 0;

    #ifdef SIM_NET_PORTS
    std::vector<System *>::iterator system_iterator = System::systemList.begin();
    System *m5_system = *system_iterator;
    assert(m5_system);
    int num_thread_contexts = m5_system->numContexts();
    #endif

    // In this loop we count the number of controllers that have the given
    // address in read only, read write and busy states.
    for (unsigned int i = 0; i < num_controllers; ++i) {
        access_perm = m_abs_cntrl_vec[i]-> getAccessPermission(line_address);
        if (access_perm == AccessPermission_Read_Only)
            num_ro++;
        else if (access_perm == AccessPermission_Read_Write)
            num_rw++;
        else if (access_perm == AccessPermission_Busy)
            num_busy++;
        else if (access_perm == AccessPermission_Backing_Store)
            // See RubySlicc_Exports.sm for details, but Backing_Store is meant
            // to represent blocks in memory *for Broadcast/Snooping protocols*,
            // where memory has no idea whether it has an exclusive copy of data
            // or not.
            num_backing_store++;
        else if (access_perm == AccessPermission_Invalid ||
                 access_perm == AccessPermission_NotPresent)
            num_invalid++;
    }
    assert(num_rw <= 1);

    // This if case is meant to capture what happens in a Broadcast/Snoop
    // protocol where the block does not exist in the cache hierarchy. You
    // only want to read from the Backing_Store memory if there is no copy in
    // the cache hierarchy, otherwise you want to try to read the RO or RW
    // copies existing in the cache hierarchy (covered by the else statement).
    // The reason is because the Backing_Store memory could easily be stale, if
    // there are copies floating around the cache hierarchy, so you want to read
    // it only if it's not in the cache hierarchy at all.
    if (num_invalid == (num_controllers - 1) && num_backing_store == 1) {
        DPRINTF(RubySystem, "only copy in Backing_Store memory, read from it\n");
        for (unsigned int i = 0; i < num_controllers; ++i) {
            access_perm = m_abs_cntrl_vec[i]->getAccessPermission(line_address);
            if (access_perm == AccessPermission_Backing_Store) {
	        #ifdef SIM_NET_PORTS
 	        MachineID id = m_abs_cntrl_vec[i]->getMachineID();
	        if( !( (id.getType() == MachineType_L1Cache) && (id.getNum() >= (num_thread_contexts+m_num_accelerators*m_num_acc_instances+m_num_TDs)) ) )
		  m_abs_cntrl_vec[i]->functionalRead(line_address, pkt);
                #else
	        m_abs_cntrl_vec[i]->functionalRead(line_address, pkt);
                #endif
                return true;
            }
        }
    } else if (num_ro > 0 || num_rw == 1) {
        // In Broadcast/Snoop protocols, this covers if you know the block
        // exists somewhere in the caching hierarchy, then you want to read any
        // valid RO or RW block.  In directory protocols, same thing, you want
        // to read any valid readable copy of the block.
        DPRINTF(RubySystem, "num_busy = %d, num_ro = %d, num_rw = %d\n",
                num_busy, num_ro, num_rw);
        // In this loop, we try to figure which controller has a read only or
        // a read write copy of the given address. Any valid copy would suffice
        // for a functional read.
        for (unsigned int i = 0;i < num_controllers;++i) {
            access_perm = m_abs_cntrl_vec[i]->getAccessPermission(line_address);
            if (access_perm == AccessPermission_Read_Only ||
                access_perm == AccessPermission_Read_Write) {
	        #ifdef SIM_NET_PORTS
 	        MachineID id = m_abs_cntrl_vec[i]->getMachineID();
	        if( !( (id.getType() == MachineType_L1Cache) && (id.getNum() >= (num_thread_contexts+m_num_accelerators*m_num_acc_instances+m_num_TDs)) ) )
		  m_abs_cntrl_vec[i]->functionalRead(line_address, pkt);
                #else
	        m_abs_cntrl_vec[i]->functionalRead(line_address, pkt);
                #endif
                return true;
            }
        }
    }

    return false;
}

// The function searches through all the buffers that exist in different
// cache, directory and memory controllers, and in the network components
// and writes the data portion of those that hold the address specified
// in the packet.
bool
RubySystem::functionalWrite(PacketPtr pkt)
{
    Address addr(pkt->getAddr());
    Address line_addr = line_address(addr);
    AccessPermission access_perm = AccessPermission_NotPresent;
    int num_controllers = m_abs_cntrl_vec.size();

    DPRINTF(RubySystem, "Functional Write request for %s\n",addr);

    uint32_t M5_VAR_USED num_functional_writes = 0;

#ifdef SIM_NET_PORTS
    std::vector<System *>::iterator system_iterator = System::systemList.begin();
    System *m5_system = *system_iterator;
    assert(m5_system);
    int num_thread_contexts = m5_system->numContexts();
#endif

    for (unsigned int i = 0; i < num_controllers;++i) {
#ifdef SIM_NET_PORTS
      MachineID id = m_abs_cntrl_vec[i]->getMachineID();
      if( !( (id.getType() == MachineType_L1Cache) && (id.getNum() >= (num_thread_contexts+m_num_accelerators*m_num_acc_instances+m_num_TDs)) ) )
        num_functional_writes +=
            m_abs_cntrl_vec[i]->functionalWriteBuffers(pkt);
#else
        num_functional_writes +=
            m_abs_cntrl_vec[i]->functionalWriteBuffers(pkt);
#endif

        access_perm = m_abs_cntrl_vec[i]->getAccessPermission(line_addr);
        if (access_perm != AccessPermission_Invalid &&
            access_perm != AccessPermission_NotPresent) {
#ifdef SIM_NET_PORTS
	  MachineID id = m_abs_cntrl_vec[i]->getMachineID();
	  if( !( (id.getType() == MachineType_L1Cache) && (id.getNum() >= (num_thread_contexts+m_num_accelerators*m_num_acc_instances+m_num_TDs)) ) )
            num_functional_writes +=
                m_abs_cntrl_vec[i]->functionalWrite(line_addr, pkt);
#else
            num_functional_writes +=
                m_abs_cntrl_vec[i]->functionalWrite(line_addr, pkt);
#endif
        }
    }

    num_functional_writes += m_network->functionalWrite(pkt);
    DPRINTF(RubySystem, "Messages written = %u\n", num_functional_writes);

    return true;
}

#ifdef CHECK_COHERENCE
// This code will check for cases if the given cache block is exclusive in
// one node and shared in another-- a coherence violation
//
// To use, the SLICC specification must call sequencer.checkCoherence(address)
// when the controller changes to a state with new permissions.  Do this
// in setState.  The SLICC spec must also define methods "isBlockShared"
// and "isBlockExclusive" that are specific to that protocol
//
void
RubySystem::checkGlobalCoherenceInvariant(const Address& addr)
{
#if 0
    NodeID exclusive = -1;
    bool sharedDetected = false;
    NodeID lastShared = -1;

    for (int i = 0; i < m_chip_vector.size(); i++) {
        if (m_chip_vector[i]->isBlockExclusive(addr)) {
            if (exclusive != -1) {
                // coherence violation
                WARN_EXPR(exclusive);
                WARN_EXPR(m_chip_vector[i]->getID());
                WARN_EXPR(addr);
                WARN_EXPR(getTime());
                ERROR_MSG("Coherence Violation Detected -- 2 exclusive chips");
            } else if (sharedDetected) {
                WARN_EXPR(lastShared);
                WARN_EXPR(m_chip_vector[i]->getID());
                WARN_EXPR(addr);
                WARN_EXPR(getTime());
                ERROR_MSG("Coherence Violation Detected -- exclusive chip with >=1 shared");
            } else {
                exclusive = m_chip_vector[i]->getID();
            }
        } else if (m_chip_vector[i]->isBlockShared(addr)) {
            sharedDetected = true;
            lastShared = m_chip_vector[i]->getID();

            if (exclusive != -1) {
                WARN_EXPR(lastShared);
                WARN_EXPR(exclusive);
                WARN_EXPR(addr);
                WARN_EXPR(getTime());
                ERROR_MSG("Coherence Violation Detected -- exclusive chip with >=1 shared");
            }
        }
    }
#endif
}
#endif

RubySystem *
RubySystemParams::create()
{
    return new RubySystem(this);
}

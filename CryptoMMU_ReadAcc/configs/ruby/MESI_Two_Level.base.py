# Copyright (c) 2006-2007 The Regents of The University of Michigan
# Copyright (c) 2009 Advanced Micro Devices, Inc.
# All rights reserved.
#
# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are
# met: redistributions of source code must retain the above copyright
# notice, this list of conditions and the following disclaimer;
# redistributions in binary form must reproduce the above copyright
# notice, this list of conditions and the following disclaimer in the
# documentation and/or other materials provided with the distribution;
# neither the name of the copyright holders nor the names of its
# contributors may be used to endorse or promote products derived from
# this software without specific prior written permission.
#
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
# "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
# LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
# A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
# OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
# SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
# LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
# DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
# THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
# (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
# OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# Authors: Brad Beckmann

import math
import m5
from m5.objects import *
from m5.defines import buildEnv
from Ruby import create_topology
from Ruby import send_evicts
from m5.util import addToPath, fatal

addToPath('../lcacc')
import Lcacc

#
# Note: the L1 Cache latency is only used by the sequencer on fast path hits
#
class L1Cache(RubyCache):
    latency = 3

#
# Note: the L2 Cache latency is not currently used
#
class L2Cache(RubyCache):
    latency = 15

def define_options(parser):
    return

def create_system(options, full_system, system, dma_ports, ruby_system):

    if buildEnv['PROTOCOL'] != 'MESI_Two_Level':
        fatal("This script requires the MESI_Two_Level protocol to be built.")

    ruby_system.num_simics_net_ports = options.num_networkports
    ruby_system.num_accelerators = options.accelerators
    ruby_system.num_TDs = options.num_tds
    ruby_system.num_acc_instances = options.num_accinstances
    ruby_system.td_tlb_size = options.td_tlb_size
    ruby_system.lcacc_tlb_size = options.lcacc_tlb_size

    cpu_sequencers = []

    #
    # The ruby network creation expects the list of nodes in the system to be
    # consistent with the NetDest list.  Therefore the l1 controller nodes must be
    # listed before the directory nodes and directory nodes before dma nodes, etc.
    #
    netport_cntrl_nodes = []
    l1_cntrl_nodes = []
    l2_cntrl_nodes = []
    dir_cntrl_nodes = []
    dma_cntrl_nodes = []

    #
    # Must create the individual controllers before the network to ensure the
    # controller constructors are called before the network constructor
    #
    l2_bits = int(math.log(options.num_l2caches, 2))
    block_size_bits = int(math.log(options.cacheline_size, 2))
    l2_cache_sets = MemorySize(options.l2_size)/options.cacheline_size/options.l2_assoc
    l2_cache_set_bits = int(math.log(l2_cache_sets,2))
    print "l2_cache_set_bits = %d" % l2_cache_set_bits

    assert(options.num_networkports == options.num_l2caches)
    num_l1_cntrls = ((options.accelerators + options.num_tds + options.num_networkports - 1)/options.num_networkports) * options.num_networkports
    print "num_l1_cntrls = %d" % num_l1_cntrls
    assert(num_l1_cntrls >= (options.accelerators + options.num_tds))

    for i in xrange(options.num_networkports):
        # First create the Ruby objects associated with
        # the CPU and Accelerator signal communication
        netport_cntrl = SimicsNetworkPortInterface_Controller(version = i,
                        transitions_per_cycle=options.ports,
                        ruby_system = ruby_system)

        exec("ruby_system.netport_cntrl%d = netport_cntrl" % i)
        netport_cntrl_nodes.append(netport_cntrl)
        # Connect the netport controller to the network
        netport_cntrl.messageOut = ruby_system.network.slave
        netport_cntrl.messageIn = ruby_system.network.master

    for i in xrange(num_l1_cntrls):
        #
        # First create the Ruby objects associated with this cpu
        #
        l1i_cache = L1Cache(size = options.l1i_size,
                            assoc = options.l1i_assoc,
                            start_index_bit = block_size_bits,
                            is_icache = True)
        l1d_cache = L1Cache(size = options.l1d_size,
                            assoc = options.l1d_assoc,
                            start_index_bit = block_size_bits,
                            is_icache = False)

        prefetcher = RubyPrefetcher.Prefetcher()

        l1_cntrl = L1Cache_Controller(version = i,
                                      L1Icache = l1i_cache,
                                      L1Dcache = l1d_cache,
                                      l2_select_num_bits = l2_bits,
				      l2_select_low_bit = (block_size_bits + l2_cache_set_bits),
                                      send_evictions = send_evicts(options),
                                      prefetcher = prefetcher,
                                      ruby_system = ruby_system,
                                      clk_domain=system.cpu[0].clk_domain,
                                      transitions_per_cycle=options.ports,
                                      enable_prefetch = False)

        cpu_seq = RubySequencer(version = i,
                                icache = l1i_cache,
                                dcache = l1d_cache,
                                clk_domain=system.cpu[0].clk_domain,
                                ruby_system = ruby_system)

        l1_cntrl.sequencer = cpu_seq
        exec("ruby_system.l1_cntrl%d = l1_cntrl" % i)

        # Add controllers and sequencers to the appropriate lists
	if len(cpu_sequencers) < options.num_cpus :
            cpu_sequencers.append(cpu_seq)
        l1_cntrl_nodes.append(l1_cntrl)

        # Connect the L1 controllers and the network
        l1_cntrl.requestFromL1Cache =  ruby_system.network.slave
        l1_cntrl.responseFromL1Cache =  ruby_system.network.slave
        l1_cntrl.unblockFromL1Cache =  ruby_system.network.slave

        l1_cntrl.requestToL1Cache =  ruby_system.network.master
        l1_cntrl.responseToL1Cache =  ruby_system.network.master


    l2_index_start = block_size_bits

    for i in xrange(options.num_l2caches):
        #
        # First create the Ruby objects associated with this cpu
        #
        l2_cache = L2Cache(size = options.l2_size,
                           assoc = options.l2_assoc,
                           start_index_bit = l2_index_start)

        l2_cntrl = L2Cache_Controller(version = i,
                                      L2cache = l2_cache,
                                      transitions_per_cycle=options.ports,
                                      ruby_system = ruby_system)

        exec("ruby_system.l2_cntrl%d = l2_cntrl" % i)
        l2_cntrl_nodes.append(l2_cntrl)

        # Connect the L2 controllers and the network
        l2_cntrl.DirRequestFromL2Cache = ruby_system.network.slave
        l2_cntrl.L1RequestFromL2Cache = ruby_system.network.slave
        l2_cntrl.responseFromL2Cache = ruby_system.network.slave

        l2_cntrl.unblockToL2Cache = ruby_system.network.master
        l2_cntrl.L1RequestToL2Cache = ruby_system.network.master
        l2_cntrl.responseToL2Cache = ruby_system.network.master


    phys_mem_size = sum(map(lambda r: r.size(), system.mem_ranges))
    assert(phys_mem_size % options.num_dirs == 0)
    mem_module_size = phys_mem_size / options.num_dirs


    # Run each of the ruby memory controllers at a ratio of the frequency of
    # the ruby system
    # clk_divider value is a fix to pass regression.
    ruby_system.memctrl_clk_domain = DerivedClockDomain(
                                          clk_domain=ruby_system.clk_domain,
                                          clk_divider=3)

    for i in xrange(options.num_dirs):
        dir_size = MemorySize('0B')
        dir_size.value = mem_module_size

        dir_cntrl = Directory_Controller(version = i,
                                         directory = RubyDirectoryMemory(
                                             version = i, size = dir_size),
                                         transitions_per_cycle = options.ports,
                                         ruby_system = ruby_system)

        exec("ruby_system.dir_cntrl%d = dir_cntrl" % i)
        dir_cntrl_nodes.append(dir_cntrl)

        # Connect the directory controllers and the network
        dir_cntrl.requestToDir = ruby_system.network.master
        dir_cntrl.responseToDir = ruby_system.network.master
        dir_cntrl.responseFromDir = ruby_system.network.slave

    for i, dma_port in enumerate(dma_ports):
        # Create the Ruby objects associated with the dma controller
        dma_seq = DMASequencer(version = i,
                               ruby_system = ruby_system,
                               slave = dma_port)

        dma_cntrl = DMA_Controller(version = i,
                                   dma_sequencer = dma_seq,
                                   transitions_per_cycle = options.ports,
                                   ruby_system = ruby_system)

        exec("ruby_system.dma_cntrl%d = dma_cntrl" % i)
        dma_cntrl_nodes.append(dma_cntrl)

        # Connect the dma controller to the network
        dma_cntrl.responseFromDir = ruby_system.network.master
        dma_cntrl.requestToDir = ruby_system.network.slave

    all_cntrls = netport_cntrl_nodes + \
                 l1_cntrl_nodes + \
                 l2_cntrl_nodes + \
                 dir_cntrl_nodes + \
                 dma_cntrl_nodes

    # Create the io controller and the sequencer
    if full_system:
        io_seq = DMASequencer(version=len(dma_ports), ruby_system=ruby_system)
        ruby_system._io_port = io_seq
        io_controller = DMA_Controller(version = len(dma_ports),
                                       dma_sequencer = io_seq,
                                       ruby_system = ruby_system)
        ruby_system.io_controller = io_controller

        # Connect the dma controller to the network
        io_controller.responseFromDir = ruby_system.network.master
        io_controller.requestToDir = ruby_system.network.slave

        all_cntrls = all_cntrls + [io_controller]

    topology = create_topology(all_cntrls, options)
    return (cpu_sequencers, dir_cntrl_nodes, topology)

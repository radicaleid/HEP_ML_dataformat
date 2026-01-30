"""
ROOT NTuple Reader for GNN4ITk data using RDataFrame
Enhanced version with explicit RNTuple support via RDF.RNTupleDS
Ensures proper RNTuple handling throughout the processing pipeline
"""

#!/usr/bin/env python3
import ROOT
import numpy as np
import pandas as pd
import time
import psutil
import os, sys, glob
from typing import Tuple, Dict, List
import argparse
from contextlib import contextmanager

start_time = time.time()
print(f"[{time.time()-start_time:.2f}] Imported necessary modules.")

@contextmanager
def root_multithreading(num_threads: int = 1, keep_enable: bool = True):
    """
    Context manager to temporarily enable/disable ROOT ImplicitMT for RDataFrame operations.
    """
    was_enabled = ROOT.IsImplicitMTEnabled()
    original_threads = ROOT.GetThreadPoolSize() if was_enabled else 0
    
    try:
        if num_threads > 1:
            if not was_enabled:
                ROOT.EnableImplicitMT(num_threads)
                print(f"✓ ROOT ImplicitMT enabled with {num_threads} threads")
            elif original_threads != num_threads:
                ROOT.DisableImplicitMT()
                ROOT.EnableImplicitMT(num_threads)
                print(f"✓ ROOT ImplicitMT adjusted to {num_threads} threads")
        else:
            if was_enabled:
                ROOT.DisableImplicitMT()
                print(f"✓ ROOT ImplicitMT disabled (single-threaded mode)")
        
        yield
        
    finally:
        current_enabled = ROOT.IsImplicitMTEnabled()
        if not keep_enable and current_enabled:
            ROOT.DisableImplicitMT()
        print(f"ROOT_MAX_THREADS: {ROOT.GetThreadPoolSize()}")
        print(f"ROOT ImplicitMT is {'enabled' if ROOT.IsImplicitMTEnabled() else 'disabled'}")


class PerformanceMonitor:
    """Simple performance monitoring class"""
    
    def __init__(self):
        self.start_time = time.time()
        self.process = psutil.Process(os.getpid())
        self.initial_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.last_time = self.start_time

    def print_stats(self, operation: str):
        """Print timing, process memory, and system-wide available memory"""
        current_time = time.time()
        current_memory = self.process.memory_info().rss / 1024 / 1024  # Process RSS in MB
        system_memory = psutil.virtual_memory()
        available_gb = system_memory.available / (1024**3) # System available in GB
        
        elapsed = current_time - self.last_time
        memory_delta = current_memory - self.initial_memory
        
        print(f"\nROOT # of threads is {ROOT.GetThreadPoolSize()}")
        print(f"=== Performance Report: {operation} ===")
        print(f"Time: {elapsed:.2f} seconds")
        print(f"Memory Delta: {memory_delta:.1f} MB")
        print(f"Process Total Memory: {current_memory:.1f} MB")
        print(f"System Available Memory: {available_gb:.2f} GB")
        print("=" * 50 + "\n")
        
        self.last_time = current_time
        self.initial_memory = current_memory


# C++ code for optimized operations
cpp_code = '''
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>

// Cluster matching with hash map for O(1) lookups + barcode padding
std::tuple<
    ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>
>
ClusterMatching(const ROOT::VecOps::RVec<int>& hit_cl1_idx, 
                const ROOT::VecOps::RVec<int>& hit_cl2_idx,
                const ROOT::VecOps::RVec<int>& cluster_indices,
                const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& CLparticleLink_barcode) {
    
    std::unordered_map<int, size_t> cluster_map;
    cluster_map.reserve(cluster_indices.size());
    for (size_t i = 0; i < cluster_indices.size(); ++i) {
        cluster_map[cluster_indices[i]] = i;
    }
    
    ROOT::VecOps::RVec<int> matched_cl1_pos, matched_cl2_pos, valid_hit_idx;
    ROOT::VecOps::RVec<int> matched_barcodes_1, matched_barcodes_2;
    
    for (size_t hit_idx = 0; hit_idx < hit_cl1_idx.size(); ++hit_idx) {
        auto cl1_it = cluster_map.find(hit_cl1_idx[hit_idx]);
        auto cl2_it = cluster_map.find(hit_cl2_idx[hit_idx]);
        
        if (cl1_it != cluster_map.end() && cl2_it != cluster_map.end()) {
            size_t cl1_pos = cl1_it->second;
            size_t cl2_pos = cl2_it->second;

            const auto& bc1 = CLparticleLink_barcode[cl1_pos];
            const auto& bc2 = CLparticleLink_barcode[cl2_pos];
            for (int b1 : bc1) {
                for (int b2 : bc2) {
                    matched_cl1_pos.push_back(cl1_pos);
                    matched_cl2_pos.push_back(cl2_pos);
                    valid_hit_idx.push_back(hit_idx);
                    matched_barcodes_1.push_back(b1);
                    matched_barcodes_2.push_back(b2);
                }
            }
        }
    }
    return std::make_tuple(matched_cl1_pos, matched_cl2_pos, valid_hit_idx, 
                          matched_barcodes_1, matched_barcodes_2);
}

ROOT::VecOps::RVec<float> CalculateR(const ROOT::VecOps::RVec<double>& x, 
                                    const ROOT::VecOps::RVec<double>& y) {
    return sqrt(x*x + y*y);
}

ROOT::VecOps::RVec<float> CalculatePhi(const ROOT::VecOps::RVec<double>& x, 
                                      const ROOT::VecOps::RVec<double>& y) {
    ROOT::VecOps::RVec<float> phi;
    phi.reserve(x.size());
    for (size_t i = 0; i < x.size(); ++i) {
        phi.push_back(atan2(y[i], x[i]));
    }
    return phi;
}

ROOT::VecOps::RVec<float> CalculateEta(const ROOT::VecOps::RVec<float>& r, 
                                      const ROOT::VecOps::RVec<double>& z) {
    ROOT::VecOps::RVec<float> eta;
    eta.reserve(r.size());
    for (size_t i = 0; i < r.size(); ++i) {
        float theta = atan2(r[i], z[i]);
        eta.push_back(-log(tan(theta/2.0)));
    }
    return eta;
}

ROOT::VecOps::RVec<float> CalculateClusterSeparation(
    const ROOT::VecOps::RVec<double>& cl1_x, const ROOT::VecOps::RVec<double>& cl1_y, 
    const ROOT::VecOps::RVec<double>& cl1_z, const ROOT::VecOps::RVec<double>& cl2_x, 
    const ROOT::VecOps::RVec<double>& cl2_y, const ROOT::VecOps::RVec<double>& cl2_z) {
    
    ROOT::VecOps::RVec<float> separation;
    separation.reserve(cl1_x.size());
    
    for (size_t i = 0; i < cl1_x.size(); ++i) {
        float dx = cl1_x[i] - cl2_x[i];
        float dy = cl1_y[i] - cl2_y[i];
        float dz = cl1_z[i] - cl2_z[i];
        separation.push_back(sqrt(dx*dx + dy*dy + dz*dz));
    }
    return separation;
}

std::tuple<
    ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>,
    ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>,
    ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<float>,
    ROOT::VecOps::RVec<float>, ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>,
    ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>
>
MergeParticlesToHits(
    const ROOT::VecOps::RVec<int>& hit_particle_ids,
    const ROOT::VecOps::RVec<int>& particle_ids,
    const ROOT::VecOps::RVec<float>& particle_pt,
    const ROOT::VecOps::RVec<float>& particle_eta,
    const ROOT::VecOps::RVec<float>& particle_px,
    const ROOT::VecOps::RVec<float>& particle_py,
    const ROOT::VecOps::RVec<float>& particle_vx,
    const ROOT::VecOps::RVec<float>& particle_vy,
    const ROOT::VecOps::RVec<float>& particle_vz,
    const ROOT::VecOps::RVec<float>& particle_radius,
    const ROOT::VecOps::RVec<float>& particle_status,
    const ROOT::VecOps::RVec<int>& particle_charge,
    const ROOT::VecOps::RVec<int>& particle_pdgid,
    const ROOT::VecOps::RVec<int>& particle_pass,
    const ROOT::VecOps::RVec<int>& particle_vProdNin,
    const ROOT::VecOps::RVec<int>& particle_vProdNout,
    const ROOT::VecOps::RVec<int>& particle_vProdStatus,
    const ROOT::VecOps::RVec<int>& particle_vProdBarcode,
    const ROOT::VecOps::RVec<int>& particle_barcode
) {
    std::unordered_map<int, size_t> particle_map;
    for (size_t i = 0; i < particle_ids.size(); ++i) {
        particle_map[particle_ids[i]] = i;
    }

    ROOT::VecOps::RVec<int> merged_particle_ids;
    ROOT::VecOps::RVec<float> merged_pt, merged_eta, merged_px, merged_py;
    ROOT::VecOps::RVec<float> merged_vx, merged_vy, merged_vz;
    ROOT::VecOps::RVec<float> merged_radius, merged_status;
    ROOT::VecOps::RVec<int> merged_charge, merged_pdgid, merged_pass;
    ROOT::VecOps::RVec<int> merged_vProdNin, merged_vProdNout, merged_vProdStatus, merged_vProdBarcode;
    ROOT::VecOps::RVec<int> merged_primary;

    size_t n_hits = hit_particle_ids.size();
    merged_particle_ids.reserve(n_hits);
    merged_pt.reserve(n_hits); merged_eta.reserve(n_hits); 
    merged_px.reserve(n_hits); merged_py.reserve(n_hits);
    merged_vx.reserve(n_hits); merged_vy.reserve(n_hits); merged_vz.reserve(n_hits);
    merged_radius.reserve(n_hits); merged_status.reserve(n_hits);
    merged_charge.reserve(n_hits); merged_pdgid.reserve(n_hits); merged_pass.reserve(n_hits);
    merged_vProdNin.reserve(n_hits); merged_vProdNout.reserve(n_hits); 
    merged_vProdStatus.reserve(n_hits); merged_vProdBarcode.reserve(n_hits);
    merged_primary.reserve(n_hits);

    for (size_t hit_idx = 0; hit_idx < n_hits; ++hit_idx) {
        int hit_particle_id = hit_particle_ids[hit_idx];

        if (hit_particle_id == 0 || particle_map.find(hit_particle_id) == particle_map.end()) {
            merged_particle_ids.push_back(0);
            merged_pt.push_back(0.0f); merged_eta.push_back(0.0f); 
            merged_px.push_back(0.0f); merged_py.push_back(0.0f);
            merged_vx.push_back(0.0f); merged_vy.push_back(0.0f); merged_vz.push_back(0.0f);
            merged_radius.push_back(0.0f); merged_status.push_back(0.0f);
            merged_charge.push_back(0); merged_pdgid.push_back(0); merged_pass.push_back(0);
            merged_vProdNin.push_back(0); merged_vProdNout.push_back(0); 
            merged_vProdStatus.push_back(0); merged_vProdBarcode.push_back(0);
            merged_primary.push_back(0);
        } else {
            size_t idx = particle_map[hit_particle_id];
            merged_particle_ids.push_back(hit_particle_id);
            merged_pt.push_back(particle_pt[idx]);
            merged_eta.push_back(particle_eta[idx]);
            merged_px.push_back(particle_px[idx]);
            merged_py.push_back(particle_py[idx]);
            merged_vx.push_back(particle_vx[idx]);
            merged_vy.push_back(particle_vy[idx]);
            merged_vz.push_back(particle_vz[idx]);
            merged_radius.push_back(particle_radius[idx]);
            merged_status.push_back(particle_status[idx]);
            merged_charge.push_back(particle_charge[idx]);
            merged_pdgid.push_back(particle_pdgid[idx]);
            merged_pass.push_back(particle_pass[idx]);
            merged_vProdNin.push_back(particle_vProdNin[idx]);
            merged_vProdNout.push_back(particle_vProdNout[idx]);
            merged_vProdStatus.push_back(particle_vProdStatus[idx]);
            merged_vProdBarcode.push_back(particle_vProdBarcode[idx]);
            merged_primary.push_back(particle_barcode[idx] < 200000 ? 1 : 0);
        }
    }

    return std::make_tuple(
        merged_particle_ids, merged_pt, merged_eta, merged_px, merged_py,
        merged_vx, merged_vy, merged_vz, merged_radius, merged_status,
        merged_charge, merged_pdgid, merged_pass,
        merged_vProdNin, merged_vProdNout, merged_vProdStatus, merged_vProdBarcode,
        merged_primary
    );
}

ROOT::VecOps::RVec<float> CalculateDistanceFromProduction(
    const ROOT::VecOps::RVec<double>& hit_x, const ROOT::VecOps::RVec<double>& hit_y, 
    const ROOT::VecOps::RVec<double>& hit_z, const ROOT::VecOps::RVec<float>& particle_vx, 
    const ROOT::VecOps::RVec<float>& particle_vy, const ROOT::VecOps::RVec<float>& particle_vz) {
    
    ROOT::VecOps::RVec<float> distances;
    distances.reserve(hit_x.size());
    
    for (size_t i = 0; i < hit_x.size(); ++i) {
        float dx = hit_x[i] - particle_vx[i];
        float dy = hit_y[i] - particle_vy[i];
        float dz = hit_z[i] - particle_vz[i];
        distances.push_back(sqrt(dx*dx + dy*dy + dz*dz));
    }
    return distances;
}

std::tuple<ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>>
BuildTrackEdges(const ROOT::VecOps::RVec<int>& hit_ids,
                const ROOT::VecOps::RVec<int>& particle_ids,
                const ROOT::VecOps::RVec<float>& hit_distances,
                const ROOT::VecOps::RVec<int>& hit_barrel_endcap,
                const ROOT::VecOps::RVec<int>& hit_layer_disk,
                const ROOT::VecOps::RVec<int>& hit_eta_module,
                const ROOT::VecOps::RVec<int>& hit_phi_module,
                const ROOT::VecOps::RVec<string>& hit_hardware) {
    
    std::map<int, std::vector<std::pair<float, int>>> particle_hits;
    
    for (size_t i = 0; i < hit_ids.size(); ++i) {
        if (particle_ids[i] != 0) {
            particle_hits[particle_ids[i]].push_back({hit_distances[i], static_cast<int>(i)});
        }
    }
    
    for (auto& [particle_id, hits] : particle_hits) {
        std::sort(hits.begin(), hits.end());
    }
    
    ROOT::VecOps::RVec<int> edge_from, edge_to;
    
    for (auto& [particle_id, hits] : particle_hits) {
        if (hits.size() < 2) continue;
        
        std::map<std::tuple<int,int,int,int,string>, std::vector<int>> module_hits;
        
        for (auto& [dist, hit_idx] : hits) {
            auto module_key = std::make_tuple(hit_barrel_endcap[hit_idx], hit_layer_disk[hit_idx], 
                                            hit_eta_module[hit_idx], hit_phi_module[hit_idx], 
                                            hit_hardware[hit_idx]);
            module_hits[module_key].push_back(hit_idx);
        }
        
        std::vector<std::vector<int>> module_groups;
        for (auto& [module, hit_indices] : module_hits) {
            module_groups.push_back(hit_indices);
        }
        
        for (size_t i = 0; i < module_groups.size() - 1; ++i) {
            for (int hit_from : module_groups[i]) {
                for (int hit_to : module_groups[i+1]) {
                    edge_from.push_back(hit_from);
                    edge_to.push_back(hit_to);
                }
            }
        }
    }
    
    return std::make_tuple(edge_from, edge_to);
}

ROOT::VecOps::RVec<bool> ApplyQualityCuts(const ROOT::VecOps::RVec<float>& pt,
                                          const ROOT::VecOps::RVec<float>& eta,
                                          const ROOT::VecOps::RVec<int>& charge,
                                          float min_pt = 0.5, float max_eta = 2.5) {
    ROOT::VecOps::RVec<bool> cuts;
    cuts.reserve(pt.size());
    for (size_t i = 0; i < pt.size(); ++i) {
        cuts.push_back(pt[i] > min_pt && abs(eta[i]) < max_eta && charge[i] != 0);
    }
    return cuts;
}

ROOT::VecOps::RVec<int> CountHitsPerParticle(const ROOT::VecOps::RVec<int>& particle_ids) {
    std::unordered_map<int, int> hit_counts;
    
    for (int particle_id : particle_ids) {
        if (particle_id != 0) {
            hit_counts[particle_id]++;
        }
    }
    
    ROOT::VecOps::RVec<int> nhits_per_hit;
    nhits_per_hit.reserve(particle_ids.size());
    
    for (int particle_id : particle_ids) {
        if (particle_id == 0) {
            nhits_per_hit.push_back(-1);
        } else {
            nhits_per_hit.push_back(hit_counts[particle_id]);
        }
    }
    
    return nhits_per_hit;
}
'''


class ComprehensiveGNN4ITkProcessor:
    """Processor with enhanced RNTuple support via RDF"""
    
    def __init__(self, filename: str, num_threads: int = 1):
        self.filename = filename
        self.num_threads = num_threads
        self.perf = PerformanceMonitor()
        self.rdf = None
        self._declare_cpp_functions()
        
    def _declare_cpp_functions(self):
        """Declare C++ functions for RDataFrame"""
        ROOT.gInterpreter.Declare(cpp_code)
        self.perf.print_stats("C++ Function Declaration")

    def load_data(self) -> ROOT.RDataFrame:
        """
        Load RNTuple data and create RDataFrame with explicit RNTupleDS support.
        Returns an RDataFrame instance backed by RNTupleDS.
        """
        if '*' in self.filename or '?' in self.filename:
            file_list = glob.glob(self.filename)
            if not file_list:
                raise FileNotFoundError(f"No files matching pattern: {self.filename}")
            
            file_list.sort()
            print(f"Found {len(file_list)} files matching pattern")
            
            # Create vector of file paths for RNTupleDS
            root_file_list = ROOT.std.vector('string')()
            for f in file_list:
                root_file_list.push_back(f)
            
            # Create RNTupleDS explicitly
            # ntuple_ds = ROOT.RDF.RNTupleDS("GNN4ITk", root_file_list)
            # Create RDataFrame from the data source
            self.rdf = ROOT.RDF.FromRNTuple("GNN4ITk", root_file_list)
        else:
            if not os.path.exists(self.filename):
                raise FileNotFoundError(f"File not found: {self.filename}")
            print(f"Loading single file: {self.filename}")
            
            # Create RNTupleDS for single file
            # ntuple_ds = ROOT.RDF.RNTupleDS("GNN4ITk", self.filename)
            self.rdf = ROOT.RDF.FromRNTuple("GNN4ITk", self.filename)
        
        # Get entry count using Count()
        n_entries = self.rdf.Count().GetValue()
        print(f"Total entries in RNTuple: {n_entries}")
        
        self.perf.print_stats("RNTuple Loading via RDF")
        return self.rdf
    
    def load_processed_data(self, input_path: str, format: str = "root") -> ROOT.RDataFrame:
        """Load previously processed data"""
        print(f"Loading previously processed data from {input_path} ({format})...")
        
        if format.lower() == "root":
            if not input_path.endswith('.root'):
                input_path = input_path + '.rntuple.root'
            if not os.path.exists(input_path):
                raise FileNotFoundError(f"Processed file not found: {input_path}")
            rdf = ROOT.RDF.FromRNTuple("GNN4ITk_processed", input_path)
            n_entries = rdf.Count().GetValue()
            print(f"✓ Loaded {n_entries} entries from {input_path}")
            self.perf.print_stats("Loaded Processed Data")
            return rdf
        else:
            raise NotImplementedError(f"Loading from {format} format not yet implemented. Please use ROOT format.")

    def merge_particles_to_hits(self, rdf: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Merge particle information to hits using RDataFrame operations"""
        merged_rdf = (
            rdf
            .Define("particle_merge_result",
                "MergeParticlesToHits("
                "particle_id, Part_barcode, "
                "Part_pt, Part_eta, Part_px, Part_py, "
                "Part_vx, Part_vy, Part_vz, "
                "Part_radius, Part_status, Part_charge, "
                "Part_pdg_id, Part_passed, "
                "Part_vProdNin, Part_vProdNout, Part_vProdStatus, Part_vProdBarcode, "
                "Part_barcode)")
            .Define("hit_particle_id",       "std::get<0>(particle_merge_result)")
            .Define("hit_particle_pt",       "std::get<1>(particle_merge_result)")
            .Define("hit_particle_eta",      "std::get<2>(particle_merge_result)")
            .Define("hit_particle_phi",      "CalculatePhi(std::get<3>(particle_merge_result), std::get<4>(particle_merge_result))")
            .Define("hit_particle_vx",       "std::get<5>(particle_merge_result)")
            .Define("hit_particle_vy",       "std::get<6>(particle_merge_result)")
            .Define("hit_particle_vz",       "std::get<7>(particle_merge_result)")
            .Define("hit_particle_radius",   "std::get<8>(particle_merge_result)")
            .Define("hit_particle_status",   "std::get<9>(particle_merge_result)")
            .Define("hit_particle_charge",   "std::get<10>(particle_merge_result)")
            .Define("hit_particle_pdgid",    "std::get<11>(particle_merge_result)")
            .Define("hit_particle_pass",     "std::get<12>(particle_merge_result)")
            .Define("hit_particle_vProdNin", "std::get<13>(particle_merge_result)")
            .Define("hit_particle_vProdNout","std::get<14>(particle_merge_result)")
            .Define("hit_particle_vProdStatus", "std::get<15>(particle_merge_result)")
            .Define("hit_particle_vProdBarcode","std::get<16>(particle_merge_result)")
            .Define("hit_particle_primary",  "std::get<17>(particle_merge_result)")
            .Define("hit_particle_nhits", "CountHitsPerParticle(hit_particle_id)")
        )
        self.perf.print_stats("Particle-Hit Merging")
        return merged_rdf

    def perform_cluster_matching(self, rdf: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Perform cluster matching using RDataFrame operations"""
        matched_rdf = (rdf
            .Define("cluster_matching_result", 
                "ClusterMatching(SPCL1_index, SPCL2_index, CLindex, CLparticleLink_barcode)")
            .Define("matched_cl1_pos", "std::get<0>(cluster_matching_result)")
            .Define("matched_cl2_pos", "std::get<1>(cluster_matching_result)")
            .Define("valid_hit_indices", "std::get<2>(cluster_matching_result)")
            .Define("particle_id_1", "std::get<3>(cluster_matching_result)")
            .Define("particle_id_2", "std::get<4>(cluster_matching_result)")
            .Define("particle_id", "ROOT::VecOps::Where(particle_id_1 == particle_id_2, particle_id_1, -1)")
            .Define("hit_id", "Take(SPindex, valid_hit_indices)")
            .Define("hit_x", "Take(SPx, valid_hit_indices)")
            .Define("hit_y", "Take(SPy, valid_hit_indices)")
            .Define("hit_z", "Take(SPz, valid_hit_indices)")
            .Define("hit_r", "CalculateR(hit_x, hit_y)")
            .Define("hit_phi", "CalculatePhi(hit_x, hit_y)")
            .Define("hit_eta", "CalculateEta(hit_r, hit_z)")
            .Define("hit_radius", "Take(SPradius, valid_hit_indices)")
            .Define("hit_overlap", "Take(SPisOverlap, valid_hit_indices)")
            .Define("hit_barrel_endcap", "Take(CLbarrel_endcap, matched_cl1_pos)")
            .Define("hit_layer_disk", "Take(CLlayer_disk, matched_cl1_pos)")
            .Define("hit_eta_module", "Take(CLeta_module, matched_cl1_pos)")
            .Define("hit_phi_module", "Take(CLphi_module, matched_cl1_pos)")
            .Define("hit_hardware", "Take(CLhardware, matched_cl1_pos)")
            .Define("hit_module_id", "Take(CLmoduleID, matched_cl1_pos)")
            # Cluster 1 properties
            .Define("hit_cluster_id_1", "Take(CLindex, matched_cl1_pos)")
            .Define("hit_cluster_x_1", "Take(CLx, matched_cl1_pos)")
            .Define("hit_cluster_y_1", "Take(CLy, matched_cl1_pos)")
            .Define("hit_cluster_z_1", "Take(CLz, matched_cl1_pos)")
            .Define("hit_cluster_norm_x_1", "Take(CLnorm_x, matched_cl1_pos)")
            .Define("hit_cluster_norm_y_1", "Take(CLnorm_y, matched_cl1_pos)")
            .Define("hit_cluster_norm_z_1", "Take(CLnorm_z, matched_cl1_pos)")
            .Define("hit_cluster_size_1", "Take(CLside, matched_cl1_pos)")
            .Define("hit_cluster_charge_count_1", "Take(CLcharge_count, matched_cl1_pos)")
            .Define("hit_cluster_count_1", "Take(CLpixel_count, matched_cl1_pos)")
            # Cluster 2 properties
            .Define("hit_cluster_id_2", "Take(CLindex, matched_cl2_pos)")
            .Define("hit_cluster_x_2", "Take(CLx, matched_cl2_pos)")
            .Define("hit_cluster_y_2", "Take(CLy, matched_cl2_pos)")
            .Define("hit_cluster_z_2", "Take(CLz, matched_cl2_pos)")
            .Define("hit_cluster_charge_count_2", "Take(CLcharge_count, matched_cl2_pos)")
            .Define("hit_cluster_count_2", "Take(CLpixel_count, matched_cl2_pos)")
        )
        self.perf.print_stats("Cluster Matching")
        return matched_rdf
    
    def add_derived_features(self, rdf: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Add all derived physics quantities using RDataFrame operations"""
        physics_rdf = (rdf
            .Define("hit_cluster_separation", 
                "CalculateClusterSeparation(hit_cluster_x_1, hit_cluster_y_1, hit_cluster_z_1, "
                "                          hit_cluster_x_2, hit_cluster_y_2, hit_cluster_z_2)")
            .Define("hit_cluster_charge_total", "hit_cluster_charge_count_1 + hit_cluster_charge_count_2")
            .Define("hit_cluster_pixel_total", "hit_cluster_count_1 + hit_cluster_count_2")
            .Define("hit_cluster_charge_asymmetry", 
                "abs(hit_cluster_charge_count_1 - hit_cluster_charge_count_2) / "
                "(hit_cluster_charge_count_1 + hit_cluster_charge_count_2 + 1e-6f)")
            .Define("hit_same_layer", "Take(CLlayer_disk, matched_cl1_pos) == Take(CLlayer_disk, matched_cl2_pos)")
            .Define("hit_barrel_region", "hit_barrel_endcap == 0")
            .Define("hit_is_noise", "hit_particle_id == 0")
            .Define("hit_is_signal", "hit_particle_id > 0")
            .Define("n_hits_total", "static_cast<int>(hit_id.size())")
            .Define("n_hits_signal", "static_cast<int>(Sum(hit_is_signal))")
            .Define("n_hits_noise", "static_cast<int>(Sum(hit_is_noise))")
            .Define("signal_fraction", "n_hits_total > 0 ? static_cast<float>(n_hits_signal) / n_hits_total : 0.0f")
            .Define("barrel_fraction", "n_hits_total > 0 ? static_cast<float>(Sum(hit_barrel_region)) / n_hits_total : 0.0f")
        )
        self.perf.print_stats("Derived Feature Calculation")
        return physics_rdf

    def build_track_edges(self, rdf: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Build track edges between hits on the same particle"""
        track_rdf = (rdf
            .Define("hit_distance_from_production",
                "CalculateDistanceFromProduction(hit_x, hit_y, hit_z, "
                "hit_particle_vx, hit_particle_vy, hit_particle_vz)")
            .Define("track_edge_result",
                "BuildTrackEdges(hit_id, hit_particle_id, hit_distance_from_production, "
                "hit_barrel_endcap, hit_layer_disk, hit_eta_module, hit_phi_module, hit_hardware)")
            .Define("track_edge_from", "std::get<0>(track_edge_result)")
            .Define("track_edge_to", "std::get<1>(track_edge_result)")
            .Define("n_track_edges", "static_cast<int>(track_edge_from.size())")
        )
        self.perf.print_stats("Track Edge Building")
        return track_rdf

    def apply_quality_cuts(self, rdf: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Apply quality cuts and filters using RDataFrame operations"""
        filtered_rdf = (rdf
            .Filter("n_hits_total > 0", "Events with hits")
            .Filter("n_hits_signal > 0", "Events with signal hits")
            .Define("quality_cuts", "ApplyQualityCuts(hit_particle_pt, hit_particle_eta, hit_particle_charge)")
        )
        self.perf.print_stats("Quality Cuts and Filtering")
        return filtered_rdf

    def create_summary_statistics(self, rdf: ROOT.RDataFrame, include_edges: bool = False) -> Dict:
        """Generate summary statistics using RDataFrame aggregations"""
        print("Generating summary statistics...")
        
        n_events = rdf.Count().GetValue()
        if n_events == 0:
            return {"n_events": 0}
        
        summary = {
            "n_events": n_events,
            "total_hits": int(rdf.Sum("n_hits_total").GetValue()),
            "total_signal_hits": int(rdf.Sum("n_hits_signal").GetValue()),
            "mean_hits": rdf.Mean("n_hits_total").GetValue(),
        }
        
        # Only include edge statistics if they were computed
        if include_edges:
            summary["total_track_edges"] = int(rdf.Sum("n_track_edges").GetValue())
            summary["mean_edges"] = rdf.Mean("n_track_edges").GetValue()
        
        print(f"=== Summary | Events: {n_events} | Total Hits: {summary['total_hits']} ===")
        return summary

    def export_data(self, rdf: ROOT.RDataFrame, output_path: str, format: str = "csv", include_edges: bool = False):
        """Export processed data from RDataFrame"""
        print(f"Exporting data to {output_path} ({format})...")
        
        # Base columns always included
        cols = [
            "run_number", "event_number", "hit_id", "hit_x", "hit_y", "hit_z", 
            "hit_r", "hit_phi", "hit_eta", "hit_particle_id", "hit_particle_pt", 
            "hit_cluster_charge_total", "n_hits_signal", "n_hits_total"]
        
        # Add edge columns only if they exist in the dataframe
        if include_edges:
            cols.extend(["track_edge_from", "track_edge_to", "n_track_edges"])
        else:
            cols.extend(["hit_particle_vx", "hit_particle_vy", "hit_particle_vz", "hit_barrel_endcap", "hit_layer_disk", "hit_eta_module", "hit_phi_module", "hit_hardware"])
        
        if format.lower() == "root":
            # comprAlgo = getattr(ROOT.RCompressionSetting.EAlgorithm, "kZLIB")
            # opts = ROOT.RDF.RSnapshotOptions("RECREATE", comprAlgo, 0, 0, 99, False)
            # opts.fOutputFormat = ROOT.RDF.ESnapshotOutputFormat.kRNTuple
            opts = ROOT.RDF.RSnapshotOptions()
            opts.fMode = "RECREATE"
            opts.fOutputFormat = ROOT.RDF.ESnapshotOutputFormat.kRNTuple 
            rsnapslot = rdf.Snapshot("GNN4ITk_processed", output_path+'.rntuple.root', cols,opts)
            print(f"Data exported to ROOT RNTuple: {output_path}.rntuple.root ({rsnapslot.Count().GetValue()}) with size {os.path.getsize(output_path+'.rntuple.root')/(1024*1024):.2f} MB")
        else:
            numpy_data = rdf.AsNumpy(cols)
            
            if format.lower() == "csv":
                pd.DataFrame(numpy_data).to_csv(output_path, index=False)
            elif format.lower() == "numpy":
                np.savez_compressed(output_path, **numpy_data)
        self.perf.print_stats(f"Data Export ({format})")

    def build_and_export_edges(self, rdf: ROOT.RDataFrame, output_file: str, format: str = "csv"):
        """
        Separate method to build track edges and export them.
        This can be called after the main pipeline processing.
        
        Args:
            rdf: Input RDataFrame (should already have hits processed)
            output_file: Path for output file with edges
            format: Output format (csv, root, numpy)
        """
        print("\n=== Building Track Edges (Optional Post-Processing) ===")
        
        # Use a fresh multithreading context for edge building
        # This helps avoid threading conflicts with previous operations
        try:
            # First, ensure we're in a clean threading state
            if ROOT.IsImplicitMTEnabled():
                print("Temporarily disabling ImplicitMT before edge building...")
                ROOT.DisableImplicitMT()
            
            # Re-enable with specified thread count
            with root_multithreading(self.num_threads):
                # Build the track edges
                print("Building track edges...")
                track_rdf = self.build_track_edges(rdf)
                
                # Force computation by triggering the graph
                n_entries = track_rdf.Count().GetValue()
                print(f"Edge building completed for {n_entries} entries")
                
                # Export with edges included
                print("Exporting data with edges...")
                self.export_data(track_rdf, output_file, format, include_edges=True)
                
                # Return summary with edge statistics
                summary = self.create_summary_statistics(track_rdf, include_edges=True)
                
                print(f"Track edges built and exported to: {output_file}")
                print(f"Total track edges: {summary.get('total_track_edges', 0)}")
                
            return track_rdf, summary
            
        except Exception as e:
            print(f"Error during edge building: {e}")
            import traceback
            traceback.print_exc()
            raise

    def process_complete_pipeline(self, output_file: str = None, 
                                 format: str = "csv") -> Tuple[ROOT.RDataFrame, Dict]:
        """Complete pipeline execution with RNTuple/RDataFrame (without track edges)
        
        Args:
            output_file: Path for output file
            format: Output format (csv, root, numpy)
        
        Returns:
            Tuple of (final RDataFrame, summary statistics dictionary)
        """
        with root_multithreading(self.num_threads):
            rdf = self.load_data()
            ROOT.RDF.Experimental.AddProgressBar(rdf);
            matched_rdf = self.perform_cluster_matching(rdf)
            merged_rdf = self.merge_particles_to_hits(matched_rdf)
            physics_rdf = self.add_derived_features(merged_rdf)
            physics_rdf.Report()
            
            # Apply quality cuts on the physics dataframe
            final_rdf = self.apply_quality_cuts(physics_rdf)
            
            # Generate summary statistics (without edges)
            summary = self.create_summary_statistics(final_rdf, include_edges=False)
            print("final_rdf count:", final_rdf.Count().GetValue())
        
        # Export data (without track edges)
        if output_file:
            self.export_data(final_rdf, output_file, format, include_edges=False)
        
        return final_rdf, summary

    def create_sample_event_display(self, rdf: ROOT.RDataFrame, event_idx: int = 0):
        """Create a sample display of processed event data"""
        print(f"\n=== Sample Event Display (Event {event_idx}) ===")
        
        sample_columns = [
            "run_number", "event_number", "n_hits_total", "n_hits_signal",
            "hit_id", "hit_x", "hit_y", "hit_z", "hit_r", "hit_phi", "hit_eta",
            "hit_particle_id", "hit_particle_pt", "hit_cluster_charge_total"
        ]
        
        # Check if track edges exist in this dataframe
        try:
            # Try to include edge columns if they exist
            rdf_test = rdf.Define("_test_edges", "n_track_edges")
            sample_columns.extend(["track_edge_from", "track_edge_to", "n_track_edges"])
            has_edges = True
        except:
            has_edges = False
        
        try:
            numpy_data = rdf.AsNumpy(sample_columns)
            
            if len(numpy_data['run_number']) > event_idx:
                print(f"Run: {numpy_data['run_number'][event_idx]}")
                print(f"Event: {numpy_data['event_number'][event_idx]}")
                print(f"Total Hits: {numpy_data['n_hits_total'][event_idx]}")
                print(f"Signal Hits: {numpy_data['n_hits_signal'][event_idx]}")
                if has_edges:
                    print(f"Track Edges: {numpy_data['n_track_edges'][event_idx]}")
                
                hit_ids = numpy_data['hit_id'][event_idx]
                if len(hit_ids) > 0:
                    print(f"\nFirst 5 hits:")
                    for i in range(min(5, len(hit_ids))):
                        print(f"  Hit {hit_ids[i]}: (x={numpy_data['hit_x'][event_idx][i]:.2f}, "
                              f"y={numpy_data['hit_y'][event_idx][i]:.2f}, "
                              f"z={numpy_data['hit_z'][event_idx][i]:.2f})")
                
                if has_edges:
                    edge_from = numpy_data['track_edge_from'][event_idx]
                    edge_to = numpy_data['track_edge_to'][event_idx]
                    if len(edge_from) > 0:
                        print(f"\nFirst 5 track edges:")
                        for i in range(min(5, len(edge_from))):
                            print(f"  Edge {i}: {edge_from[i]} -> {edge_to[i]}")
            else:
                print("No events available for display")
                
        except Exception as e:
            print(f"Error creating event display: {e}")
        
        print("=" * 50)


def main():
    """Main processing function"""
    parser = argparse.ArgumentParser(
        description='Process GNN4ITk RNTuple data with RDataFrame')
    parser.add_argument('--input_file', '-i', required=True, 
                       help='Input RNTuple file path (supports wildcards)')
    parser.add_argument('--output', '-o', help='Output file path')
    parser.add_argument('--format', '-f', choices=['csv', 'root', 'numpy'], 
                       default='csv', help='Output format (default: csv)')
    parser.add_argument('--threads', '-t', type=int, default=4, 
                       help='Number of threads for processing (default: 4)')
    # parser.add_argument('--edge-threads', type=int, default=None,
    #                    help='Number of threads for edge building (default: min(threads, 4) for stability)')
    parser.add_argument('--build-edges', action='store_true',
                       help='Build track edges (separate post-processing step)')
    parser.add_argument('--edges-output', help='Output file for track edges (if different from main output)')
    parser.add_argument('--show-sample', action='store_true',
                       help='Display sample event data')
    
    args = parser.parse_args()
    
    try:
        print(f"[{time.time()-start_time:.2f}] Initializing RNTuple processor...")
        # verbosity = ROOT.RLogScopedVerbosity(ROOT.Detail.RDF.RDFLogChannel(), ROOT.ELogLevel.kLogDebug+10)
        processor = ComprehensiveGNN4ITkProcessor(args.input_file, num_threads=args.threads)
        
        # Check if main output already exists (for ROOT format only)
        output_exists = False
        if args.output and args.format.lower() == 'root':
            output_path = args.output if args.output.endswith('.root') else args.output + '.rntuple.root'
            if os.path.exists(output_path):
                output_exists = True
                print(f"[{time.time()-start_time:.2f}] ✓ Found existing output file: {output_path}")
        
        # Load existing data or process from scratch
        if args.build_edges:
            edges_output = args.edges_output if args.edges_output else args.output+'_with_edges'
            print(f"\n[{time.time()-start_time:.2f}] Building track edges (optional post-processing)...")

            processor.num_threads = args.threads
            with root_multithreading(args.threads):
                if output_exists:
                    print(f"[{time.time()-start_time:.2f}] Loading existing processed data for edge building...")
                    final_rdf = processor.load_processed_data(args.output, args.format)
                    stats = processor.create_summary_statistics(final_rdf, include_edges=False)
                else:
                    print(f"[{time.time()-start_time:.2f}] Running complete processing pipeline before edge building...")
                    final_rdf, stats = processor.process_complete_pipeline(
                        output_file=args.output, format=args.format)
                try:
                    final_rdf, edge_stats = processor.build_and_export_edges(final_rdf, edges_output, args.format)
                    stats.update(edge_stats)  # Merge edge statistics into main stats
                except Exception as e:
                    print(f"WARNING: Edge building failed: {e}")
                    print("Main pipeline results are still available in the output file.")
                    print("You may want to retry edge building with fewer threads (--edge-threads 1 or 2)")
            
        else:
            print(f"[{time.time()-start_time:.2f}] Running complete processing pipeline...")
            final_rdf, stats = processor.process_complete_pipeline(
                output_file=args.output, format=args.format)


        # Optionally build track edges as a separate step
        # if args.build_edges:
        #     edges_output = args.edges_output if args.edges_output else args.output.replace('.', '_with_edges.')
        #     print(f"\n[{time.time()-start_time:.2f}] Building track edges (optional post-processing)...")
            
        #     # For edge building, use configurable thread count (default: fewer threads for stability)
        #     edge_threads = args.edge_threads if args.edge_threads is not None else min(args.threads, 4)
        #     if edge_threads != args.threads:
        #         print(f"Note: Using {edge_threads} threads for edge building (configured for stability)")
            
        #     # Update processor thread count for edge building
        #     processor.num_threads = edge_threads
            
        #     try:
        #         edge_rdf, edge_stats = processor.build_and_export_edges(final_rdf, edges_output, args.format)
        #         stats.update(edge_stats)  # Merge edge statistics into main stats
        #     except Exception as e:
        #         print(f"WARNING: Edge building failed: {e}")
        #         print("Main pipeline results are still available in the output file.")
        #         print("You may want to retry edge building with fewer threads (--edge-threads 1 or 2)")
        #         if not args.output:
        #             return 1  # Exit with error if no main output was saved
        
        if args.show_sample:
            print(f"[{time.time()-start_time:.2f}] Generating sample event display...")
            with root_multithreading(args.threads):
                processor.create_sample_event_display(final_rdf)
        
        print(f"\nProcessing completed successfully!")
        print(f"Processed {stats['n_events']} events")
        print(f"Total hits: {stats['total_hits']}")
        if 'total_track_edges' in stats:
            print(f"Total track edges: {stats['total_track_edges']}")
        
        if args.output:
            print(f"[{time.time()-start_time:.2f}] Results exported to: {args.output}")
        if args.build_edges and args.edges_output:
            print(f"[{time.time()-start_time:.2f}] Track edges exported to: {args.edges_output}")
            
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    print(f"[{time.time()-start_time:.2f}] All done!")
    return 0


if __name__ == "__main__":
    print(f"[{time.time()-start_time:.2f}] Starting GNN4ITk RNTuple processing...")
    exit(main())
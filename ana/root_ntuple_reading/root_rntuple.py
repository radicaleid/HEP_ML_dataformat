#!/usr/bin/env python3
"""
ROOT NTuple Reader for GNN4ITk data using RDataFrame
Extended version that prepares ALL variables needed for graph building
Combines particle, hit, and cluster information with maximum performance
"""

import ROOT
import numpy as np
import pandas as pd
import time
import psutil
import os,sys,glob
from typing import Tuple, Dict, List
import argparse

# # Add parent directory to path for importing perfmon_utils
# script_dir = os.path.dirname(os.path.abspath(__file__))
# parent_dir = os.path.dirname(script_dir)
# sys.path.insert(0, parent_dir)

# from perfmon_utils import PerformanceMonitor

# Enable ROOT's implicit multi-threading for better performance
# ROOT.EnableImplicitMT(4)  # Adjust based on your CPU cores
class PerformanceMonitor:
    """Simple performance monitoring class"""
    
    def __init__(self):
        self.start_time = time.time()
        self.process = psutil.Process(os.getpid())
        self.initial_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        self.last_time = self.start_time
    
    def print_stats(self, operation: str):
        """Print timing and memory statistics"""
        current_time = time.time()
        current_memory = self.process.memory_info().rss / 1024 / 1024  # MB
        
        elapsed = current_time - self.last_time
        memory_delta = current_memory - self.initial_memory
        
        print(f"=== Performance Report: {operation} ===")
        print(f"Time: {elapsed:.2f} seconds")
        print(f"Memory Delta: {memory_delta:.1f} MB")
        print(f"Total Memory: {current_memory:.1f} MB")
        print("=" * 50)
        
        self.last_time = current_time
        self.initial_memory = current_memory
        
# Define comprehensive cluster matching and physics calculation functions in ROOT C++
ROOT.gInterpreter.Declare('''
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include <tuple>
#include <algorithm>
#include <cmath>

// Cluster matching with hash map for O(1) lookups + barcode padding
std::tuple<
    ROOT::VecOps::RVec<int>,                   // matched_cl1_pos
    ROOT::VecOps::RVec<int>,                   // matched_cl2_pos
    ROOT::VecOps::RVec<int>,                   // valid_hit_idx
    ROOT::VecOps::RVec<int>,                   // matched_barcodes_1 (flattened)
    ROOT::VecOps::RVec<int>                    // matched_barcodes_2 (flattened)
>
ClusterMatching(const ROOT::VecOps::RVec<int>& hit_cl1_idx, 
                            const ROOT::VecOps::RVec<int>& hit_cl2_idx,
                            const ROOT::VecOps::RVec<int>& cluster_indices,
                            const ROOT::VecOps::RVec<ROOT::VecOps::RVec<int>>& CLparticleLink_barcode) {
    
    // Build hash map for O(1) cluster index lookups
    std::unordered_map<int, size_t> cluster_map;
    cluster_map.reserve(cluster_indices.size());
    
    for (size_t i = 0; i < cluster_indices.size(); ++i) {
        cluster_map[cluster_indices[i]] = i;
    }
    
    // Pre-allocate result vectors
    ROOT::VecOps::RVec<int> matched_cl1_pos, matched_cl2_pos, valid_hit_idx;
    ROOT::VecOps::RVec<int> matched_barcodes_1, matched_barcodes_2;
    
    matched_cl1_pos.reserve(hit_cl1_idx.size());
    matched_cl2_pos.reserve(hit_cl2_idx.size());
    valid_hit_idx.reserve(hit_cl1_idx.size());
    matched_barcodes_1.reserve(hit_cl1_idx.size());
    matched_barcodes_2.reserve(hit_cl2_idx.size());
    
    // Fast matching using hash map
    for (size_t hit_idx = 0; hit_idx < hit_cl1_idx.size(); ++hit_idx) {
        auto cl1_it = cluster_map.find(hit_cl1_idx[hit_idx]);
        auto cl2_it = cluster_map.find(hit_cl2_idx[hit_idx]);
        
        if (cl1_it != cluster_map.end() && cl2_it != cluster_map.end()) {
            size_t cl1_pos = cl1_it->second;
            size_t cl2_pos = cl2_it->second;

            // Merge and pad barcodes from both clusters
            const auto& bc1 = CLparticleLink_barcode[cl1_pos];
            const auto& bc2 = CLparticleLink_barcode[cl2_pos];
            for (size_t i = 0; i < bc1.size(); ++i) {
                for (size_t j = 0; j < bc2.size(); ++j) {
                    matched_cl1_pos.push_back(cl1_pos);
                    matched_cl2_pos.push_back(cl2_pos);
                    valid_hit_idx.push_back(hit_idx);
                    matched_barcodes_1.push_back(bc1[i]);
                    matched_barcodes_2.push_back(bc2[j]);
                }
            }
        }
    }
    
    return std::make_tuple(matched_cl1_pos, matched_cl2_pos, valid_hit_idx, matched_barcodes_1, matched_barcodes_2);
}

// Physics calculation functions
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
    const ROOT::VecOps::RVec<double>& cl1_x, const ROOT::VecOps::RVec<double>& cl1_y, const ROOT::VecOps::RVec<double>& cl1_z,
    const ROOT::VecOps::RVec<double>& cl2_x, const ROOT::VecOps::RVec<double>& cl2_y, const ROOT::VecOps::RVec<double>& cl2_z) {
    
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

// Particle merging function - merges particle info to hits based on particle ID
std::tuple<
    ROOT::VecOps::RVec<int>,    // merged_particle_ids
    ROOT::VecOps::RVec<float>,  // pt
    ROOT::VecOps::RVec<float>,  // eta
    ROOT::VecOps::RVec<float>,  // px
    ROOT::VecOps::RVec<float>,  // py
    ROOT::VecOps::RVec<float>,  // vx
    ROOT::VecOps::RVec<float>,  // vy
    ROOT::VecOps::RVec<float>,  // vz
    ROOT::VecOps::RVec<float>,  // radius
    ROOT::VecOps::RVec<float>,  // status (float in your file, could be int)
    ROOT::VecOps::RVec<int>,    // charge
    ROOT::VecOps::RVec<int>,    // pdgId
    ROOT::VecOps::RVec<int>,    // pass
    ROOT::VecOps::RVec<int>,    // vProdNIn
    ROOT::VecOps::RVec<int>,    // vProdNOut
    ROOT::VecOps::RVec<int>,    // vProdStatus
    ROOT::VecOps::RVec<int>,    // vProdBarcode
    ROOT::VecOps::RVec<int>     // primary flag
>
MergeParticlesToHits(
    const ROOT::VecOps::RVec<int>        &hit_particle_ids,
    const ROOT::VecOps::RVec<int>        &particle_ids,
    const ROOT::VecOps::RVec<float>      &particle_pt,
    const ROOT::VecOps::RVec<float>      &particle_eta,
    const ROOT::VecOps::RVec<float>      &particle_px,
    const ROOT::VecOps::RVec<float>      &particle_py,
    const ROOT::VecOps::RVec<float>      &particle_vx,
    const ROOT::VecOps::RVec<float>      &particle_vy,
    const ROOT::VecOps::RVec<float>      &particle_vz,
    const ROOT::VecOps::RVec<float>      &particle_radius,
    const ROOT::VecOps::RVec<float>      &particle_status,
    const ROOT::VecOps::RVec<int>        &particle_charge,
    const ROOT::VecOps::RVec<int>        &particle_pdgid,
    const ROOT::VecOps::RVec<int>        &particle_pass,
    const ROOT::VecOps::RVec<int>        &particle_vProdNin,
    const ROOT::VecOps::RVec<int>        &particle_vProdNout,
    const ROOT::VecOps::RVec<int>        &particle_vProdStatus,
    const ROOT::VecOps::RVec<int>        &particle_vProdBarcode,
    const ROOT::VecOps::RVec<int>        &particle_barcode // still needed for "primary" flag
) {
    // Build particle lookup map
    std::unordered_map<int, size_t> particle_map;
    for (size_t i = 0; i < particle_ids.size(); ++i) {
        particle_map[particle_ids[i]] = i;
    }

    // Initialize output vectors
    ROOT::VecOps::RVec<int> merged_particle_ids;
    ROOT::VecOps::RVec<float> merged_pt, merged_eta, merged_px, merged_py;
    ROOT::VecOps::RVec<float> merged_vx, merged_vy, merged_vz;
    ROOT::VecOps::RVec<float> merged_radius, merged_status;
    ROOT::VecOps::RVec<int> merged_charge, merged_pdgid, merged_pass;
    ROOT::VecOps::RVec<int> merged_vProdNin, merged_vProdNout, merged_vProdStatus, merged_vProdBarcode;
    ROOT::VecOps::RVec<int> merged_primary;

    // Reserve
    size_t n_hits = hit_particle_ids.size();
    merged_particle_ids.reserve(n_hits);
    merged_pt.reserve(n_hits); merged_eta.reserve(n_hits); merged_px.reserve(n_hits); merged_py.reserve(n_hits);
    merged_vx.reserve(n_hits); merged_vy.reserve(n_hits); merged_vz.reserve(n_hits);
    merged_radius.reserve(n_hits); merged_status.reserve(n_hits);
    merged_charge.reserve(n_hits); merged_pdgid.reserve(n_hits); merged_pass.reserve(n_hits);
    merged_vProdNin.reserve(n_hits); merged_vProdNout.reserve(n_hits); merged_vProdStatus.reserve(n_hits); merged_vProdBarcode.reserve(n_hits);
    merged_primary.reserve(n_hits);

    // Merge particle info to hits
    for (size_t hit_idx = 0; hit_idx < n_hits; ++hit_idx) {
        int hit_particle_id = hit_particle_ids[hit_idx];

        if (hit_particle_id == 0 || particle_map.find(hit_particle_id) == particle_map.end()) {
            // noise hit
            merged_particle_ids.push_back(0);
            merged_pt.push_back(0.0f); merged_eta.push_back(0.0f); merged_px.push_back(0.0f); merged_py.push_back(0.0f);
            merged_vx.push_back(0.0f); merged_vy.push_back(0.0f); merged_vz.push_back(0.0f);
            merged_radius.push_back(0.0f); merged_status.push_back(0.0f);
            merged_charge.push_back(0); merged_pdgid.push_back(0); merged_pass.push_back(0);
            merged_vProdNin.push_back(0); merged_vProdNout.push_back(0); merged_vProdStatus.push_back(0); merged_vProdBarcode.push_back(0);
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

// Calculate distance from production vertex
ROOT::VecOps::RVec<float> CalculateDistanceFromProduction(
    const ROOT::VecOps::RVec<double>& hit_x, const ROOT::VecOps::RVec<double>& hit_y, const ROOT::VecOps::RVec<double>& hit_z,
    const ROOT::VecOps::RVec<float>& particle_vx, const ROOT::VecOps::RVec<float>& particle_vy, const ROOT::VecOps::RVec<float>& particle_vz) {
    
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

// Build track edges from particle hits
std::tuple<ROOT::VecOps::RVec<int>, ROOT::VecOps::RVec<int>>
BuildTrackEdges(const ROOT::VecOps::RVec<int>& hit_ids,
                const ROOT::VecOps::RVec<int>& particle_ids,
                const ROOT::VecOps::RVec<float>& hit_distances,
                const ROOT::VecOps::RVec<int>& hit_barrel_endcap,
                const ROOT::VecOps::RVec<int>& hit_layer_disk,
                const ROOT::VecOps::RVec<int>& hit_eta_module,
                const ROOT::VecOps::RVec<int>& hit_phi_module,
                const ROOT::VecOps::RVec<string>& hit_hardware) {
    
    // Group hits by particle and module
    std::map<int, std::vector<std::pair<float, int>>> particle_hits; // particle_id -> (distance, hit_index)
    
    for (size_t i = 0; i < hit_ids.size(); ++i) {
        if (particle_ids[i] != 0) { // Only signal hits
            particle_hits[particle_ids[i]].push_back({hit_distances[i], static_cast<int>(i)});
        }
    }
    
    // Sort hits by distance from production vertex for each particle
    for (auto& [particle_id, hits] : particle_hits) {
        std::sort(hits.begin(), hits.end());
    }
    
    ROOT::VecOps::RVec<int> edge_from, edge_to;
    
    // Build edges between consecutive hits on same particle, grouped by modules
    for (auto& [particle_id, hits] : particle_hits) {
        if (hits.size() < 2) continue;
        
        // Group hits by module
        std::map<std::tuple<int,int,int,int,string>, std::vector<int>> module_hits; // (barrel_endcap, layer_disk, eta_module, phi_module, hardware) -> hit_indices
        
        for (auto& [dist, hit_idx] : hits) {
            auto module_key = std::make_tuple(hit_barrel_endcap[hit_idx], hit_layer_disk[hit_idx], 
                                            hit_eta_module[hit_idx], hit_phi_module[hit_idx], hit_hardware[hit_idx]);
            module_hits[module_key].push_back(hit_idx);
        }
        
        // Build edges between consecutive modules
        std::vector<std::vector<int>> module_groups;
        for (auto& [module, hit_indices] : module_hits) {
            module_groups.push_back(hit_indices);
        }
        
        // Connect consecutive modules
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

// Quality cuts function
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

// Count hits per particle
ROOT::VecOps::RVec<int> CountHitsPerParticle(const ROOT::VecOps::RVec<int>& particle_ids) {
    std::unordered_map<int, int> hit_counts;
    
    // Count hits for each particle
    for (int particle_id : particle_ids) {
        if (particle_id != 0) {
            hit_counts[particle_id]++;
        }
    }
    
    // Map back to hit array
    ROOT::VecOps::RVec<int> nhits_per_hit;
    nhits_per_hit.reserve(particle_ids.size());
    
    for (int particle_id : particle_ids) {
        if (particle_id == 0) {
            nhits_per_hit.push_back(-1); // Noise hits
        } else {
            nhits_per_hit.push_back(hit_counts[particle_id]);
        }
    }
    
    return nhits_per_hit;
}
''')

class ComprehensiveGNN4ITkProcessor:
    """Extended processor class for complete GNN4ITk graph data preparation"""
    
    def __init__(self, filename: str):
        self.filename = filename
        self.perf = PerformanceMonitor()
        self.df = None
        
    def load_data(self) -> ROOT.RDataFrame:
        """Load ROOT file(s) and create RDataFrame with wildcard support"""
        
        # Handle wildcard patterns
        if '*' in self.filename or '?' in self.filename:
            file_list = glob.glob(self.filename)
            if not file_list:
                raise FileNotFoundError(f"No files found matching pattern: {self.filename}")
            
            file_list.sort()
            print(f"Found {len(file_list)} files matching pattern: {self.filename}")
            
            # Convert to ROOT vector for better compatibility
            root_file_list = ROOT.std.vector('string')()
            for f in file_list:
                root_file_list.push_back(f)
                
            self.df = ROOT.RDataFrame("GNN4ITk", root_file_list)
        else:
            # Single file case
            if not os.path.exists(self.filename):
                raise FileNotFoundError(f"File not found: {self.filename}")
            print(f"Loading single file: {self.filename}")
            self.df = ROOT.RDataFrame("GNN4ITk", self.filename)
        
        print(f"Total events in dataset: {self.df.Count().GetValue()}")
        self.perf.print_stats("File Loading")
        return self.df

    def merge_particles_to_hits(self, df: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Merge particle information to hits using RDataFrame operations"""
        
        merged_df = (
            df
            # Step 1: Merge particle info to hits
            .Define("particle_merge_result",
                "MergeParticlesToHits("
                "particle_id, Part_barcode, "
                "Part_pt, Part_eta, Part_px, Part_py, "
                "Part_vx, Part_vy, Part_vz, "
                "Part_radius, Part_status, Part_charge, "
                "Part_pdg_id, Part_passed, "
                "Part_vProdNin, Part_vProdNout, Part_vProdStatus, Part_vProdBarcode, "
                "Part_barcode)")  # keep barcode for primary flag
            
            # Step 2: Extract merged particle data
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

            # Step 3: Count hits per particle
            .Define("hit_particle_nhits", "CountHitsPerParticle(hit_particle_id)")
            
            # Step 4: Calculate distance from production vertex
            .Define("hit_distance_from_production", "CalculateDistanceFromProduction(SPx, SPy, SPz, hit_particle_vx, hit_particle_vy, hit_particle_vz)")
            # .Define("hit_particle_nhits", "Take(particle_nhits, valid_hit_indices)")
            # .Define("hit_distance_from_production", "Take(hit_R, valid_hit_indices)")
        )
        
        self.perf.print_stats("Particle-Hit Merging")
        return merged_df

    def perform_cluster_matching(self, df: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Perform cluster matching using hash maps"""
      
        matched_df = (df
            # Step 1: cluster matching
            .Define("cluster_matching_result", 
                   "ClusterMatching(SPCL1_index, SPCL2_index, CLindex, CLparticleLink_barcode)")
            
            # Step 2: Extract matching results
            .Define("matched_cl1_pos", "std::get<0>(cluster_matching_result)")
            .Define("matched_cl2_pos", "std::get<1>(cluster_matching_result)")
            .Define("valid_hit_indices", "std::get<2>(cluster_matching_result)")
            .Define("particle_id_1", "std::get<3>(cluster_matching_result)")
            .Define("particle_id_2", "std::get<4>(cluster_matching_result)")
            .Define("particle_id", "ROOT::VecOps::Where(particle_id_1 == particle_id_2, particle_id_1, -1)")
            
            # Step 3: Create matched hit data using Take for vectorized access
            .Define("hit_id", "Take(SPindex, valid_hit_indices)")
            .Define("hit_x", "Take(SPx, valid_hit_indices)")
            .Define("hit_y", "Take(SPy, valid_hit_indices)")
            .Define("hit_z", "Take(SPz, valid_hit_indices)")
            .Define("hit_r", "CalculateR(hit_x, hit_y)")
            .Define("hit_phi", "CalculatePhi(hit_x, hit_y)")
            .Define("hit_eta", "CalculateEta(hit_r, hit_z)")
            .Define("hit_radius", "Take(SPradius, valid_hit_indices)")
            .Define("hit_overlap", "Take(SPisOverlap, valid_hit_indices)")
            
            # Hit detector information
            .Define("hit_barrel_endcap", "Take(CLbarrel_endcap, matched_cl1_pos)")
            .Define("hit_layer_disk", "Take(CLlayer_disk, matched_cl1_pos)")
            .Define("hit_eta_module", "Take(CLeta_module, matched_cl1_pos)")
            .Define("hit_phi_module", "Take(CLphi_module, matched_cl1_pos)")
            .Define("hit_hardware", "Take(CLhardware, matched_cl1_pos)")
            .Define("hit_module_id", "Take(CLmoduleID, matched_cl1_pos)")
            
            # Step 4: Create matched cluster 1 data
            .Define("hit_cluster_id_1", "Take(CLindex, matched_cl1_pos)")
            .Define("hit_cluster_x_1", "Take(CLx, matched_cl1_pos)")
            .Define("hit_cluster_y_1", "Take(CLy, matched_cl1_pos)")
            .Define("hit_cluster_z_1", "Take(CLz, matched_cl1_pos)")
            .Define("hit_cluster_norm_x_1", "Take(CLnorm_x, matched_cl2_pos)")
            .Define("hit_cluster_norm_y_1", "Take(CLnorm_y, matched_cl2_pos)")
            .Define("hit_cluster_norm_z_1", "Take(CLnorm_z, matched_cl2_pos)")
            .Define("hit_cluster_size_1", "Take(CLside, matched_cl2_pos)")
            .Define("hit_cluster_r_1", "CalculateR(hit_cluster_x_1, hit_cluster_y_1)")
            .Define("hit_cluster_phi_1", "CalculatePhi(hit_cluster_x_1, hit_cluster_y_1)")
            .Define("hit_cluster_eta_1", "CalculateEta(hit_cluster_r_1, hit_cluster_z_1)")
            .Define("hit_cluster_charge_count_1", "Take(CLcharge_count, matched_cl2_pos)")
            .Define("hit_cluster_count_1", "Take(CLpixel_count, matched_cl2_pos)")
            .Define("hit_cluster_localDir0_1", "Take(CLloc_direction1, matched_cl2_pos)")
            .Define("hit_cluster_localDir1_1", "Take(CLloc_direction2, matched_cl2_pos)")
            .Define("hit_cluster_localDir2_1", "Take(CLloc_direction3, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir0_1", "Take(CLJan_loc_direction1, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir1_1", "Take(CLJan_loc_direction2, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir2_1", "Take(CLJan_loc_direction3, matched_cl2_pos)")
            .Define("hit_cluster_loc_eta_1", "Take(CLloc_eta, matched_cl1_pos)")
            .Define("hit_cluster_loc_phi_1", "Take(CLloc_phi, matched_cl1_pos)")
            .Define("hit_cluster_glob_eta_1", "Take(CLglob_eta, matched_cl1_pos)")
            .Define("hit_cluster_glob_phi_1", "Take(CLglob_phi, matched_cl1_pos)")
            .Define("hit_cluster_barrel_endcap_1", "Take(CLbarrel_endcap, matched_cl1_pos)")
            .Define("hit_cluster_layer_disk_1", "Take(CLlayer_disk, matched_cl1_pos)")
            .Define("hit_cluster_eta_module_1", "Take(CLeta_module, matched_cl1_pos)")
            .Define("hit_cluster_phi_module_1", "Take(CLphi_module, matched_cl1_pos)")
            .Define("hit_cluster_hardware_1", "Take(CLhardware, matched_cl1_pos)")
            .Define("hit_cluster_module_id_1", "Take(CLmoduleID, matched_cl1_pos)")
            
            # Step 5: Create matched cluster 2 data
            .Define("hit_cluster_id_2", "Take(CLindex, matched_cl2_pos)")
            .Define("hit_cluster_x_2", "Take(CLx, matched_cl2_pos)")
            .Define("hit_cluster_y_2", "Take(CLy, matched_cl2_pos)")
            .Define("hit_cluster_z_2", "Take(CLz, matched_cl2_pos)")
            .Define("hit_cluster_norm_x_2", "Take(CLnorm_x, matched_cl2_pos)")
            .Define("hit_cluster_norm_y_2", "Take(CLnorm_y, matched_cl2_pos)")
            .Define("hit_cluster_norm_z_2", "Take(CLnorm_z, matched_cl2_pos)")
            .Define("hit_cluster_size_2", "Take(CLside, matched_cl2_pos)")
            .Define("hit_cluster_r_2", "CalculateR(hit_cluster_x_2, hit_cluster_y_2)")
            .Define("hit_cluster_phi_2", "CalculatePhi(hit_cluster_x_2, hit_cluster_y_2)")
            .Define("hit_cluster_eta_2", "CalculateEta(hit_cluster_r_2, hit_cluster_z_2)")
            .Define("hit_cluster_charge_count_2", "Take(CLcharge_count, matched_cl2_pos)")
            .Define("hit_cluster_count_2", "Take(CLpixel_count, matched_cl2_pos)")
            .Define("hit_cluster_loc_eta_2", "Take(CLloc_eta, matched_cl2_pos)")
            .Define("hit_cluster_loc_phi_2", "Take(CLloc_phi, matched_cl2_pos)")
            .Define("hit_cluster_localDir0_2", "Take(CLloc_direction1, matched_cl2_pos)")
            .Define("hit_cluster_localDir1_2", "Take(CLloc_direction2, matched_cl2_pos)")
            .Define("hit_cluster_localDir2_2", "Take(CLloc_direction3, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir0_2", "Take(CLJan_loc_direction1, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir1_2", "Take(CLJan_loc_direction2, matched_cl2_pos)")
            .Define("hit_cluster_lengthDir2_2", "Take(CLJan_loc_direction3, matched_cl2_pos)")
            .Define("hit_cluster_glob_eta_2", "Take(CLglob_eta, matched_cl2_pos)")
            .Define("hit_cluster_glob_phi_2", "Take(CLglob_phi, matched_cl2_pos)")
            .Define("hit_cluster_eta_angle_2", "Take(CLeta_angle, matched_cl2_pos)")
            .Define("hit_cluster_phi_angle_2", "Take(CLphi_angle, matched_cl2_pos)")
            .Define("hit_cluster_barrel_endcap_2", "Take(CLbarrel_endcap, matched_cl2_pos)")
            .Define("hit_cluster_layer_disk_2", "Take(CLlayer_disk, matched_cl2_pos)")
            .Define("hit_cluster_eta_module_2", "Take(CLeta_module, matched_cl1_pos)")
            .Define("hit_cluster_phi_module_2", "Take(CLphi_module, matched_cl1_pos)")
            .Define("hit_cluster_hardware_2", "Take(CLhardware, matched_cl1_pos)")
            .Define("hit_cluster_module_id_2", "Take(CLmoduleID, matched_cl1_pos)")
        )
        
        self.perf.print_stats("Cluster Matching")
        return matched_df
    
    def add_derived_features(self, df: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Add all derived physics quantities using vectorized operations"""
        
        physics_df = (df
            # Cluster separation metrics
            .Define("hit_cluster_separation", 
                   "CalculateClusterSeparation(hit_cluster_x_1, hit_cluster_y_1, hit_cluster_z_1, "
                   "                          hit_cluster_x_2, hit_cluster_y_2, hit_cluster_z_2)")
            
            # Combined cluster properties
            .Define("hit_cluster_charge_total", "hit_cluster_charge_count_1 + hit_cluster_charge_count_2")
            .Define("hit_cluster_pixel_total", "hit_cluster_count_1 + hit_cluster_count_2")
            .Define("hit_cluster_charge_mean", "(hit_cluster_charge_count_1 + hit_cluster_charge_count_2) / 2.0f")
            
            # Cluster charge asymmetry (safer calculation)
            .Define("hit_cluster_charge_asymmetry", 
                   "abs(hit_cluster_charge_count_1 - hit_cluster_charge_count_2) / (hit_cluster_charge_count_1 + hit_cluster_charge_count_2 + 1e-6f)")
            
            # Layer and detector region information
            .Define("hit_same_layer", "hit_cluster_layer_disk_1 == hit_cluster_layer_disk_2")
            .Define("hit_barrel_region", "hit_cluster_barrel_endcap_1 == 0")
            .Define("hit_endcap_region", "hit_cluster_barrel_endcap_1 != 0")
            
            # Hit-level derived quantities
            .Define("hit_is_noise", "hit_particle_id == 0")
            .Define("hit_is_signal", "hit_particle_id > 0")
            
            # Event-level quantities
            .Define("n_hits_total", "static_cast<int>(hit_id.size())")
            .Define("n_hits_signal", "static_cast<int>(Sum(hit_is_signal))")
            .Define("n_hits_noise", "static_cast<int>(Sum(hit_is_noise))")
            .Define("signal_fraction", "n_hits_total > 0 ? static_cast<float>(n_hits_signal) / n_hits_total : 0.0f")
            
            # Detector region statistics
            .Define("n_hits_barrel", "n_hits_total > 0 ? static_cast<int>(Sum(hit_barrel_region)) : 0")
            .Define("n_hits_same_layer", "n_hits_total > 0 ? static_cast<int>(Sum(hit_same_layer)) : 0")
            .Define("barrel_fraction", "n_hits_total > 0 ? static_cast<float>(n_hits_barrel) / n_hits_total : 0.0f")
            
            # Global event statistics
            .Define("total_cluster_charge", "n_hits_total > 0 ? Sum(hit_cluster_charge_total) : 0.0f")
            .Define("mean_cluster_charge", "n_hits_total > 0 ? Mean(hit_cluster_charge_total) : 0.0f")
            .Define("mean_hit_radius", "n_hits_total > 0 ? Mean(hit_r) : 0.0f")
            .Define("mean_cluster_separation", "n_hits_total > 0 ? Mean(hit_cluster_separation) : 0.0f")
        )
        
        self.perf.print_stats("Derived Feature Calculation")
        return physics_df

    def build_track_edges(self, df: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Build track edges between hits on the same particle"""
        
        track_df = (df
            # Filter signal hits only for track building
            .Define("signal_hit_mask", "hit_particle_id > 0")
            
            # Build track edges using particle and module information
            .Define("track_edge_result",
                   "BuildTrackEdges(hit_id, hit_particle_id, hit_distance_from_production, "
                   "hit_barrel_endcap, hit_layer_disk, hit_eta_module, hit_phi_module, hit_hardware)")
            
            # Extract track edges
            .Define("track_edge_from", "std::get<0>(track_edge_result)")
            .Define("track_edge_to", "std::get<1>(track_edge_result)")
            
            # Calculate number of track edges
            .Define("n_track_edges", "static_cast<int>(track_edge_from.size())")
        )
        
        self.perf.print_stats("Track Edge Building")
        return track_df

    def apply_quality_cuts(self, df: ROOT.RDataFrame) -> ROOT.RDataFrame:
        """Apply quality cuts and filters"""
        
        filtered_df = (df
            # Basic event filters
            .Filter("n_hits_total > 0", "Events with hits")
            .Filter("n_hits_signal > 0", "Events with signal hits")
            
            # Physics cuts
            .Define("quality_cuts", "ApplyQualityCuts(hit_particle_pt, hit_particle_eta, hit_particle_charge)")
            .Define("n_hits_passing_cuts", "n_hits_total > 0 ? static_cast<int>(Sum(quality_cuts)) : 0")
            
            # Additional quality filters
            .Filter("mean_cluster_separation > 0.01f", "Minimum cluster separation")
            .Filter("total_cluster_charge < 1e6f", "Maximum charge cut")
            .Filter("n_track_edges > 0", "Events with track edges")
            
            # Count events passing all cuts
            .Define("passes_all_cuts", "1")
        )
        
        self.perf.print_stats("Quality Cuts and Filtering")
        return filtered_df

    def create_summary_statistics(self, df: ROOT.RDataFrame) -> Dict:
        """Generate comprehensive summary statistics"""
        
        print("Generating summary statistics...")
        
        # Take a snapshot of key statistics
        stats_df = (df
            .Define("event_stats", "1")  # Just to trigger evaluation
        )
        
        # Get basic counts
        n_events = stats_df.Count().GetValue()
        
        if n_events == 0:
            return {"n_events": 0, "warning": "No events passed filters"}
        
        # Calculate summary statistics using RDataFrame aggregations
        summary_stats = {
            "n_events": n_events,
            "mean_hits_per_event": stats_df.Mean("n_hits_total").GetValue(),
            "mean_signal_hits_per_event": stats_df.Mean("n_hits_signal").GetValue(),
            "mean_noise_hits_per_event": stats_df.Mean("n_hits_noise").GetValue(),
            "mean_track_edges_per_event": stats_df.Mean("n_track_edges").GetValue(),
            "mean_barrel_fraction": stats_df.Mean("barrel_fraction").GetValue(),
            "mean_signal_fraction": stats_df.Mean("signal_fraction").GetValue(),
            "mean_cluster_charge": stats_df.Mean("mean_cluster_charge").GetValue(),
            "mean_cluster_separation": stats_df.Mean("mean_cluster_separation").GetValue(),
            "total_hits": int(stats_df.Sum("n_hits_total").GetValue()),
            "total_signal_hits": int(stats_df.Sum("n_hits_signal").GetValue()),
            "total_track_edges": int(stats_df.Sum("n_track_edges").GetValue()),
        }
        
        # Print summary
        print(f"\n=== Processing Summary ===")
        for key, value in summary_stats.items():
            if isinstance(value, float):
                print(f"{key}: {value:.3f}")
            else:
                print(f"{key}: {value}")
        print("=" * 30)
        
        self.perf.print_stats("Summary Statistics")
        return summary_stats

    def export_data(self, df: ROOT.RDataFrame, output_path: str, format: str = "csv"):
        """Export processed data to various formats"""
        
        print(f"Exporting data to {output_path} in {format} format...")
        
        # Define the columns to export for graph building
        hit_features = [
            "hit_id", "hit_x", "hit_y", "hit_z", "hit_r", "hit_phi", "hit_eta",
            "hit_barrel_endcap", "hit_layer_disk", "hit_eta_module", "hit_phi_module",
            "hit_hardware", "hit_module_id", "hit_radius", "hit_overlap"
        ]
        
        particle_features = [
            "hit_particle_id", "hit_particle_pt", "hit_particle_eta", "hit_particle_phi",
            "hit_particle_vx", "hit_particle_vy", "hit_particle_vz", "hit_particle_charge",
            "hit_particle_primary", "hit_particle_nhits"
        ]
        
        cluster_features = [
            "hit_cluster_id_1", "hit_cluster_x_1", "hit_cluster_y_1", "hit_cluster_z_1",
            "hit_cluster_r_1", "hit_cluster_phi_1", "hit_cluster_eta_1",
            "hit_cluster_charge_count_1", "hit_cluster_count_1",
            "hit_cluster_id_2", "hit_cluster_x_2", "hit_cluster_y_2", "hit_cluster_z_2",
            "hit_cluster_r_2", "hit_cluster_phi_2", "hit_cluster_eta_2",
            "hit_cluster_charge_count_2", "hit_cluster_count_2"
        ]
        
        derived_features = [
            "hit_cluster_separation", "hit_cluster_charge_total", "hit_cluster_charge_asymmetry",
            "hit_same_layer", "hit_barrel_region", "hit_is_noise", "hit_distance_from_production"
        ]
        
        track_features = [
            "track_edge_from", "track_edge_to"
        ]
        
        event_features = [
            "run_number", "event_number", "n_hits_total", "n_hits_signal", "n_track_edges",
            "signal_fraction", "barrel_fraction", "total_cluster_charge"
        ]
        
        all_features = hit_features + particle_features + cluster_features + derived_features + track_features + event_features

        # Export as numpy arrays
        numpy_data = df.AsNumpy(all_features)
        print("format.lower()",format.lower())
        if format.lower() == "csv":
            # Convert to pandas DataFrame
            df_pd = pd.DataFrame(numpy_data)
            print(df_pd.head(2))
            # Save to CSV
            df_pd.to_csv(output_path, index=False)
            
        elif format.lower() == "root":
            # Export as ROOT file
            df.Snapshot("GNN4ITk_processed", output_path, all_features, opts)
            
        elif format.lower() == "numpy":
            np.savez_compressed(output_path, **numpy_data)
            
        else:
            raise ValueError(f"Unsupported export format: {format}")
            
        self.perf.print_stats(f"Data Export ({format})")

    def process_complete_pipeline(self, output_file: str = None, export_format: str = "numpy") -> Tuple[ROOT.RDataFrame, Dict]:
        """Run the complete processing pipeline with all graph variables"""
        
        print("Starting comprehensive GNN4ITk processing pipeline...")
        
        # Load data
        df = self.load_data()

        # Step 1: Perform cluster matching
        matched_df = self.perform_cluster_matching(df)
        
        # Step 2: Merge particle information to hits
        merged_df = self.merge_particles_to_hits(matched_df)
        
        # Step 3: Add all derived features
        physics_df = self.add_derived_features(merged_df)
        
        # Step 4: Build track edges
        track_df = self.build_track_edges(physics_df)
        
        # Step 5: Apply quality cuts
        final_df = self.apply_quality_cuts(track_df)
        
        # Step 6: Generate summary statistics
        stats = self.create_summary_statistics(final_df)
        
        # Step 7: Export data if requested
        if output_file:
            self.export_data(final_df, output_file, export_format)
        
        return final_df, stats

    def create_sample_event_display(self, df: ROOT.RDataFrame, event_idx: int = 0):
        """Create a sample display of processed event data"""
        
        print(f"\n=== Sample Event Display (Event {event_idx}) ===")
        
        # Convert to numpy for easy inspection
        sample_columns = [
            "run_number", "event_number", "n_hits_total", "n_hits_signal", "n_track_edges",
            "hit_id", "hit_x", "hit_y", "hit_z", "hit_r", "hit_phi", "hit_eta",
            "hit_particle_id", "hit_particle_pt", "hit_cluster_charge_total",
            "track_edge_from", "track_edge_to"
        ]
        
        try:
            numpy_data = df.AsNumpy(sample_columns)
            
            if len(numpy_data['run_number']) > event_idx:
                print(f"Run: {numpy_data['run_number'][event_idx]}")
                print(f"Event: {numpy_data['event_number'][event_idx]}")
                print(f"Total Hits: {numpy_data['n_hits_total'][event_idx]}")
                print(f"Signal Hits: {numpy_data['n_hits_signal'][event_idx]}")
                print(f"Track Edges: {numpy_data['n_track_edges'][event_idx]}")
                
                # Show first few hits
                hit_ids = numpy_data['hit_id'][event_idx]
                if len(hit_ids) > 0:
                    print(f"\nFirst 5 hits:")
                    for i in range(min(5, len(hit_ids))):
                        print(f"  Hit {hit_ids[i]}: (x={numpy_data['hit_x'][event_idx][i]:.2f}, "
                              f"y={numpy_data['hit_y'][event_idx][i]:.2f}, "
                              f"z={numpy_data['hit_z'][event_idx][i]:.2f})")
                
                # Show first few track edges
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
    
    parser = argparse.ArgumentParser(description='Process GNN4ITk ROOT NTuple data with complete graph variables')
    parser.add_argument('--input_file', '-i', help='Input ROOT file path')
    parser.add_argument('--output', '-o', help='Output file path')
    parser.add_argument('--format', '-f', choices=['csv', 'root', 'numpy'], default='csv',
                       help='Output format (default: numpy)')
    parser.add_argument('--threads', '-t', type=int, default=4, 
                       help='Number of threads for processing (default: 4)')
    parser.add_argument('--show-sample', action='store_true',
                       help='Display sample event data')
    
    args = parser.parse_args()
    
    # Set number of threads
    ROOT.EnableImplicitMT(args.threads)
    print(f"Using {args.threads} threads for processing")
    
    try:
        # Initialize processor
        processor = ComprehensiveGNN4ITkProcessor(args.input_file)
        
        # Run complete pipeline
        final_df, stats = processor.process_complete_pipeline(args.output, args.format)
        
        # Show sample event if requested
        if args.show_sample:
            processor.create_sample_event_display(final_df)
        
        print(f"\nProcessing completed successfully!")
        print(f"Processed {stats['n_events']} events")
        print(f"Total hits: {stats['total_hits']}")
        print(f"Total track edges: {stats['total_track_edges']}")
        
        if args.output:
            print(f"Results exported to: {args.output}")
            
    except Exception as e:
        print(f"Error during processing: {e}")
        import traceback
        traceback.print_exc()
        return 1
    
    return 0


if __name__ == "__main__":
    exit(main())
#include <ROOT/RDataFrame.hxx>
#include <ROOT/RVec.hxx>
#include <TFile.h>
#include <TStopwatch.h>
#include <iostream>
#include <unordered_map>
#include <vector>

using namespace ROOT;
using namespace ROOT::VecOps;

// Performance monitoring class
class PerformanceMonitor {
private:
    TStopwatch timer;
    size_t initial_memory;
    
public:
    PerformanceMonitor() { 
        timer.Start(); 
        initial_memory = GetCurrentMemoryUsage();
    }
    
    void PrintStats(const std::string& operation) {
        timer.Stop();
        size_t current_memory = GetCurrentMemoryUsage();
        std::cout << "=== Performance Report: " << operation << " ===" << std::endl;
        std::cout << "Time: " << timer.RealTime() << " seconds" << std::endl;
        std::cout << "Memory Delta: " << (current_memory - initial_memory) / 1024 / 1024 << " MB" << std::endl;
        std::cout << "======================================" << std::endl;
        timer.Reset();
        timer.Start();
        initial_memory = current_memory;
    }
    
private:
    size_t GetCurrentMemoryUsage() {
        // Simple memory usage estimation (platform-specific implementation needed)
        return 0; // Placeholder
    }
};

// Optimized cluster matching function using hash maps for O(1) lookup
auto OptimizedClusterMatching = [](const RVec<int>& hit_cl1_idx, const RVec<int>& hit_cl2_idx,
                                   const RVec<int>& cluster_indices) {
    
    // Pre-build hash map for O(1) cluster lookups
    std::unordered_map<int, size_t> cluster_map;
    cluster_map.reserve(cluster_indices.size());
    
    for (size_t i = 0; i < cluster_indices.size(); ++i) {
        cluster_map[cluster_indices[i]] = i;
    }
    
    // Result vectors - reserve space for performance
    RVec<int> matched_cl1_positions, matched_cl2_positions, hit_indices;
    matched_cl1_positions.reserve(hit_cl1_idx.size());
    matched_cl2_positions.reserve(hit_cl1_idx.size());
    hit_indices.reserve(hit_cl1_idx.size());
    
    // Fast cluster matching
    for (size_t hit_idx = 0; hit_idx < hit_cl1_idx.size(); ++hit_idx) {
        auto cl1_it = cluster_map.find(hit_cl1_idx[hit_idx]);
        auto cl2_it = cluster_map.find(hit_cl2_idx[hit_idx]);
        
        if (cl1_it != cluster_map.end() && cl2_it != cluster_map.end()) {
            matched_cl1_positions.push_back(cl1_it->second);
            matched_cl2_positions.push_back(cl2_it->second);
            hit_indices.push_back(hit_idx);
        }
    }
    
    return std::make_tuple(matched_cl1_positions, matched_cl2_positions, hit_indices);
};

// Main processing function
void ProcessGNN4ITkData(const std::string& filename) {
    PerformanceMonitor perf;
    
    // Open file and create RDataFrame
    std::cout << "Opening file: " << filename << std::endl;
    RDataFrame df("GNN4ITk", filename);
    perf.PrintStats("File Opening");
    
    // Enable multi-threading for better performance
    ROOT::EnableImplicitMT(4); // Adjust based on your CPU cores
    
    // Step 1: Create optimized particle DataFrame with renamed columns
    auto particle_df = df.Define("subevent", "Part_event_number")
                         .Define("barcode", "Part_barcode")
                         .Define("px", "Part_px")
                         .Define("py", "Part_py")
                         .Define("pz", "Part_pz")
                         .Define("pt", "Part_pt")
                         .Define("eta", "Part_eta")
                         .Define("vx", "Part_vx")
                         .Define("vy", "Part_vy")
                         .Define("vz", "Part_vz")
                         .Define("radius", "Part_radius")
                         .Define("status", "Part_status")
                         .Define("charge", "Part_charge")
                         .Define("pdgId", "Part_pdg_id")
                         .Define("pass", "Part_passed")
                         .Define("vProdNIn", "Part_vProdNin")
                         .Define("vProdNOut", "Part_vProdNout")
                         .Define("vProdStatus", "Part_vProdStatus")
                         .Define("vProdBarcode", "Part_vProdBarcode");
    
    perf.PrintStats("Particle DataFrame Creation");
    
    // Step 2: Create optimized hit DataFrame
    auto hit_df = df.Define("hit_id", "SPindex")
                    .Define("x", "SPx")
                    .Define("y", "SPy")
                    .Define("z", "SPz")
                    .Define("cluster_index_1", "SPCL1_index")
                    .Define("cluster_index_2", "SPCL2_index")
                    .Define("isOverlap", "SPisOverlap")
                    .Define("hit_radius", "SPradius")
                    .Define("covr", "SPcovr")
                    .Define("covz", "SPcovz")
                    .Define("hl_topstrip", "SPhl_topstrip")
                    .Define("hl_botstrip", "SPhl_botstrip");
    
    perf.PrintStats("Hit DataFrame Creation");
    
    // Step 3: Create optimized cluster DataFrame with all relevant fields
    auto cluster_df = df.Define("cluster_id", "CLindex")
                        .Define("hardware", "CLhardware")
                        .Define("cluster_x", "CLx")
                        .Define("cluster_y", "CLy")
                        .Define("cluster_z", "CLz")
                        .Define("barrel_endcap", "CLbarrel_endcap")
                        .Define("layer_disk", "CLlayer_disk")
                        .Define("eta_module", "CLeta_module")
                        .Define("phi_module", "CLphi_module")
                        .Define("side", "CLside")
                        .Define("module_id", "CLmoduleID")
                        .Define("count", "CLpixel_count")
                        .Define("charge_count", "CLcharge_count")
                        .Define("loc_eta", "CLloc_eta")
                        .Define("loc_phi", "CLloc_phi")
                        .Define("localDir0", "CLloc_direction1")
                        .Define("localDir1", "CLloc_direction2")
                        .Define("localDir2", "CLloc_direction3")
                        .Define("lengthDir0", "CLJan_loc_direction1")
                        .Define("lengthDir1", "CLJan_loc_direction2")
                        .Define("lengthDir2", "CLJan_loc_direction3")
                        .Define("glob_eta", "CLglob_eta")
                        .Define("glob_phi", "CLglob_phi")
                        .Define("eta_angle", "CLeta_angle")
                        .Define("phi_angle", "CLphi_angle")
                        .Define("norm_x", "CLnorm_x")
                        .Define("norm_y", "CLnorm_y")
                        .Define("norm_z", "CLnorm_z")
                        .Define("cl_subevent", "CLparticleLink_eventIndex")
                        .Define("cl_barcode", "CLparticleLink_barcode");
    
    perf.PrintStats("Cluster DataFrame Creation");
    
    // Step 4: Optimized cluster matching using custom function
    auto combined_df = df.Define("cluster_matching_result", 
                                [](const RVec<int>& hit_cl1_idx, const RVec<int>& hit_cl2_idx,
                                   const RVec<int>& cluster_indices) {
                                    return OptimizedClusterMatching(hit_cl1_idx, hit_cl2_idx, cluster_indices);
                                }, {"SPCL1_index", "SPCL2_index", "CLindex"})
                        
                        // Extract matched indices
                        .Define("matched_cl1_pos", "std::get<0>(cluster_matching_result)")
                        .Define("matched_cl2_pos", "std::get<1>(cluster_matching_result)")
                        .Define("valid_hit_indices", "std::get<2>(cluster_matching_result)")
                        
                        // Create matched hit data
                        .Define("matched_hit_ids", "Take(SPindex, valid_hit_indices)")
                        .Define("matched_hit_x", "Take(SPx, valid_hit_indices)")
                        .Define("matched_hit_y", "Take(SPy, valid_hit_indices)")
                        .Define("matched_hit_z", "Take(SPz, valid_hit_indices)")
                        
                        // Create matched cluster data
                        .Define("matched_cl1_x", "Take(CLx, matched_cl1_pos)")
                        .Define("matched_cl1_y", "Take(CLy, matched_cl1_pos)")
                        .Define("matched_cl1_z", "Take(CLz, matched_cl1_pos)")
                        .Define("matched_cl1_charge", "Take(CLcharge_count, matched_cl1_pos)")
                        .Define("matched_cl1_count", "Take(CLpixel_count, matched_cl1_pos)")
                        
                        .Define("matched_cl2_x", "Take(CLx, matched_cl2_pos)")
                        .Define("matched_cl2_y", "Take(CLy, matched_cl2_pos)")
                        .Define("matched_cl2_z", "Take(CLz, matched_cl2_pos)")
                        .Define("matched_cl2_charge", "Take(CLcharge_count, matched_cl2_pos)")
                        .Define("matched_cl2_count", "Take(CLpixel_count, matched_cl2_pos)");
    
    perf.PrintStats("Cluster Matching");
    
    // Step 5: Add derived physics quantities
    auto physics_df = combined_df
        // Hit-level derived quantities
        .Define("hit_r", "sqrt(matched_hit_x*matched_hit_x + matched_hit_y*matched_hit_y)")
        .Define("hit_phi", "atan2(matched_hit_y, matched_hit_x)")
        .Define("hit_theta", "atan2(hit_r, matched_hit_z)")
        .Define("hit_eta", "-log(tan(hit_theta/2.0))")
        
        // Cluster separation metrics
        .Define("cl_separation", "sqrt((matched_cl1_x-matched_cl2_x)*(matched_cl1_x-matched_cl2_x) + "
                                 "(matched_cl1_y-matched_cl2_y)*(matched_cl1_y-matched_cl2_y) + "
                                 "(matched_cl1_z-matched_cl2_z)*(matched_cl1_z-matched_cl2_z))")
        
        // Combined cluster properties
        .Define("total_cluster_charge", "matched_cl1_charge + matched_cl2_charge")
        .Define("total_cluster_count", "matched_cl1_count + matched_cl2_count");
    
    perf.PrintStats("Physics Quantities");
    
    // Step 6: Performance optimized filtering and selection
    auto filtered_df = physics_df
        .Filter("matched_hit_ids.size() > 0", "Non-empty events")
        .Filter("total_cluster_charge > 0", "Valid cluster charge")
        .Define("n_valid_hits", "matched_hit_ids.size()");
    
    perf.PrintStats("Filtering");
    
    // Step 7: Create summary statistics
    auto stats = filtered_df
        .Define("event_summary", [](ULong64_t entry, const RVec<int>& hit_ids, 
                                   const RVec<float>& charges, const RVec<float>& separations) {
            return std::make_tuple(entry, hit_ids.size(), 
                                 charges.size() > 0 ? *std::max_element(charges.begin(), charges.end()) : 0.0f,
                                 separations.size() > 0 ? Mean(separations) : 0.0f);
        }, {"rdfentry_", "matched_hit_ids", "total_cluster_charge", "cl_separation"});
    
    perf.PrintStats("Summary Statistics");
    
    // Step 8: Demonstrate data access patterns
    std::cout << "\n=== Data Access Examples ===" << std::endl;
    
    // Count total events
    auto total_events = stats.Count();
    std::cout << "Total events processed: " << *total_events << std::endl;
    
    // Average hits per event
    auto avg_hits = stats.Mean("n_valid_hits");
    std::cout << "Average hits per event: " << *avg_hits << std::endl;
    
    // Maximum cluster charge
    auto max_charge = filtered_df.Max("total_cluster_charge");
    std::cout << "Maximum total cluster charge: " << *max_charge << std::endl;
    
    perf.PrintStats("Data Access");
    
    // Optional: Save processed data
    // filtered_df.Snapshot("ProcessedData", "output.root", 
    //                     {"matched_hit_x", "matched_hit_y", "matched_hit_z", 
    //                      "matched_cl1_charge", "matched_cl2_charge", "cl_separation"});
    
    std::cout << "\n=== Processing Complete ===" << std::endl;
}

// Example usage and performance testing
int main(int argc, char* argv[]) {
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " <ntuple_file.root>" << std::endl;
        return 1;
    }
    
    std::string filename = argv[1];
    
    try {
        ProcessGNN4ITkData(filename);
    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}

// Additional utility functions for advanced analysis
namespace AnalysisUtils {
    
    // Function to create physics-based cuts
    auto CreatePhysicsCuts = [](float min_pt = 1.0, float max_eta = 2.5, int min_clusters = 2) {
        return [min_pt, max_eta, min_clusters](const RVec<float>& pt, const RVec<float>& eta, 
                                              const RVec<int>& n_clusters) {
            RVec<bool> cuts;
            cuts.reserve(pt.size());
            for (size_t i = 0; i < pt.size(); ++i) {
                cuts.push_back(pt[i] > min_pt && abs(eta[i]) < max_eta && n_clusters[i] >= min_clusters);
            }
            return cuts;
        };
    };
    
    // Function for detector geometry calculations
    auto CalculateDetectorRegion = [](const RVec<int>& barrel_endcap, const RVec<int>& layer_disk) {
        RVec<std::string> regions;
        regions.reserve(barrel_endcap.size());
        for (size_t i = 0; i < barrel_endcap.size(); ++i) {
            if (barrel_endcap[i] == 0) {
                regions.push_back("Barrel_L" + std::to_string(layer_disk[i]));
            } else {
                regions.push_back("Endcap_D" + std::to_string(layer_disk[i]));
            }
        }
        return regions;
    };
}
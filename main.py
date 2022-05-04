import os
from pymetabolomics import *

mzML_dir = "/home/axel/Nextcloud/workspace/MetabolomicsWorkflowMayer/mzML"
mzML_dir = "/home/axel/Nextcloud/workspace/Tests/CentogeneWorkflowTest/mzML"
results = "results"
Helper.reset_directory(results)

FeatureFinderMetabo.run(mzML_dir, os.path.join(results, "FFM"),
                        {"noise_threshold_int": 10000.0,
                        "mass_error_ppm": 10.0,
                        "remove_single_traces": "true"})

MapAligner.run(os.path.join(results, "FFM"), os.path.join(results, "FFM_aligned"), os.path.join(results, "Trafo"),
                {"max_num_peaks_considered": -1,
                "superimposer:mz_pair_max_distance": 0.05,
                "pairfinder:distance_MZ:max_difference": 10.0,
                "pairfinder:distance_MZ:unit": "ppm"})

FeatureLinker.run(os.path.join(results, "FFM_aligned"), os.path.join(results, "ConsensusMaps", "FFM.consensusXML"))

DataFrames.create_consensus_table(os.path.join(results, "ConsensusMaps", "FFM.consensusXML"), 
                                os.path.join(results, "ConsensusMaps", "FFM.tsv"))

FeatureMapHelper.split_consensus_map(os.path.join(results, "ConsensusMaps", "FFM.consensusXML"),
                                    os.path.join(results, "ConsensusMaps", "FFM_complete.consensusXML"),
                                    os.path.join(results, "ConsensusMaps", "FFM_missing.consensusXML"))
FeatureMapHelper.consensus_to_feature_maps(os.path.join(results, "ConsensusMaps", "FFM_complete.consensusXML"),
                                        os.path.join(results, "FFM"),
                                        os.path.join(results, "FFM_complete"))

MapAligner.run(mzML_dir, os.path.join(results, "MzML_aligned"), os.path.join(results, "Trafo"))

library = FeatureFinderMetaboIdent.load_library(os.path.join(results, "ConsensusMaps", "FFM.consensusXML"), os.path.join(results, "FFMID_library.tsv"))
FeatureFinderMetaboIdent.run(os.path.join(results, "MzML_aligned"), os.path.join(results, "FFMID"), library)

FeatureMapHelper.merge_feature_maps(os.path.join(results, "FeatureMaps_merged"), os.path.join(results, "FFM_complete"), os.path.join(results, "FFMID"))

# # AdductDecharger (optional)
# decharged_dir = reset_directory(os.path.join(result_dir, "FFMID_decharged"))
MetaboliteAdductDecharger.run(os.path.join(results, "FeatureMaps_merged"), os.path.join(results, "FeatureMaps_decharged"),
                        {"potential_adducts": [b"H:+:0.5",b"Na:+:0.3", b"H-1O-1:+:0.2"],
                        "charge_min": 1,
                        "charge_max": 3,
                        "max_neutrals": 2})

# # IDMapper (optional)
# mapped_dir = reset_directory(os.path.join(result_dir, "FFMID_id_mapped"))
# IDMapper.run(mzML_aligned_dir, decharged_dir, mapped_dir)
import os
from core import *
from helpers import *
from dataframes import *

mzML_dir = "/home/axel/Nextcloud/workspace/MetabolomicsWorkflowMayer/mzML"
mzML_dir = "/home/axel/Nextcloud/workspace/Tests/CentogeneWorkflowTest/mzML"

results = Helper.reset_directory("results")
interim = Helper.reset_directory(os.path.join(results, "interim"))

FeatureFinderMetabo.run(mzML_dir, os.path.join(interim, "FFM"),
                        {"noise_threshold_int": 10000.0,
                        "mass_error_ppm": 10.0,
                         "remove_single_traces": "true"})

MapAligner.run(os.path.join(interim, "FFM"), os.path.join(interim, "FFM_aligned"), os.path.join(interim, "Trafo"),
               {"max_num_peaks_considered": -1,
                "superimposer:mz_pair_max_distance": 0.05,
                "pairfinder:distance_MZ:max_difference": 10.0,
                "pairfinder:distance_MZ:unit": "ppm"})

FeatureLinker.run(os.path.join(interim, "FFM_aligned"),
                  os.path.join(interim,  "FFM.consensusXML"))

DataFrames.create_consensus_table(os.path.join(interim,  "FFM.consensusXML"),
                                  os.path.join(interim,  "FFM.tsv"))

FeatureMapHelper.split_consensus_map(os.path.join(interim,  "FFM.consensusXML"),
                                     os.path.join(
                                         interim,  "FFM_complete.consensusXML"),
                                     os.path.join(interim,  "FFM_missing.consensusXML"))
FeatureMapHelper.consensus_to_feature_maps(os.path.join(interim,  "FFM_complete.consensusXML"),
                                           os.path.join(interim, "FFM"),
                                           os.path.join(interim, "FFM_complete"))

MapAligner.run(mzML_dir, os.path.join(interim, "MzML_aligned"),
               os.path.join(interim, "Trafo"))

library = FeatureFinderMetaboIdent.load_library(os.path.join(interim,  "FFM.consensusXML"),
                                                os.path.join(interim, "FFMID_library.tsv"))
FeatureFinderMetaboIdent.run(os.path.join(
    interim, "MzML_aligned"), os.path.join(interim, "FFMID"), library)

FeatureMapHelper.merge_feature_maps(os.path.join(interim, "FeatureMaps_merged"), os.path.join(
    interim, "FFM_complete"), os.path.join(interim, "FFMID"))

MetaboliteAdductDecharger.run(os.path.join(interim, "FeatureMaps_merged"), os.path.join(interim, "FeatureMaps_decharged"),
                              {"potential_adducts": [b"H:+:0.5", b"Na:+:0.3", b"H-1O-1:+:0.2"],
                               "charge_min": 1,
                               "charge_max": 3,
                               "max_neutrals": 2})

MapID.run(os.path.join(interim, "MzML_aligned"), os.path.join(
    interim, "FeatureMaps_decharged"), os.path.join(interim, "FeatureMaps_ID_mapped"))

FeatureLinker.run(os.path.join(interim, "FeatureMaps_merged"),
                  os.path.join(interim, "FeatureMatrix.consensusXML"))

DataFrames.create_consensus_table(os.path.join(
    interim, "FeatureMatrix.consensusXML"), os.path.join(results, "FeatureMatrix.tsv"))

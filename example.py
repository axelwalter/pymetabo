import os
from core import *
from helpers import *
from dataframes import *
from sirius import *

mzML_dir = "../data/mzML"

results = Helper().reset_directory("results")
interim = Helper().reset_directory(os.path.join(results, "interim"))

PrecursorCorrector().to_highest_intensity(mzML_dir, os.path.join(interim, "mzML_PCpeak"))

FeatureFinderMetabo().run(os.path.join(interim, "mzML_PCpeak"), os.path.join(interim, "FFM"),
                        {"noise_threshold_int": 10000.0,
                        "mass_error_ppm": 10.0,
                        "remove_single_traces": "true"})

MapAligner().run(os.path.join(interim, "FFM"), os.path.join(interim, "FFM_aligned"),
                os.path.join(interim, "Trafo"),
                {"max_num_peaks_considered": -1,
                "superimposer:mz_pair_max_distance": 0.05,
                "pairfinder:distance_MZ:max_difference": 10.0,
                "pairfinder:distance_MZ:unit": "ppm"})

PrecursorCorrector().to_nearest_feature(os.path.join(interim, "mzML_PCpeak"), 
                                        os.path.join(interim, "mzML_PCfeature"),
                                        os.path.join(interim, "FFM"))

MapAligner().run(os.path.join(interim, "mzML_PCfeature"), os.path.join(interim, "MzML_aligned"),
                os.path.join(interim, "Trafo"))

# FeatureLinker().run(os.path.join(interim, "FFM_aligned"),
#                 os.path.join(interim,  "FFM.consensusXML"))

# DataFrames().create_consensus_table(os.path.join(
#     interim, "FFM.consensusXML"), os.path.join(interim, "FFM_consensus.tsv"))

# FeatureMapHelper().FFMID_libraries_for_missing_features(os.path.join(interim,  "FFM.consensusXML"),
#                                                     os.path.join(interim,  "FFMID_libraries"))

# FeatureFinderMetaboIdent().run(os.path.join(interim, "MzML_aligned"),
#                             os.path.join(interim,  "FFMID"),
#                             os.path.join(interim,  "FFMID_libraries"))

# FeatureMapHelper().merge_feature_maps(os.path.join(interim, "FeatureMaps_merged"), os.path.join(
#     interim, "FFM"), os.path.join(interim, "FFMID"))

MetaboliteAdductDecharger().run(os.path.join(interim, "FFM"), 
                                os.path.join(interim, "FeatureMaps_decharged"),
                                {"potential_adducts": [b"H:+:0.4",b"Na:+:0.2",b"NH4:+:0.2", 
                                                        b"H-1O-1:+:0.1", b"H-3O-2:+:0.1"],
                                "charge_min": 1,
                                "charge_max": 3,
                                "max_neutrals": 2})

MapID().run(os.path.join(interim, "MzML_aligned"), os.path.join(
    interim, "FeatureMaps_decharged"), os.path.join(interim, "FeatureMaps_ID_mapped"))

FeatureLinker().run(os.path.join(interim, "FeatureMaps_ID_mapped"),
                os.path.join(interim, "FeatureMatrix.consensusXML"))

DataFrames().create_consensus_table(os.path.join(
    interim, "FeatureMatrix.consensusXML"), os.path.join(results, "FeatureMatrix.tsv"))

Sirius().run(os.path.join(interim, "MzML_aligned"), os.path.join(interim, "FeatureMaps_decharged"),
            os.path.join(results, "Sirius"), params={
                "preprocessing:filter_by_num_masstraces": 2,
                "preprocessing:feature_only": "true",
                "sirius:profile": "orbitrap",
                "sirius:db": "none",
                "sirius:ions_considered": "[M+H]+, [M-H2O+H]+, [M+Na]+, [M+NH4]+",
                "sirius:elements_enforced": "CHN[15]OS[4]Cl[2]P[2]",
                "project:processors": 2,
                "fingerid:db": "BIO"
            })
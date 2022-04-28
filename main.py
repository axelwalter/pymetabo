import os
from modules.helpers import reset_directory, load_feature_maps
import modules.FeatureFinderMetabo as FeatureFinderMetabo
import modules.MapAligner as MapAligner
import modules.FeatureLinker as FeatureLinker
import modules.DataFrames as DataFrames
import modules.FeatureFinderMetaboIdent as FeatureFinderMetaboIdent
import modules.AdductDecharger as AdductDecharger
import modules.IDMapper as IDMapper
import modules.FeatureMapFunctions as FeatureMapFunctions

mzML_dir = "/home/axel/Nextcloud/workspace/Tests/CentogeneWorkflowTest/mzML"
result_dir = "results"

reset_directory(result_dir)

# FeatureFinderMetabo
ffm_dir = reset_directory(os.path.join(result_dir, "FFM"))
for mzML_file in os.listdir(mzML_dir):
    if not mzML_file.endswith("mzML"):
        continue
    FeatureFinderMetabo.run(os.path.join(mzML_dir, mzML_file), os.path.join(ffm_dir, mzML_file[:-4] + "featureXML"))

# MapAligner for Feature Maps
ffm_aligned_dir = reset_directory(os.path.join(result_dir, "FFM_aligned"))
trafo_dir = reset_directory(os.path.join(result_dir, "Trafo"))
MapAligner.align_feature_maps(load_feature_maps(ffm_dir), ffm_aligned_dir, trafo_dir)

# FeatureLinker
ffm_consensus_dir = reset_directory(os.path.join(result_dir, "FFM_consensus"))
FeatureLinker.run(load_feature_maps(ffm_aligned_dir), os.path.join(ffm_consensus_dir, "FFM.consensusXML"))
DataFrames.create_table_FFM(os.path.join(ffm_consensus_dir, "FFM.consensusXML"), os.path.join(ffm_consensus_dir, "FFM.tsv"))

# Reconstruct FeatureMaps from ConsensusMap with Features that are detected in ALL samples
# only the missing features will be re-quantified
FeatureMapFunctions.split_consensus_map(os.path.join(ffm_consensus_dir, "FFM.consensusXML"),
                                        os.path.join(ffm_consensus_dir, "FFM_complete.consensusXML"),
                                        os.path.join(ffm_consensus_dir, "FFM_missing.consensusXML"))
ffm_complete = reset_directory(os.path.join(result_dir, "FFM_complete"))
FeatureMapFunctions.consensus_to_feature_maps(os.path.join(ffm_consensus_dir, "FFM_complete.consensusXML"), 
                                            load_feature_maps(os.path.join(result_dir, "FFM_aligned")), ffm_complete)
# Map Aligner for mzML files
mzML_aligned_dir = reset_directory(os.path.join(result_dir, "MzML_aligned"))
MapAligner.align_mzML(mzML_dir, mzML_aligned_dir, trafo_dir)

# FeatureFinderMetaboIdent (Requantification)
ffmid_lib_dir = reset_directory(os.path.join(result_dir, "FFMID_library"))
ffmid_dir = reset_directory(os.path.join(result_dir, "FFMID"))
library = DataFrames.load_FFMID_library(os.path.join(ffm_consensus_dir, "FFM_missing.consensusXML"),
                                        os.path.join(ffmid_lib_dir, "library_requantification.tsv"))
FeatureFinderMetaboIdent.run(ffmid_dir, mzML_aligned_dir, library)

# Merge FeatureMaps from FFM and FFMID
fm_merged_dir = reset_directory(os.path.join(result_dir, "FeatureMaps_merged"))
FeatureMapFunctions.merge_feature_maps(fm_merged_dir, ffm_complete, ffmid_dir)

# AdductDecharger (optional)
decharged_dir = reset_directory(os.path.join(result_dir, "FFMID_decharged"))
AdductDecharger.run(ffmid_dir, decharged_dir)

# IDMapper (optional)
mapped_dir = reset_directory(os.path.join(result_dir, "FFMID_id_mapped"))
IDMapper.run(mzML_aligned_dir, decharged_dir, mapped_dir)
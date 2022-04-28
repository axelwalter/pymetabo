import os
from pyopenms import *

def split_consensus_map(consensusXML_file, consensusXML_complete_file, consensusXML_missing_file):
    consensus_map = ConsensusMap()
    ConsensusXMLFile().load(consensusXML_file, consensus_map)

    headers = consensus_map.getColumnHeaders()

    complete = ConsensusMap(consensus_map)
    complete.clear(False)
    missing = ConsensusMap(consensus_map)
    missing.clear(False)

    for cf in consensus_map:
        if len(cf.getFeatureList()) < len(headers):
            missing.push_back(cf)
        else:
            complete.push_back(cf)
            
    ConsensusXMLFile().store(consensusXML_complete_file, complete)
    ConsensusXMLFile().store(consensusXML_missing_file, missing)

# takes a (filtered) ConsensusMap and reconstructs the FeatureMaps
def consensus_to_feature_maps(consensusXML_file, feature_maps, reconstructed_featureXML_dir):
    consensus_map = ConsensusMap()
    ConsensusXMLFile().load(consensusXML_file, consensus_map)
    to_keep_ids = [item for sublist in [[feature.getUniqueId() for feature in cf.getFeatureList()] for cf in consensus_map] for item in sublist]
    for fm in feature_maps:
        fm_filterd = FeatureMap(fm)
        fm_filterd.clear(False)
        for f in fm:
            if f.getUniqueId() in to_keep_ids:
                fm_filterd.push_back(f)
        FeatureXMLFile().store(os.path.join(reconstructed_featureXML_dir, os.path.basename(fm_filterd.getMetaValue("spectra_data")[0].decode())[:-4] + "featureXML"), fm_filterd)

def merge_feature_maps(fm_merged_dir, featureXML_complete_dir, featureXML_FFMID):
    for file_ffm in os.listdir(featureXML_complete_dir):
        for file_ffmid in os.listdir(featureXML_FFMID):
            if file_ffm == file_ffmid:
                fm_ffm = FeatureMap()
                FeatureXMLFile().load(os.path.join(featureXML_complete_dir, file_ffm), fm_ffm)
                fm_ffmid = FeatureMap()
                FeatureXMLFile().load(os.path.join(featureXML_FFMID, file_ffm), fm_ffmid)
                for f in fm_ffmid:
                    fm_ffm.push_back(f)
                FeatureXMLFile().store(os.path.join(fm_merged_dir, file_ffm), fm_ffm)
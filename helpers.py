import os
import csv
import shutil
from pyopenms import *

class Helper:
    def reset_directory(self, path: str):
        try:
            shutil.rmtree(path)
            os.mkdir(path)
        except OSError:
            os.mkdir(path)
        return path

    def load_feature_maps(self, path: str):
        feature_maps = []
        for file in os.listdir(path):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(path, file), feature_map)
            feature_maps.append(feature_map)
        return feature_maps

class FeatureMapHelper:
    def split_consensus_map(self, consensusXML_file, consensusXML_complete_file="", consensusXML_missing_file=""):
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
        if consensusXML_complete_file.endswith("consensusXML"):
            ConsensusXMLFile().store(consensusXML_complete_file, complete)
        if consensusXML_missing_file.endswith("consensusXML"):
            ConsensusXMLFile().store(consensusXML_missing_file, missing)

    # takes a (filtered) ConsensusMap (usually the complete) and reconstructs the FeatureMaps
    def consensus_to_feature_maps(self, consensusXML_file, original_featureXML_dir, reconstructed_featureXML_dir):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        to_keep_ids = [item for sublist in [[feature.getUniqueId(
        ) for feature in cf.getFeatureList()] for cf in consensus_map] for item in sublist]
        Helper().reset_directory(reconstructed_featureXML_dir)
        for file in os.listdir(original_featureXML_dir):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(original_featureXML_dir, file), feature_map)
            fm_filterd = FeatureMap(feature_map)
            fm_filterd.clear(False)
            for feature in feature_map:
                if feature.getUniqueId() in to_keep_ids:
                    fm_filterd.push_back(feature)
            FeatureXMLFile().store(os.path.join(reconstructed_featureXML_dir, os.path.basename(
                fm_filterd.getMetaValue("spectra_data")[0].decode())[:-4] + "featureXML"), fm_filterd)

    def merge_feature_maps(self, featureXML_merged_dir, featureXML_dir_a, featureXML_dir_b):
        Helper().reset_directory(featureXML_merged_dir)
        for file_ffm in os.listdir(featureXML_dir_a):
            for file_ffmid in os.listdir(featureXML_dir_b):
                if file_ffm == file_ffmid:
                    fm_ffm = FeatureMap()
                    FeatureXMLFile().load(os.path.join(featureXML_dir_a, file_ffm), fm_ffm)
                    fm_ffmid = FeatureMap()
                    FeatureXMLFile().load(os.path.join(featureXML_dir_b, file_ffm), fm_ffmid)
                    for f in fm_ffmid:
                        fm_ffm.push_back(f)
                    FeatureXMLFile().store(os.path.join(featureXML_merged_dir, file_ffm), fm_ffm)
    
    # create a library for each FeatureXML file with features that are missing in this particular sample
    # in the given ConsensusXMLFile
    def FFMID_libraries_for_missing_features(self, consensusXML_file, library_dir):
        Helper().reset_directory(library_dir)
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        headers = {i:header.filename for i, header in enumerate(consensus_map.getColumnHeaders().values())}
        # write header in all library files here
        all_lib_files = [os.path.join(library_dir, headers[i][:-5] + ".tsv") for i in headers.keys()]
        for file in all_lib_files:
            with open(file, "a", newline="") as tsvfile:
                writer = csv.writer(tsvfile, delimiter="\t")
                writer.writerow(["CompoundName", "SumFormula", "Mass", "Charge", "RetentionTime",
                                "RetentionTimeRange", "IsoDistribution"])
        for i, cf in enumerate(consensus_map):
            lib_files = [file for i, file in enumerate(all_lib_files) if i not in 
                    [f.getMapIndex() for f in cf.getFeatureList()]]
            if lib_files:
                row = ["f_"+str(i), "", cf.getMZ()*cf.getCharge()-(cf.getCharge()*1.007825),
                        cf.getCharge(), cf.getRT(), 0, 0]
                for file in lib_files:
                    with open(file, "a", newline="") as tsvfile:
                        writer = csv.writer(tsvfile, delimiter="\t")
                        writer.writerow(row)



import os
from pyopenms import *

def run(feature_maps, consensusXML_file):
    feature_grouper = FeatureGroupingAlgorithmKD()

    consensus_map = ConsensusMap()
    file_descriptions = consensus_map.getColumnHeaders()

    for i, feature_map in enumerate(feature_maps):
        file_description = file_descriptions.get(i, ColumnHeader())
        file_description.filename = os.path.basename(feature_map.getMetaValue("spectra_data")[0].decode())
        file_description.size = feature_map.size()
        file_descriptions[i] = file_description

    feature_grouper.group(feature_maps, consensus_map)
    consensus_map.setColumnHeaders(file_descriptions)

    consensus_map.setUniqueIds()
    ConsensusXMLFile().store(consensusXML_file, consensus_map)
    print(consensus_map.size())


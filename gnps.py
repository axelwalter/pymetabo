from pyopenms import *
from .helpers import Helper
from pathlib import Path

class GNPSExport:
    def run(self, consensusXML_file, aligned_mzML_dir, mgf_file, quantification_file, metavalue_file, pair_table_file):
        mzML_files = [f for f in Path(aligned_mzML_dir).iterdir() if f.is_file() and str(f).endswith("mzML")]
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        filtered_map = ConsensusMap(consensus_map)
        filtered_map.clear(False)
        for feature in consensus_map:
            if feature.getPeptideIdentifications():
                filtered_map.push_back(feature)

        consensusXML_file = "filtered.consensusXML"
        ConsensusXMLFile().store(consensusXML_file, filtered_map)

        # for FFBM
        GNPSMGFFile().store(String(consensusXML_file), [file.encode() for file in mzML_files], String(mgf_file))
        GNPSQuantificationFile().store(consensus_map, quantification_file)
        GNPSMetaValueFile().store(consensus_map, metavalue_file)

        # for IIMN
        IonIdentityMolecularNetworking().annotateConsensusMap(consensus_map)
        IonIdentityMolecularNetworking().writeSupplementaryPairTable(consensus_map, pair_table_file)
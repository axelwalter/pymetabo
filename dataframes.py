from pyopenms import *

class DataFrames:
    def create_consensus_table(self, consensusXML_file, table_file):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        df = consensus_map.get_df().drop(["sequence"], axis=1)
        for cf in consensus_map:
            if cf.metaValueExists("best ion"):
                df["adduct"] = [cf.getMetaValue("best ion") for cf in consensus_map]
                break
        for cf in consensus_map:
            if cf.metaValueExists("label"):
                df["name"] = [cf.getMetaValue("label") for cf in consensus_map]
                break
        df.to_csv(table_file, sep="\t")
        return df
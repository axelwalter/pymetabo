from matplotlib.pyplot import table
from pyopenms import *
import pandas as pd

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
        if table_file.endswith("tsv"):
            df.to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.to_feather(table_file)
        return df
    
    def FFMID_chroms_to_df(self, featureXML_file, table_file):
        fm = FeatureMap()
        FeatureXMLFile().load(featureXML_file, fm)
        chroms = {}
        for f in fm:
            for i, sub in enumerate(f.getSubordinates()):
                name = f.getMetaValue('label') + "_" + str(i+1)
                chroms[name + "_int"] = [int(y[1]) for y in sub.getConvexHulls()[0].getHullPoints()]
                chroms[name + "_RT"] = [x[0] for x in sub.getConvexHulls()[0].getHullPoints()]
        df = pd.DataFrame({ key:pd.Series(value) for key, value in chroms.items() })
        if table_file.endswith("tsv"):
            df.to_csv(featureXML_file[:-10]+"tsv", sep="\t")
        elif table_file.endswith("ftr"):
            df.to_feather(table_file)
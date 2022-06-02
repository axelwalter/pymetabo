from matplotlib.pyplot import table
from pyopenms import *
import pandas as pd
import os

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
            df.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)
        return df
    
    def FFMID_chroms_to_df(self, featureXML_file, table_file, time_unit = "seconds"):
        time_factor = 1
        if time_unit == "minutes":
            time_factor = 60
        fm = FeatureMap()
        FeatureXMLFile().load(featureXML_file, fm)
        chroms = {}
        for f in fm:
            for i, sub in enumerate(f.getSubordinates()):
                name = f.getMetaValue('label') + "_" + str(i+1)
                chroms[name + "_int"] = [int(y[1]) for y in sub.getConvexHulls()[0].getHullPoints()]
                chroms[name + "_RT"] = [x[0]/time_factor for x in sub.getConvexHulls()[0].getHullPoints()]
        df = pd.DataFrame({ key:pd.Series(value) for key, value in chroms.items() })
        if table_file.endswith("tsv"):
            df.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

    def FFMID_auc_to_df(self, featureXML_file, table_file):
        fm = FeatureMap()
        FeatureXMLFile().load(featureXML_file, fm)
        aucs = {}
        for f in fm:
            aucs[f.getMetaValue('label')] = [int(f.getIntensity())]
        df = pd.DataFrame({ key:pd.Series(value) for key, value in aucs.items() })
        if table_file.endswith("tsv"):
            df.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

    def FFMID_auc_combined_to_df(self, df_auc_file, table_file):
        if df_auc_file.endswith("tsv"):
            df = pd.read_csv(df_auc_file, sep="\t")
        elif df_auc_file.endswith("ftr"):
            df = pd.read_feather(df_auc_file).drop(columns=["index"])
        aucs_condensed = {}
        for a in set([c.split("#")[0] for c in df.columns]):
            aucs_condensed[a] = 0
            for b in [b for b in df.columns if ((a+"#" in b and b.startswith(a)) or a == b)]:
                aucs_condensed[a] += df[b][0]
        df_combined = pd.DataFrame({ key:pd.Series(value) for key, value in aucs_condensed.items() })
        if table_file.endswith("tsv"):
            df_combined.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df_combined.reset_index().to_feather(table_file)

    def get_auc_summary(self, df_files, table_file):
        # get a list of auc dataframe file paths (df_files), combine them into a summary (consensus) df
        dfs = []
        indeces = []
        for file in df_files:
            if file.endswith("tsv"):
                df = pd.read_csv(file, sep="\t")
            elif file.endswith("ftr"):
                df = pd.read_feather(file).drop(columns=["index"])
            print(df.columns)
            dfs.append(df)
            indeces.append(os.path.basename(file)[:-4])
        summary = pd.concat(dfs)
        summary = summary.set_index(pd.Series(indeces))
        if table_file.endswith("tsv"):
            summary.reset_index().to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            summary.reset_index().to_feather(table_file)

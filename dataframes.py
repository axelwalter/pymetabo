from matplotlib.pyplot import table
from pyopenms import *
import pandas as pd
import os
from pathlib import Path

class DataFrames:

    def what():
        return 1
    def create_consensus_table(self, consensusXML_file, table_file, sirius_ms_dir=""):
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
        if "adduct" in df.columns:
            df.index = [f"{round(mz, 4)}@{int(rt)}@{adduct}" for mz, rt, adduct in zip(df["mz"].tolist(), df["RT"].tolist(), df["adduct"].tolist())]
        else:    
            df.index = [f"{round(mz, 4)}@{int(rt)}" for mz, rt in zip(df["mz"].tolist(), df["RT"].tolist())]
        not_sample = [c for c in df.columns if c not in ["mz", "RT", "charge", "adduct", "name", "quality"]]
        df[not_sample] = df[not_sample].applymap(lambda x: int(round(x, 0)) if isinstance(x, (int, float)) else x)

        # annotate original feature Ids which are in the Sirius .ms files
        if sirius_ms_dir:
            ms_files = [Path(sirius_ms_dir, file) for file in os.listdir(sirius_ms_dir)]
            map = {Path(value.filename).stem: key for key, value in consensus_map.getColumnHeaders().items()}
            for file in ms_files:
                if file.exists():
                    key = map[file.stem]
                    id_list = []
                    content = file.read_text()
                    for cf in consensus_map:
                        # get a map with map index and feature id for each consensus feature -> get the feature id key exists
                        f_map = {fh.getMapIndex(): fh.getUniqueId() for fh in cf.getFeatureList()}
                        if key in f_map.keys():
                            f_id = str(f_map[key])
                        else:
                            f_id = ""
                        if f_id and f_id in content:
                            id_list.append(f_id)
                        else:
                            id_list.append("")
                    df[file.stem+"_SiriusID"] = id_list


        if table_file.endswith("tsv"):
            df.to_csv(table_file, sep="\t")
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
        empty = []
        for file in df_files:
            if file.endswith("tsv"):
                df = pd.read_csv(file, sep="\t")
            elif file.endswith("ftr"):
                df = pd.read_feather(file)
                if "index" in df.columns:
                    df.index = df["index"]
                    df = df.drop(columns=["index"])
            sample_name = os.path.basename(file)[:-4].split("AUC")[0]
            if df.empty:
                empty.append(sample_name)
            else:
                dfs.append(df)
                indeces.append(sample_name)
        df = pd.concat(dfs)
        df = df.set_index(pd.Series(indeces))
        df = df.transpose()
        df = df.fillna(0)
        for sample in empty:
            df[sample] = 0
        df = df.applymap(lambda x: int(round(x, 0)) if isinstance(x, (int, float)) else x)
        df.sort_index(axis=1, inplace=True)
        if table_file.endswith("tsv"):
            df.to_csv(table_file, sep="\t")
        elif table_file.endswith("ftr"):
            df.reset_index().to_feather(table_file)

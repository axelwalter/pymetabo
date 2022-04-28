import os
from pyopenms import *
import pandas as pd

def create_table_FFM(consensusXML_file, table_file):
    consensus_map = ConsensusMap()
    ConsensusXMLFile().load(consensusXML_file, consensus_map)
    consensus_map.get_df().drop(["sequence"], axis=1).to_csv(table_file, sep="\t")

def load_FFMID_library(consensusXML_file, library_file):
    consensus_map = ConsensusMap()
    ConsensusXMLFile().load(consensusXML_file, consensus_map)
    #Import the consensus tsv table and keep only the columns: RT, mz and charge
    library = consensus_map.get_df()[['RT','mz', "charge"]]
    #convert the mz and RT columns to floats and charge to integer for calculations
    library["charge"] = pd.to_numeric(library["charge"], downcast="integer")
    library["mz"] = pd.to_numeric(library["mz"], downcast="float")
    library["RT"] = pd.to_numeric(library["RT"], downcast="float")
    library= library.rename(columns={"RT": "RetentionTime", "charge":"Charge"})
    #Add a columns named "Mass" and calculate the neutral masses from the charge and mz:
    library["Mass"] = 0.0
    for i in library.index:
        library.at[i, "Mass"] = (library.loc[i, "mz"] * library.loc[i, "Charge"]) - (library.loc[i, "Charge"] * 1.007825)
    #drop the mz column
    library= library.drop(columns= "mz")
    library["Charge"] = [[c] for c in library["Charge"]]
    library["RetentionTime"] = [[rt] for rt in library["RetentionTime"]]
    #add the rest of the columns required for the MetaboliteIdentificationTable and fill with zeros or blanks, except the "Compound Name"
    #which, since they are all unknown, can be filled with f_#
    library['CompoundName'] = [i for i in range(0,len(library))]
    library['CompoundName'] = "f_" + library['CompoundName'].astype(str)
    library["SumFormula"] = ""
    library["RetentionTimeRange"]= [[0.0] for _ in range(len(library.index))]
    library["IsoDistribution"]= [[0.0] for _ in range(len(library.index))]
    library = library[["CompoundName","SumFormula", "Mass","Charge","RetentionTime","RetentionTimeRange", "IsoDistribution"]]
    library.to_csv(library_file, sep="\t")
    library_reformatted = []
    for _, row in library.iterrows():
        library_reformatted.append(FeatureFinderMetaboIdentCompound(row["CompoundName"], row["SumFormula"], row["Mass"], row["Charge"], row["RetentionTime"], row["RetentionTimeRange"], row["IsoDistribution"]))
    return library_reformatted
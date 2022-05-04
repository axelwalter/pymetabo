import os
import shutil
import csv
import pandas as pd
from pyopenms import *

# TODO FeatureFinderMetaboIdent.create_template_library 
# TODO instead of directory (eg. with FeatureXML files, accept also list of filenames)

class Helper:
    def reset_directory(path: str):
        try:
            shutil.rmtree(path)
            os.mkdir(path)
        except OSError:
            os.mkdir(path)

    def load_feature_maps(path: str):
        feature_maps = []
        for file in os.listdir(path):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(path, file), feature_map)
            feature_maps.append(feature_map)
        return feature_maps

class FeatureFinderMetabo:
    def run(mzML, featureXML, params = {}, q_threshold = 0):
        if os.path.isdir(mzML):
            mzML_files = [os.path.join(mzML, file) for file in os.listdir(mzML)]
        else:
            mzML_files = [mzML]
        if not featureXML.endswith(".featureXML"):
            Helper.reset_directory(featureXML)
        for mzML_file in mzML_files:
            exp = MSExperiment()
            MzMLFile().load(mzML_file, exp)
            exp.sortSpectra(True)

            mtd = MassTraceDetection()
            mtd_par = mtd.getDefaults()
            for key, value in params.items():
                if key.encode() in mtd_par.keys():
                    mtd_par.setValue(key, value) 
            mtd.setParameters(mtd_par)
            mass_traces = []
            mtd.run(exp, mass_traces, 0) # input MSExperiment, empty list for detected mass traces, max_size (if not 0, sets the maximum number of mass traces)

            epd = ElutionPeakDetection()
            epd_par = epd.getDefaults()
            for key, value in params.items():
                if key.encode() in epd_par.keys():
                    epd_par.setValue(key, value) 
            epd.setParameters(epd_par)
            elution_peaks = []
            epd.detectPeaks(mass_traces, elution_peaks) # list with all detected mass traces, list of mass traces that represent an elution peak

            ffm = FeatureFindingMetabo()
            ffm_par = ffm.getDefaults() 
            for key, value in params.items():
                if key.encode() in ffm_par.keys():
                    ffm_par.setValue(key, value)
            ffm.setParameters(ffm_par)
            feature_map = FeatureMap()
            feat_chrom = []
            ffm.run(elution_peaks, feature_map, feat_chrom) # elution peaks, empty FeatureMap, empty list for feature chromatograms

            feature_map.setUniqueIds()
            feature_map.setPrimaryMSRunPath([mzML_file.encode()])

            feature_map_filtered = FeatureMap(feature_map)
            feature_map_filtered.clear(False)

            if q_threshold:
                for f in feature_map:
                        if f.getOverallQuality() > q_threshold: # 0.0005 was good
                            feature_map_filtered.push_back(f)
                print('Features before quality filter: ' + str(feature_map.size()))
                print('Features after quality filter: ' + str(feature_map_filtered.size()))
                feature_map = feature_map_filtered
            if os.path.isdir(featureXML):
                FeatureXMLFile().store(os.path.join(featureXML, os.path.basename(mzML_file)[:-4] + "featureXML"), feature_map)
            else:
                FeatureXMLFile().store(featureXML, feature_map)

class MapAligner:
    def run(input_files, aligned_dir, trafo_dir, params = {}):
        aligner = MapAlignmentAlgorithmPoseClustering()
        aligner_par= aligner.getDefaults()
        for key, value in params.items():
            if key.encode() in aligner_par.keys():
                aligner_par.setValue(key, value)
        aligner.setParameters(aligner_par)
        inputs = os.listdir(input_files)
        if inputs and inputs[0].endswith("featureXML"):
            Helper.reset_directory(aligned_dir)
            Helper.reset_directory(trafo_dir)
            feature_maps = Helper.load_feature_maps(input_files)
            transformations = {} # store TransformationDescriptions for MapAlignmentTransformer of MSExperiments during Requantification
            ref_index = feature_maps.index(sorted(feature_maps, key=lambda x:x.size())[-1])
            aligner.setReference(feature_maps[ref_index])
            print("Map Alignment reference map: ", feature_maps[ref_index].getMetaValue("spectra_data")[0].decode())

            for feature_map in feature_maps[:ref_index] + feature_maps[ref_index+1:]:
                trafo = TransformationDescription()
                aligner.align(feature_map, trafo) # store information on aligmentment in TransformationDescription, RTs in FeatureMap not modified at this point
                transformer = MapAlignmentTransformer() 
                transformer.transformRetentionTimes(feature_map, trafo, True) # FeatureMap, TransformationDescription, bool: keep original RTs as meta value
                transformations[feature_map.getMetaValue("spectra_data")[0].decode()] = trafo
                TransformationXMLFile().store(os.path.join(trafo_dir, os.path.basename(feature_map.getMetaValue("spectra_data")[0].decode())[:-4] + "trafoXML"), trafo)

            for feature_map in feature_maps:
                print(feature_map.size())
                FeatureXMLFile().store(os.path.join(aligned_dir, os.path.basename(feature_map.getMetaValue("spectra_data")[0].decode())[:-4] + "featureXML"), feature_map)

        elif inputs and inputs[0].endswith("mzML"):
            Helper.reset_directory(aligned_dir)
            for file in os.listdir(input_files):
                exp = MSExperiment()
                MzMLFile().load(os.path.join(input_files, file), exp)
                exp.sortSpectra(True)
                if file[:-4] + "trafoXML" not in os.listdir(trafo_dir):
                    MzMLFile().store(os.path.join(aligned_dir, file), exp)
                    continue
                transformer = MapAlignmentTransformer()
                trafo_description = TransformationDescription()
                TransformationXMLFile().load(os.path.join(trafo_dir, file[:-4] + "trafoXML"), trafo_description, False)
                transformer.transformRetentionTimes(exp, trafo_description, True)
                MzMLFile().store(os.path.join(aligned_dir, file), exp)

class FeatureLinker:
    def run(featureXML_dir, consensusXML_file, params = {}):
        feature_maps = Helper.load_feature_maps(featureXML_dir)
        feature_grouper = FeatureGroupingAlgorithmKD()

        feature_grouper_params = feature_grouper.getDefaults() 
        for key, value in params.items():
            if key.encode() in feature_grouper_params.keys():
                feature_grouper_params.setValue(key, value)
        feature_grouper.setParameters(feature_grouper_params)

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
        Helper.reset_directory(os.path.dirname(consensusXML_file))
        ConsensusXMLFile().store(consensusXML_file, consensus_map)
        print(f"ConsensusMap size: {consensus_map.size()}")

class DataFrames:
    def create_consensus_table(consensusXML_file, table_file):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        df = consensus_map.get_df().drop(["sequence"], axis=1)
        df.to_csv(table_file, sep="\t")
        return df

class FeatureFinderMetaboIdent:
    def load_library(input_file, library_file = ""):
        # input file can be a consensusXML or tsv file
        if input_file.endswith("consensusXML"):
            consensus_map = ConsensusMap()
            ConsensusXMLFile().load(input_file, consensus_map)
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
            if library_file:
                library.to_csv(library_file, sep="\t")
            metabo_table = []
            for _, row in library.iterrows():
                metabo_table.append(FeatureFinderMetaboIdentCompound(row["CompoundName"], row["SumFormula"], row["Mass"], row["Charge"], row["RetentionTime"], row["RetentionTimeRange"], row["IsoDistribution"]))
            return metabo_table
        elif input_file.endswith("tsv"):
            metabo_table = []
            with open(input_file, 'r') as tsv_file:
                tsv_reader = csv.reader(tsv_file, delimiter="\t")
                next(tsv_reader) # skip header
                for row in tsv_reader:
                    metabo_table.append(FeatureFinderMetaboIdentCompound(
                        row[0], # name
                        row[1], # sum formula
                        float(row[2]), # mass
                        [int(charge) for charge in row[3].split(',')], # charges
                        [float(rt) for rt in row[4].split(',')], # RTs
                        [float(rt_range) for rt_range in row[5].split(',')], # RT ranges
                        [float(iso_distrib) for iso_distrib in row[6].split(',')] # isotope distributions
                    ))
            return metabo_table

    def create_template_library(file_path):
        if file_path.endswith("tsv"):
            pass

    def run(mzML, featureXML, library, params = {}):
        if os.path.isdir(mzML):
            mzML_files = [os.path.join(mzML, file) for file in os.listdir(mzML)]
        else:
            mzML_files = [mzML]
        if not featureXML.endswith(".featureXML"): # -> it is a directory
            Helper.reset_directory(featureXML)
        for mzML_file in mzML_files:
            exp = MSExperiment()
            MzMLFile().load(mzML_file, exp)
            ffmid = FeatureFinderAlgorithmMetaboIdent()
            ffmid.setMSData(exp)
            feature_map = FeatureMap()
            ffmid_params = ffmid.getDefaults()
            for key, value in params.items():
                if key.encode() in ffmid_params.keys():
                    ffmid_params.setValue(key, value)
            ffmid.setParameters(ffmid_params)
            # run the FeatureFinderMetaboIdent with the metabo_table and aligned mzML file path
            ffmid.run(library, feature_map, mzML_file)
            print(feature_map.size())
            feature_map.setUnassignedPeptideIdentifications([])
            feature_map.setProteinIdentifications([])
            if os.path.isdir(featureXML):
                FeatureXMLFile().store(os.path.join(featureXML, os.path.basename(mzML_file)[:-4] + "featureXML"), feature_map)
            else:
                FeatureXMLFile().store(featureXML, feature_map)

class FeatureMapHelper:
    def split_consensus_map(consensusXML_file, consensusXML_complete_file="", consensusXML_missing_file=""):
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
    def consensus_to_feature_maps(consensusXML_file, original_featureXML_dir, reconstructed_featureXML_dir):
        consensus_map = ConsensusMap()
        ConsensusXMLFile().load(consensusXML_file, consensus_map)
        to_keep_ids = [item for sublist in [[feature.getUniqueId() for feature in cf.getFeatureList()] for cf in consensus_map] for item in sublist]
        Helper.reset_directory(reconstructed_featureXML_dir)
        for file in os.listdir(original_featureXML_dir):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(original_featureXML_dir, file), feature_map)
            fm_filterd = FeatureMap(feature_map)
            fm_filterd.clear(False)
            for feature in feature_map:
                if feature.getUniqueId() in to_keep_ids:
                    fm_filterd.push_back(feature)
            FeatureXMLFile().store(os.path.join(reconstructed_featureXML_dir, os.path.basename(fm_filterd.getMetaValue("spectra_data")[0].decode())[:-4] + "featureXML"), fm_filterd)

    def merge_feature_maps(featureXML_merged_dir, featureXML_dir_a, featureXML_dir_b):
        Helper.reset_directory(featureXML_merged_dir)
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

class MetaboliteAdductDecharger:
    def run(fm_dir, fm_decharged_dir, params = {}):
        Helper.reset_directory(fm_decharged_dir)
        for file in os.listdir(fm_dir):
            feature_map = FeatureMap()
            FeatureXMLFile().load(os.path.join(fm_dir, file), feature_map)
            mfd = MetaboliteFeatureDeconvolution()
            mdf_par = mfd.getDefaults()
            for key, value in params.items():
                if key.encode() in mdf_par.keys():
                    mdf_par.setValue(key, value)
            mfd.setParameters(mdf_par)

            feature_map_decharged = FeatureMap()
            mfd.compute(feature_map, feature_map_decharged, ConsensusMap(), ConsensusMap())
            FeatureXMLFile().store(os.path.join(fm_decharged_dir, file), feature_map_decharged)
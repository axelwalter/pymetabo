import os
from pyopenms import *

def align_feature_maps(feature_maps, aligned_dir, trafo_dir):
    transformations = {} # store TransformationDescriptions for MapAlignmentTransformer of MSExperiments during Requantification
    aligner = MapAlignmentAlgorithmPoseClustering()
    aligner_par= aligner.getDefaults()
    aligner_par.setValue("max_num_peaks_considered", -1)
    aligner_par.setValue("superimposer:mz_pair_max_distance", 0.05)
    aligner_par.setValue("pairfinder:distance_MZ:max_difference", 10.0)
    aligner_par.setValue("pairfinder:distance_MZ:unit", "ppm")
    aligner.setParameters(aligner_par)

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
    
    return transformations

def align_mzML(mzML_dir, aligned_dir, trafo_dir):
    for file in os.listdir(mzML_dir):
        exp = MSExperiment()
        MzMLFile().load(os.path.join(mzML_dir, file), exp)
        exp.sortSpectra(True)
        if file[:-4] + "trafoXML" not in os.listdir(trafo_dir):
            MzMLFile().store(os.path.join(aligned_dir, file), exp)
            continue
        transformer = MapAlignmentTransformer()
        trafo_description = TransformationDescription()
        TransformationXMLFile().load(os.path.join(trafo_dir, file[:-4] + "trafoXML"), trafo_description, False)
        transformer.transformRetentionTimes(exp, trafo_description, True)
        MzMLFile().store(os.path.join(aligned_dir, file), exp)
from pyopenms import *

def run(mzML_file, featureXML_file):
    exp = MSExperiment()
    MzMLFile().load(mzML_file, exp)
    exp.sortSpectra(True)

    mtd = MassTraceDetection()
    mtd_par = mtd.getDefaults()
    mtd_par.setValue("mass_error_ppm", 10.0) 
    mtd_par.setValue("noise_threshold_int", 10000.0)
    mtd.setParameters(mtd_par)
    mass_traces = []
    mtd.run(exp, mass_traces, 0) # input MSExperiment, empty list for detected mass traces, max_size (if not 0, sets the maximum number of mass traces)

    epd = ElutionPeakDetection()
    epd_par = epd.getDefaults()
    epd_par.setValue("width_filtering", "fixed")
    epd.setParameters(epd_par)
    elution_peaks = []
    epd.detectPeaks(mass_traces, elution_peaks) # list with all detected mass traces, list of mass traces that represent an elution peak

    ffm = FeatureFindingMetabo()
    ffm_par = ffm.getDefaults() 
    ffm_par.setValue("isotope_filtering_model", "none")
    ffm_par.setValue("remove_single_traces", "true")
    ffm_par.setValue("mz_scoring_by_elements", "false")
    ffm_par.setValue("report_convex_hulls", "true")
    ffm.setParameters(ffm_par)
    feature_map = FeatureMap()
    feat_chrom = []
    ffm.run(elution_peaks, feature_map, feat_chrom) # elution peaks, empty FeatureMap, empty list for feature chromatograms

    feature_map.setUniqueIds()
    feature_map.setPrimaryMSRunPath([mzML_file.encode()])

    feature_map_filtered = FeatureMap(feature_map)
    feature_map_filtered.clear(False)


    for f in feature_map:
        if f.getOverallQuality() > 0.001: # 0.0005 was good
            feature_map_filtered.push_back(f)


    print('Features before quality filter: ' + str(feature_map.size()))
    print('Features after quality filter: ' + str(feature_map_filtered.size()))

    feature_map = feature_map_filtered

    FeatureXMLFile().store(featureXML_file, feature_map)
import os
import shutil
from pyopenms import *
from .helpers import Helper
from pathlib import Path

class Sirius:
    def run(self, mzML_dir, featureXML_dir, sirius_dir, sirius_executable_path, only_export_ms_file, params={}):
        Helper().reset_directory(sirius_dir)

        feature_files = sorted(os.listdir(featureXML_dir))
        mzml_files = sorted(os.listdir(mzML_dir))
        no_convex_hulls_dir = Helper().reset_directory(os.path.join(sirius_dir, "no_convex_hulls"))
        if not only_export_ms_file:
            formula_dir = Helper().reset_directory(os.path.join(sirius_dir, "formulas"))
            structure_dir = Helper().reset_directory(os.path.join(sirius_dir, "structures"))
        ms_dir = Helper().reset_directory(os.path.join(sirius_dir, "sirius_files"))

        for mzML_file in mzml_files:
            exp = MSExperiment()
            MzMLFile().load(os.path.join(mzML_dir, mzML_file), exp)
            exp.sortSpectra(True)
    
            for feature_file in feature_files:
                if mzML_file[:-5] == feature_file[:-11]:
                    fm = FeatureMap()
                    FeatureXMLFile().load(os.path.join(featureXML_dir, feature_file), fm)
                    # remove subordinates
                    fm_no_sub = FeatureMap(fm)
                    fm_no_sub.clear(False)
                    for f in fm:
                        f.setConvexHulls([])
                        f.setSubordinates([])
                        fm_no_sub.push_back(f)
                    sirius_feature_path = os.path.join(no_convex_hulls_dir, feature_file)
                    FeatureXMLFile().store(sirius_feature_path, fm_no_sub)

                    sirius_algo = SiriusAdapterAlgorithm()
                    sirius_algo_par = sirius_algo.getDefaults()
                    for key, value in params.items():
                        if key.encode() in sirius_algo_par.keys():
                            sirius_algo_par.setValue(key, value)
                    sirius_algo.setParameters(sirius_algo_par)
                        
                    fm_info = FeatureMapping_FeatureMappingInfo()
                    mapping = FeatureMapping_FeatureToMs2Indices() 
                    sirius_algo.preprocessingSirius(sirius_feature_path,
                                                    exp,
                                                    fm_info,
                                                    mapping)
                    sirius_algo.logFeatureSpectraNumber(sirius_feature_path,
                                                        mapping,
                                                        exp)
                    msfile = SiriusMSFile()
                    sirius_tmp = SiriusTemporaryFileSystemObjects(3) # debug level
                    msfile.store(exp,
                                String(sirius_tmp.getTmpMsFile()),
                                mapping, 
                                sirius_algo.isFeatureOnly(),
                                sirius_algo.getIsotopePatternIterations(), 
                                sirius_algo.isNoMasstraceInfoIsotopePattern(), 
                                []) # empty compound info
                    
                    # save only files with MS2 data (others will be empty)
                    if os.path.isfile(sirius_tmp.getTmpMsFile()) and os.path.getsize(sirius_tmp.getTmpMsFile()) > 0:
                        shutil.copy(sirius_tmp.getTmpMsFile(), os.path.join(ms_dir, mzML_file[:-4]+"ms"))

                    if only_export_ms_file:
                        continue

                    out_csifingerid = os.path.join(structure_dir, mzML_file[:-5] +".mzTab")
                    subdirs = sirius_algo.callSiriusQProcess(String(sirius_tmp.getTmpMsFile()),
                                                            String(sirius_tmp.getTmpOutDir()),
                                                            String(sirius_executable_path),
                                                            String(out_csifingerid),
                                                            False)
                    sirius_result = MzTab()
                    SiriusMzTabWriter.read(subdirs,
                                        mzML_file,
                                        sirius_algo.getNumberOfSiriusCandidates(),
                                        sirius_result)

                    sirius_file = os.path.join(formula_dir, mzML_file[:-5] +".mzTab")
                    MzTabFile().store(sirius_file, sirius_result)

        # finally delete the temporary FeatureXML file without convex hulls
        shutil.rmtree(no_convex_hulls_dir)


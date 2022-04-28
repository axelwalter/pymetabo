import os
from pyopenms import *

def run(ffmid_dir, mzML_dir, library):
    for file in os.listdir(mzML_dir):
        exp = MSExperiment()
        MzMLFile().load(os.path.join(mzML_dir, file), exp)
        ffmid = FeatureFinderAlgorithmMetaboIdent()
        ffmid.setMSData(exp)
        fm = FeatureMap()
        params = ffmid.getDefaults()
        params[b'extract:mz_window'] = 10.0 # 5 ppm
        # params[b'extract:rt_window'] = 20.0 # 20 seconds
        params[b'detect:peak_width'] = 60.0 # 3 seconds
        # params[b'signal_to_noise'] = 3.0 # 0.8
        params[b'model:type'] = 'symmetric'
        # params[b'check:width'] = 50
        ffmid.setParameters(params)
        # run the FeatureFinderMetaboIdent with the metabo_table and aligned mzML file path
        ffmid.run(library, fm, file)
        print(fm.size())
        fm.setUnassignedPeptideIdentifications([])
        fm.setProteinIdentifications([])
        FeatureXMLFile().store(os.path.join(ffmid_dir, file[:-4] + "featureXML"), fm)
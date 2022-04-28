import os
from pyopenms import *

def run(mzML_dir, fm_dir, fm_mapped_dir):
    use_centroid_rt= False
    use_centroid_mz= True

    mapper = IDMapper()

    for mzML_file in os.listdir(mzML_dir):
        exp = MSExperiment()
        MzMLFile().load(os.path.join(mzML_dir, mzML_file), exp)
        for feature_file in os.listdir(fm_dir):
            fm = FeatureMap()
            FeatureXMLFile().load(os.path.join(fm_dir, feature_file), fm)
            if feature_file[:-10] == mzML_file[:-4]:
                peptide_ids = []
                protein_ids = []
                mapper.annotate(fm, peptide_ids, protein_ids, use_centroid_rt, use_centroid_mz, exp)
                FeatureXMLFile().store(os.path.join(fm_mapped_dir, feature_file), fm)
import os
from pyopenms import *

def run(fm_dir, fm_decharged_dir):
    for file in os.listdir(fm_dir):
        feature_map = FeatureMap()
        FeatureXMLFile().load(os.path.join(fm_dir, file), feature_map)
        mfd = MetaboliteFeatureDeconvolution()
        mdf_par = mfd.getDefaults()
        mdf_par.setValue("potential_adducts", [b"H:+:0.5",b"Na:+:0.3", b"H-1O-1:+:0.2"])
        mdf_par.setValue("charge_min", 1, "Minimal possible charge")
        mdf_par.setValue("charge_max", 1, "Maximal possible charge")
        mdf_par.setValue("charge_span_max", 1)
        mdf_par.setValue("max_neutrals", 1)
        mfd.setParameters(mdf_par)

        feature_map_decharged = FeatureMap()
        mfd.compute(feature_map, feature_map_decharged, ConsensusMap(), ConsensusMap())
        FeatureXMLFile().store(os.path.join(fm_decharged_dir, file), feature_map_decharged)
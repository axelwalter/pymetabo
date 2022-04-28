import os
import shutil
from pyopenms import FeatureMap, FeatureXMLFile

def reset_directory(path: str):
    try:
        shutil.rmtree(path)
        os.mkdir(path)
    except OSError:
        os.mkdir(path)
    return path

def load_feature_maps(path: str):
    feature_maps = []
    for file in os.listdir(path):
        feature_map = FeatureMap()
        FeatureXMLFile().load(os.path.join(path, file), feature_map)
        feature_maps.append(feature_map)
    return feature_maps
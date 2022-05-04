# Metabolomics workflow with pyOpenMS


- TODO FeatureFinderMetaboIdent.create_template_library
- TODO instead of directory (eg. with FeatureXML files, accept also list of filenames)
- TODO single files for MetaboliteAdductDecharger and IDMapper

This module contains pyOpenMS tools for metabolomics wrapped as classes that can be used to create clean workflows.
File input and output is handled by passing either single files or directories with multiple files. Currently passing pyOpenMS
objects directly without writing them to file is not yet supported.

The module contains example files to demonstrate the use case.

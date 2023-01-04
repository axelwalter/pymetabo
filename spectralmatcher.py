from pyopenms import *

class SpectralMatcher:
    def run(self, database, mgf_file, output_mztab):
        # Load the MGF file to an MSExperiment format:
        exp = MSExperiment()
        MascotGenericFile().load(mgf_file, exp)

        # Perform spectral matching with a library in MGF format that is located under "resources":
        speclib = MSExperiment()
        MascotGenericFile().load(database, speclib)
        mztab= MzTab()
        out_merged= ""
        MSMS_match= MetaboliteSpectralMatching()
        MSMS_match_par = MSMS_match.getDefaults()
        MSMS_match_par.setValue('merge_spectra', 'false')
        MSMS_match.setParameters(MSMS_match_par)
        MSMS_match.run(exp, speclib, mztab,  String(out_merged))
        MzTabFile().store(output_mztab, mztab)

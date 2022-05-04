import pandas as pd
import streamlit as st
import os
import plotly.express as px
from modules.helpers import reset_directory, load_feature_maps
import modules.FeatureFinderMetabo as FFM
import modules.MapAligner as MapAligner
import modules.FeatureLinker as FeatureLinker
import modules.DataFrames as DataFrames
import modules.FeatureFinderMetaboIdent as FFMID
import modules.AdductDecharger as AdductDecharger
import modules.IDMapper as IDMapper
import modules.FeatureMapFunctions as FeatureMapFunctions
st.set_page_config(layout="wide")

st.header("pyMetabolomics")

st.markdown("**Specify input and result folders:**")
mzML_dir = st.text_input("Path to input mzML file folder.",
                         "/home/axel/Nextcloud/workspace/MetabolomicsWorkflowMayer/mzML")
result_dir = st.text_input("Path to result folder.", "results")

st.markdown("**Feature Detection (FeatureFinderMetabo):**")

st.markdown("**Feature Map Alignment:**")

st.markdown("**Feature Linking:**")

st.markdown("**Requantification (FeatureFinderMetaboIdent):**")
ffmid = st.checkbox("requantify missing values", value=True)
if ffmid:
    st.date_input("pick a date")

if st.button("Run Workflow"):
    reset_directory(result_dir)
    prog = st.progress(0)
    # FeatureFinderMetabo
    ffm_dir = reset_directory(os.path.join(result_dir, "FFM"))
    for i, mzML_file in enumerate(os.listdir(mzML_dir)):
        if not mzML_file.endswith("mzML"):
            continue
        FFM.run(os.path.join(mzML_dir, mzML_file), os.path.join(
            ffm_dir, mzML_file[:-4] + "featureXML"))
        prog.progress((i+1)/len(os.listdir(mzML_dir)))
    st.success('Feature detection by FFM done!')

    # MapAligner for Feature Maps
    ffm_aligned_dir = reset_directory(os.path.join(result_dir, "FFM_aligned"))
    trafo_dir = reset_directory(os.path.join(result_dir, "Trafo"))
    MapAligner.align_feature_maps(load_feature_maps(
        ffm_dir), ffm_aligned_dir, trafo_dir)
    st.success('Map alignment of FFM feature maps done!')

    # FeatureLinker
    ffm_consensus_dir = reset_directory(
        os.path.join(result_dir, "FFM_consensus"))
    FeatureLinker.run(load_feature_maps(ffm_aligned_dir),
                      os.path.join(ffm_consensus_dir, "FFM.consensusXML"))
    DataFrames.create_table_FFM(os.path.join(
        ffm_consensus_dir, "FFM.consensusXML"), os.path.join(ffm_consensus_dir, "FFM.tsv"))
    st.success('Feature Linking of FFM feature maps done!')

    ffm_table = pd.read_csv(os.path.join(
        ffm_consensus_dir, "FFM.tsv"), sep="\t")
    st.markdown(
        "**Intermediate result:** Consensus features from FeatureFinderMetabo.")
    st.write(ffm_table)
    col1, col2 = st.columns(2)
    df = ffm_table.drop(columns=["mz", "RT", "charge", "id", "quality"])
    na_values = {}
    for column in df.columns:
        na_values[column] = df[column].value_counts()[0]
    with col1:
        st.bar_chart(pd.Series(na_values, name="missing values"))
    with col2:
        st.metric("highest intensity", max(df.max(axis=1)))
        st.metric("median quality", ffm_table["quality"].median())
        st.metric("mean quality", ffm_table["quality"].mean())
    fig = px.scatter(ffm_table, x="RT", y="mz", hover_name="quality")
    st.plotly_chart(fig)

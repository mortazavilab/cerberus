import streamlit as st
import cerberus
import matplotlib.pyplot as plt
import tempfile
from utils import *

def load_data(fname, og_fname):
    try:
        ca = load_cerberus_from_path(fname)
        st.session_state.ca = ca
        st.session_state.ca_path = og_fname

        # for each triplet set, get the valid gnames,
        # marker colors, and marker sizes
        st.session_state.gnames = {}
        st.session_state.num_cols = {}
        st.session_state.all_cols = {}

        for triplet_set in ca.triplets.source.unique():
            temp = ca.triplets.loc[ca.triplets.source==triplet_set]
            st.session_state.gnames[triplet_set] = sorted(
                temp.gname
                .dropna()
                .astype(str)
                .unique()
                .tolist()
            )

            st.session_state.num_cols[triplet_set] = sorted(
                temp
                .select_dtypes(include="number")
                .dropna(axis=1, how="all")
                .columns
                .tolist()
            )

            st.session_state.all_cols[triplet_set] = sorted(
                temp
                .dropna(axis=1, how='all')
                .columns
                .tolist()
            )

    except Exception as e:
        st.error(f"Failed to load data: {e}")
        st.stop()

@st.cache_resource(show_spinner=False)
def load_cerberus_from_path(path: str):
    return cerberus.read(path)

def render_landing_tab():
    """
    Upload data
    """

    uploaded_file = st.file_uploader(
        "Upload CerberusAnnotation",
        type=["h5"]
    )

    if uploaded_file is not None:
        if st.session_state.get("ca_path") != uploaded_file.name:
            print('UWUWUUWUWUUWUWUWU')
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5") as tmp:
                tmp.write(uploaded_file.getbuffer())
                tmp_path = tmp.name

            load_data(tmp_path, uploaded_file.name)

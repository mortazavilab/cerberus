import streamlit as st
import cerberus
import matplotlib.pyplot as plt
import tempfile
from utils import *

def load_data(fname):
    try:
        ca = load_cerberus_from_path(fname)
        st.session_state.ca = ca
        st.session_state.ca_path = fname


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
            with tempfile.NamedTemporaryFile(delete=False, suffix=".h5") as tmp:
                tmp.write(uploaded_file.getbuffer())
                tmp_path = tmp.name

            load_data(tmp_path)

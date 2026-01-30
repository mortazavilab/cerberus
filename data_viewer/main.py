import os
os.environ["STREAMLIT_SERVER_HEADLESS"] = "true"
os.environ["STREAMLIT_LOG_LEVEL"] = "debug"


import streamlit as st
import pandas as pd
import os
import sys
import argparse
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import matplotlib as mpl
import cerberus

SESSION_VARS = ['ca', 'ca_path', 'prev_triplet_set']

__version__ = '0.0.1'


d = os.path.dirname(__file__)

st.set_page_config(
    page_title="Cerberus & TriSET triplets viewer",
    page_icon=f"{d}/cerberus_logo.png",
    layout="wide",
    initial_sidebar_state="expanded"
)

def init_session_state():
    defaults = {c: None for c in SESSION_VARS}

    for k, v in defaults.items():
        st.session_state.setdefault(k, v)

def reset_data():
    for k in SESSION_VARS:
        st.session_state[k] = None


def main():
    init_session_state()

    st.sidebar.image(f"{d}/cerberus_logo.png", width=300)

    tab_landing, tab_simplex = st.tabs([
        "Data upload",
        "Simplex view",
    ])

    from tabs.landing_tab import render_landing_tab
    from tabs.simplex_tab import render_simplex_tab

    with tab_landing:
        render_landing_tab()


    with tab_simplex:
        if st.session_state.ca is None:
            st.info("No data loaded yet.")
        else:
            render_simplex_tab()

if __name__ == "__main__":
    main()

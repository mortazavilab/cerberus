import streamlit as st
import cerberus
import matplotlib.pyplot as plt
from utils import *

def simplex_options():

    st.markdown("### Simplex view controls")

    plotting_options = {}

    plotting_options['triplet_set'] = st.selectbox(
        "Triplet set",
        options=list(st.session_state.gnames.keys()),
        key="triplet_set"
    )


    # gene name picker
    plotting_options['gname'] = st.selectbox(
        "Gene name",
        options=st.session_state.gnames[
            plotting_options['triplet_set']
        ],
        key="gname",
        disabled=len(st.session_state.gnames[
            plotting_options['triplet_set']
        ]) <= 1
    )

    st.markdown("#### Plot type (select at least 1)")

    plotting_options['scatter'] = st.checkbox("Scatter")
    plotting_options['density'] = st.checkbox('Density')

    st.markdown("#### Marker style options")

    # marker size picker
    plotting_options['marker_size'] = st.selectbox(
        "Marker size",
        options=st.session_state.num_cols[
            plotting_options['triplet_set']
        ],
        key="marker_size"
    )

    # marker color pickser
    plotting_options['marker_color'] = st.selectbox(
        'Marker color',
        options=st.session_state.all_cols[
            plotting_options['triplet_set']
        ],
        key='marker_color'
    )

    plotting_options['legend'] = st.checkbox('Legend')


    # advanced options
    adv_expander = st.expander("Advanced")
    with adv_expander:

        st.markdown("#### Threshold options")

        # sector boundaries
        plotting_options['tss_thresh'] = st.slider(
                      "TSS-high threshold",
                      0.0, 1.0, 0.5)
        plotting_options['tes_thresh'] = st.slider(
                      "TES-high threshold",
                      0.0, 1.0, 0.5)
        plotting_options['spl_thresh'] = st.slider(
                      "Splicing-high threshold",
                      0.0, 1.0, 0.5)

        st.markdown("#### Other")

        # log continuous marker colors

        # change density kernel size (also default should be lower)
        plotting_options['density_scale'] = st.slider(
                      "Density resolution",
                      1, 50, 25,
                      help='Lower number = lower resolution')

    # actually make the plot
    run = st.button('Plot')

    return plotting_options, run

def render_simplex_plot(po):
    """
    Render plot
    """

    ca = st.session_state.ca

    if not po['scatter'] and not po['density']:
        st.write('Please check at least one of "Scatter" or "Density"')
        return

    # if gname and triplet_set:
    if po['gname']:
        gid = ca.triplets.loc[ca.triplets.gname==po['gname']].gid.values[0]
        # genecards_link = f'https://www.genecards.org/Search/Keyword?queryString=%20%5Ball%5D%20%20(%20%20{gid}%20%20)%20&advanced=true'
        st.header(f"{po['triplet_set']} triplets for {po['gname']} ({gid})")
    else:
        st.header(f"{po['triplet_set']} triplets")

    legend_scatter = legend_density = False
    if po['legend'] and po['scatter']: legend_scatter = True
    if po['legend'] and po['density']: legend_density = True

    # use colors in here - auto pulls from 'hex_color'
    # when marker color is either sample or sample_display
    if po['marker_color'] in ['sample', 'sample_display'] :
        temp = ca.triplets[['sample', 'hex_color']]
        temp = (
            temp
            .loc[temp['sample'].notnull()]
            .drop_duplicates()
            .set_index('sample')
        )
        cmap = temp.squeeze().to_dict()
    else: cmap = None

    df = ca.plot_simplex(
        subset={'source': po['triplet_set']},
        gene=po['gname'],
        scatter=po['scatter'],

        # density stuff
        density=po['density'],
        density_scale=po['density_scale'],

        log_density=True,
        density_cmap='Purples',

        # marker size stuff
        size=po['marker_size'],
        log_size=True,

        size_scale=.75,

        # marker color
        hue=po['marker_color'],
        cmap=cmap,

        sectors=True,
        legend=legend_scatter,
        density_cbar=legend_density,

        # sector boundaries
        sect_alpha=po['tss_thresh'],
        sect_beta=po['spl_thresh'],
        sect_gamma=po['tes_thresh']
    )

    st.pyplot(plt.gcf())
    plt.clf()
    plt.close('all')

    data_expander = st.expander("Data table")
    with data_expander:
        st.dataframe(df, hide_index=True)

def render_simplex_tab():
    """
    Original simplex dorito plot
    """
    with st.sidebar:
        plotting_options, run = simplex_options()

    if run:
        render_simplex_plot(plotting_options)

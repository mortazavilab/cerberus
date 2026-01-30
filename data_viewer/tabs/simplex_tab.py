import streamlit as st
import cerberus
import matplotlib.pyplot as plt
from utils import *

def render_simplex_tab():
    """
    Original simplex dorito plot
    """

    ca = st.session_state.ca

    with st.sidebar.expander("Simplex view controls", expanded=False):

        st.markdown("### Data subset options")

        triplet_set = st.selectbox(
            "Triplet set",
            options=ca.triplets.source.unique().tolist(),
            key="triplet_set"
        )

    if st.session_state.prev_triplet_set != triplet_set:

        temp_triplet_set = ca.triplets.loc[
            ca.triplets.source == triplet_set
        ]

        gnames = (
            temp_triplet_set.gname
                .dropna()
                .astype(str)
                .unique()
                .tolist()
        )

        gname = st.selectbox(
            "Gene name",
            options=sorted(gnames),
            key="gname",
            disabled=len(gnames) <= 1
        )

        num_cols = (
            temp_triplet_set
                .select_dtypes(include="number")
                .dropna(axis=1, how="all")
                .columns
                .tolist()
        )

        marker_size = st.selectbox(
            "Marker size",
            options=sorted(num_cols),
            key="marker_size"
        )

        # which triplet set did we just analyze
        st.session_state.prev_triplet_set = triplet_set



    # # options
    # with st.sidebar.expander("Simplex view controls",
    #     expanded=False):
    #     st.markdown("### Simplex view options")
    #     st.markdown("---")
    #
    #     ca = st.session_state.ca
    #
    #     st.markdown("#### Data subset options")
    #
    #     # triplet set picker
    #     triplet_set = st.selectbox(
    #         label='Triplet set',
    #         options=ca.triplets.source.unique().tolist()
    #     )
    #
    #     temp_triplet_set = ca.triplets.loc[ca.triplets.source==triplet_set]
    #
    #     # gene picker -- limit to genes in this triplet set -- precomputing
    #     # and storing in a dict would make this much quicker
    #     gnames = (temp_triplet_set
    #         .gname
    #         .dropna()
    #         .astype(str)
    #         .unique()
    #         .tolist()
    #     )
    #
    #     def_gname = None
    #     if len(gnames) > 1:
    #         gnames = sorted(gnames)
    #         if 'ELN' in gnames:
    #             def_gname = 'ELN'
    #         else:
    #             def_gname = gnames[0]
    #         disabled = False
    #     else:
    #         disabled = True
    #
    #     gname = st.selectbox(label='Gene name', options=gnames, disabled=disabled)
    #     # gname = check_multiselect(gname, 'gene name')
    #
    #     st.markdown("#### Plot type (select at least 1)")
    #     scatter = st.checkbox("Scatter")
    #     density = st.checkbox('Density')
    #
    #     st.markdown("#### Marker style options")
    #
    #
    #     # marker size
    #     num_cols = (
    #         temp_triplet_set
    #         .select_dtypes(include="number")
    #         .dropna(axis=1, how='all')
    #         .columns
    #         .tolist()
    #     )
    #     marker_size = st.selectbox(label='Marker size', options=sorted(num_cols))
    #     # marker_size = check_multiselect(marker_size, 'marker size')
    #
    #     all_cols = (
    #         temp_triplet_set
    #         .dropna(axis=1, how='all')
    #         .columns
    #         .tolist()
    #     )
    #     marker_color = st.selectbox(label='Marker color', options=all_cols)
    #     legend = st.checkbox('Legend')
    #
    #     # marker_color = check_multiselect(marker_color, 'marker color')
    #
    #     st.markdown("---")
    #     st.markdown("### Advanced")
    #
    #
    #     # adv_expander = st.expander("Advanced")
    #     # with adv_expander:
    #
    #     st.markdown("#### Threshold options")
    #
    #
    #     # sector boundaries
    #     tss_thresh = st.slider(
    #                   "TSS-high threshold",
    #                   0.0, 1.0, 0.5)
    #     tes_thresh = st.slider(
    #                   "TES-high threshold",
    #                   0.0, 1.0, 0.5)
    #     spl_thresh = st.slider(
    #                   "Splicing-high threshold",
    #                   0.0, 1.0, 0.5)
    #
    #     st.markdown("#### Other")
    #
    #     # log continuous marker colors
    #
    #     # change density kernel size (also default should be lower)
    #     density_scale = st.slider(
    #                   "Density resolution",
    #                   1, 50, 25,
    #                   help='Lower number = lower resolution')
    #
    #     # actually make the plot
    #     run = st.button('Plot')

    if run:
        if not scatter and not density:
            st.write('Please check at least one of "Scatter" or "Density"')
            run = False

    if run:

        # if gname and triplet_set:
        if gname:
            gid = ca.triplets.loc[ca.triplets.gname==gname].gid.values[0]
            # genecards_link = f'https://www.genecards.org/Search/Keyword?queryString=%20%5Ball%5D%20%20(%20%20{gid}%20%20)%20&advanced=true'
            st.header(f"{triplet_set} triplets for {gname} ({gid})")
        else:
            st.header(f"{triplet_set} triplets")

        legend_scatter = legend_density = False
        if legend and scatter: legend_scatter = True
        if legend and density: legend_density = True

        # use colors in here
        if marker_color in ['sample', 'sample_display'] :
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
            subset={'source': triplet_set},
            gene=gname,
            scatter=scatter,

            # density stuff
            density=density,
            # density_scale=density_scale,
            density_scale=10,

            log_density=True,
            density_cmap='Purples',

            # marker size stuff
            size=marker_size,
            log_size=True,

            size_scale=.75,

            # marker color
            hue=marker_color,
            cmap=cmap,

            sectors=True,
            legend=legend_scatter,
            density_cbar=legend_density,

            # sector boundaries
            sect_alpha=tss_thresh,
            sect_beta=spl_thresh,
            sect_gamma=tes_thresh
        )

        st.pyplot(plt.gcf())
        plt.clf()
        plt.close('all')

        data_expander = st.expander("Data table")
        with data_expander:
            st.dataframe(df, hide_index=True)

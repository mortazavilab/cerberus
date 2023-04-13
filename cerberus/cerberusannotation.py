from cerberus.cerberus import *
import swan_vis as swan
import sparse
import pandas as pd
import ternary
import matplotlib.pyplot as plt
import matplotlib as mpl
from sklearn import preprocessing
from collections import defaultdict
import matplotlib.ticker as tck


class CerberusAnnotation():
    def __init__(self):
        self.tss = None
        self.tes = None
        self.ic = None

        self.tss_map = None
        self.tes_map = None
        self.t_map = None

        self.tss_sources = []
        self.tes_sources = []
        self.ic_sources = []
        self.all_sources = []

        self.triplets = pd.DataFrame()
        self.triplet_sources = []

    def get_sources(self, df):
        """
        Get list of source names that were used to generate the triplet features
        in df

        Parameters:
            df (pandas DataFrame): DataFrame of TSSs, TESs, or ICs

        Returns:
            s (list of str): List of sources from df
        """
        s = [s.split(',') for s in df.source.unique().tolist()]
        s = list(set(flatten(s)))
        return s

    def get_all_sources(self):
        """
        Get sources that are present in tss, ic, and tes. Essentially gives
        the list of sources that are from full-length transcripts.
        """
        all_s = list(set(self.tss_sources)&\
                     set(self.tes_sources)&\
                     set(self.ic_sources))
        return all_s

    def set_tss(self, tss):
        self.tss = tss
        self.tss_sources = self.get_sources(self.tss)
        self.all_sources = self.get_all_sources()

    def set_tes(self, tes):
        self.tes = tes
        self.tes_sources = self.get_sources(self.tes)
        self.all_sources = self.get_all_sources()

    def set_ic(self, ic):
        self.ic = ic
        self.ic_sources = self.get_sources(self.ic)
        self.all_sources = self.get_all_sources()


################################################################################
############### Pertaining to triplet computations #############################
################################################################################

    def add_sg_info(self, df, sg=None):
        """
        Add information about transcripts from corresponding SwanGraph

        Parameters:
            df (pandas DataFrame): DF of transcripts
            sg (swan_vis SwanGraph): SwanGraph with transcript IDs that match
                those in the CerberusAnnotation

        Returns:
            df (pandas DataFrame): DF with SwanGraph t_df info added
        """
        if sg is not None:
            temp = sg.t_df[['gid', 'gname']].copy(deep=True)
            df.rename({'gene_id': 'gid'}, axis=1, inplace=True)
            temp.rename({'gid':'gene_id'}, axis=1, inplace=True)
            temp = add_stable_gid(temp)
            temp.reset_index(drop=True)
            temp.drop_duplicates(inplace=True)
            temp.rename({'gene_id': 'gid'}, axis=1, inplace=True)
            df = df.merge(temp, how='left', on='gid')

        return df

    def add_triplets(self, df, source=None):
        """
        Add a triplet dataframe to the existing triplet register
        in `self.triplets`.

        Parameters:
            df (pandas DataFrame): DF output from `get_source_triplets`, or
                `get_subset_triplets`
            source (str): Name of source. Required if `source` column not in df
        """

        if not isinstance(df, pd.DataFrame):
            df = pd.DataFrame()

        # make sure the source or sources haven't been added already
        if len(self.triplet_sources) > 0:
            if 'source' in df.columns:
                if source in self.triplet_sources:
                    raise ValueError('Source {} already present. Choose a new name.'.format(source))

            else:
                sources = df.source.unique().tolist()
                for source in sources:
                    if source in self.triplet_sources:
                        raise ValueError('Source {} already present. Choose a new name.'.format(source))

        if 'source' not in df.columns:
                df['source'] = source
        self.triplets = pd.concat([self.triplets, df], axis=0)

        if 'source' in df.columns:
            self.triplet_sources += df.source.unique().tolist()
        else:
            self.triplet_sources.append(source)

    def get_source_triplets(self, **kwargs):
        """
        Compute triplets for all sources. No filtering is done based on
        expression etc.

        Returns:
            trip_df (pandas DataFrame): DF indicating how many tss, tes, and ic
                features are found in each source
        """
        sources = self.all_sources + ['all']
        trip_df = pd.DataFrame()
        for source in sources:

            if source != 'all':

                # make sure that source is in each df
                if source not in self.all_sources:
                    raise ValueError('Source {} not in cerberus annotation'.format(source))

                tss = self.tss.loc[self.tss.source.str.contains(source)].copy(deep=True)
                tes = self.tes.loc[self.tes.source.str.contains(source)].copy(deep=True)
                ic = self.ic.loc[self.ic.source.str.contains(source)].copy(deep=True)

            else:
                tss = self.tss.copy(deep=True)
                tes = self.tes.copy(deep=True)
                ic = self.ic.copy(deep=True)

            # compute total # of unique ends / ics for each df
            tss = count_instances(tss, 'tss')
            tes = count_instances(tes, 'tes')
            ic = count_instances(ic, 'ic')

            df = tss
            df = df.merge(tes, on='gene_id', how='outer')
            df = df.merge(ic, on='gene_id', how='outer')
            df['source'] = source

            # add n_isos, if possible
            if source in self.t_map.source.unique().tolist() or source == 'all':
                if source != 'all':
                    iso = self.t_map.loc[self.t_map.source == source].copy(deep=True)
                else:
                    iso = self.t_map.copy(deep=True)
                iso = iso[['transcript_id', 'gene_id']].drop_duplicates()
                iso = iso.groupby('gene_id').count().reset_index()
                iso.rename({'transcript_id': 'n_iso'}, axis=1, inplace=True)
                df = df.merge(iso, how='left', on='gene_id')

            trip_df = pd.concat([trip_df, df])

        # add splicing ratio and compute simplex coords
        trip_df = compute_splicing_ratio(trip_df)
        trip_df = compute_simplex_coords(trip_df)
        trip_df = assign_sector(trip_df)

        # rename gid col
        trip_df.rename({'gene_id': 'gid'}, axis=1, inplace=True)

        trip_df = self.add_sg_info(trip_df, **kwargs)

        return trip_df

    def get_expressed_triplets(self, sg, obs_col, min_tpm,
                               source=None, subset=None):
        """
        Compute triplets for each sample for a given metadata
        column and minimum tpm

        Parameters:
            sg (swan_vis SwanGraph): SwanGraph w/ abundance information for each
                transcript
            obs_col (str): Column in `sg.adata.obs` to compute detection on
                the basis of
            min_tpm (float): Minimum TPM to consider a transcript expressed
            subset (pandas DataFrame): Table with `obs_col` category and which
                tids to consider from that category

        Returns:
            trips (pandas DataFrame): DF w/ # tss, ic, tes expressed in each
                category
        """

        gb = obs_col

        # get sample, transcript level TPM values and format
        t_df = swan.calc_tpm(sg.adata, obs_col='dataset', how='max')
        t_df = t_df.sparse.to_dense()
        t_df = t_df.merge(sg.adata.obs[gb], how='left', left_index=True, right_index=True)
        t_df.reset_index(drop=True, inplace=True)
        t_df.set_index(gb)
        t_df = t_df.groupby(gb, observed=True).max()
        t_df = t_df.T

        # get transcript detection table and format
        t_df = (t_df >= min_tpm)

        # names of each dataset
        dataset_cols = t_df.columns.tolist()

        # do some reformatting, add gene_id
        t_df.reset_index(inplace=True)
        t_df.rename({'index': 'tid'}, axis=1, inplace=True)
        t_df['gene_id'] = get_gid_from_tids(t_df.tid.tolist())

        # limit only to isoforms that are in the subset set
        if isinstance(subset, pd.DataFrame):
            temp = t_df.copy(deep=True)
            temp['gene_id'] = get_stable_gid(temp, 'gene_id')
            temp = temp.melt(id_vars=['tid', 'gene_id'], var_name=gb, value_name='detected')
            subset['gene_id'] = get_stable_gid(subset, 'gid')
            subset = subset[['gene_id', gb, 'tid']]
            temp = temp.merge(subset, how='inner', on=['gene_id', gb, 'tid'])
            temp = temp.pivot(index=['tid', 'gene_id'], columns=[gb], values='detected').reset_index().fillna(False)
            t_df = t_df[['tid']]
            temp = temp.merge(t_df, how='left', on='tid')
            t_df = temp.copy(deep=True)

        t_df['gid_stable'] = get_stable_gid(t_df, col='gene_id')
        t_df.drop('gene_id', axis=1, inplace=True)
        t_df.rename({'gid_stable': 'gene_id'}, axis=1, inplace=True)

        # count the unique tss, ic, tes, and isoforms from each expressed
        # isoform in each sample
        trip_df = pd.DataFrame()
        for c in dataset_cols:
            # ensures that we're only looking at isoforms >= 1 TPM in each sample
            temp = self.get_subset_triplets(t_df.loc[t_df[c], 'tid'].tolist())
            temp[gb] = c
            trip_df = pd.concat([trip_df, temp])

        # add gene tpm and format
        if sg is not None:
            g_df = swan.calc_tpm(sg.gene_adata, obs_col=gb)
            g_df = g_df.sparse.to_dense()
            g_df = g_df.T
            g_df.reset_index(inplace=True)
            g_df.rename({'index': 'gene_id'}, axis=1, inplace=True)
            g_df = add_stable_gid(g_df)
            g_df.set_index('gene_id', inplace=True)
            g_df = g_df.melt(var_name=gb, value_name='gene_tpm', ignore_index=False)
            g_df.set_index(gb, inplace=True, append=True)
            trip_df.set_index(['gene_id', gb], inplace=True)
            trip_df = trip_df.merge(g_df, how='left', left_index=True, right_index=True)

        # compute splicing ratio and simplex coords
        trip_df = compute_splicing_ratio(trip_df)
        trip_df = compute_simplex_coords(trip_df)
        trip_df = assign_sector(trip_df)

        trip_df.reset_index(inplace=True)
        trip_df = self.add_sg_info(trip_df, sg)

        if source:
            trip_df['source'] = source

        return trip_df

    def get_subset_triplets(self, tids, source=None, **kwargs):
        """
        Get the set of triplets for each gene based on a list of input isoforms

        Parameters:
            tids (list of str): List of tids
            source (str): Name to use for triplets source
        """
        df = get_feats_from_tids(tids)
        df.drop('transcript_triplet', axis=1, inplace=True)
        df = df.groupby('gene_id').nunique()
        df.rename({'tss_id': 'n_tss',
                     'ic_id': 'n_ic',
                     'tes_id': 'n_tes',
                     'transcript_id': 'n_iso'}, axis=1, inplace=True)
        df.reset_index(inplace=True)

        # compute splicing ratio and simplex coords
        df = compute_splicing_ratio(df)
        df = compute_simplex_coords(df)
        df = assign_sector(df)

        df = self.add_sg_info(df, **kwargs)

        if source:
            df['source'] = source

        return df

################################################################################
########################### Writer #############################################
################################################################################
    def write(self, fname):
        """
        Write cerberus annotation object to an h5 file.

        Parameters:
            fname (str): Name / path of h5 file to save to
        """
        write_h5(self.ic,
                 self.tss,
                 self.tes,
                 fname,
                 self.tss_map,
                 self.tes_map,
                 self.t_map,
                 self.triplets)

###############################################################################
########################### Plotting ############################################
###############################################################################
    def plot_simplex(self,

                # arguments for the subset of data shown
                subset=None, # used in main
                gene=None, # used in main

                # arguments for what plotting style is used
                density=False, # used in main, passed to scatter
                sectors=False, # used in main
                scatter=True, # used in main
                jitter=False, # used in main

                density_scale=1, # used in main
                density_cbar=False, # used in main
                density_vmax=None, # not used in main but needs explicit passing
                title=None, # used in main
                top='splicing_ratio', # used in main

                size_scale=1,

                fname=None,
                **kwargs):
        """
        Plot a simplex plot from CerberusAnnotation triplets

        Parameters:
            subset (dict of lists): List mapping counts column names
                to values in said columns to include in the data
                ie {'source': 'v40', 'sample': ['brain', 'heart']}
            gene (str): Gene name to subset on
            density (bool): Plot density plot
            sectors (bool): Include lines that bound each sector
            scatter (bool): Plot scatter plot
            jitter (bool): Jitter scattered points
            density_scale (float): Number of triangles to plot per side
            density_cbar (bool): Include reference colorbar
            density_vmax (float): Value to serve as max. in density colorbar
            title (str): Title to display on plot
            size_scale (float): How large to make figure
            fname (str): File to save figure as

            **kwargs: Arguments to pass to TODO

        Returns:
            temp (pandas DataFrame): Points plotted in this plot
        """

        scale = True

        # figure size handling
        size_dict = {}

        # heights / widths
        size_dict['height'] = 10
        size_dict['width1'] = 13.6
        size_dict['width2'] = 17.2
        size_dict['width3'] = 10.716
        # size_dict['width3'] = 11.2


        # linewidths
        size_dict['linewidth'] = 3
        size_dict['linewidth2'] = 2

        # marker sizes
        size_dict['big_source_marker'] = 200
        size_dict['big_marker'] = 100
        size_dict['small_marker'] = 20
        size_dict['marker_max'] = 300
        size_dict['marker_min'] = 10

        # jitter
        size_dict['sigma'] = (1/250)*scale*density_scale

        # ticks
        size_dict['tick_linewidth'] = 2

        # fontsizes
        size_dict['tick_size'] = 30
        size_dict['cbar_size'] = 16
        size_dict['title_size'] = 36.6
        size_dict['label_size'] = 36.6
        size_dict['cbar_tick_size'] = 16

        size_dict['title_pad'] = 50

        for key, item in size_dict.items():
            # area scale
            if 'marker' in key:
                size_dict[key] = item*size_scale*size_scale
            # linear scale
            else:
                size_dict[key] = item*size_scale


        # offsets don't need to be scaled
        size_dict['left_label_offset'] = 0.15
        size_dict['right_label_offset'] = 0.18
        size_dict['bottom_label_offset'] = 0.1
        size_dict['tick_offset'] = 0.025

        #### subset dataset and transform numbers as needed ####
        temp = self.triplets.copy(deep=True)

        # plot settings
        mpl.rcParams['font.family'] = 'Arial'
        mpl.rcParams['pdf.fonttype'] = 42

        # if we have a gene name, limit to those entries
        if gene:
            temp, gene = subset_df_on_gene(temp, gene)

        # if we have a list of allowed sources, limit to those entries
        if subset:
            temp = subset_df(temp, subset)

        # scale and assign which columns to use
        c = dict()
        if scale:
            # temp['total'] = temp.n_tss+temp.n_tes+temp[top]
            # temp['tss_ratio'] = temp.n_tss/temp.total
            # temp['tes_ratio'] = temp.n_tes/temp.total
            # temp['spl_ratio'] = temp[top]/temp.total

            # temp = compute_simplex_coords(temp, top)

            c['a'] = 'tss_ratio'
            c['b'] = 'spl_ratio'
            c['c'] = 'tes_ratio'
        else:
            c['a'] = 'n_tss'
            c['b'] = top
            c['c'] = 'n_tes'

        if scale == True:
            scale = 1
            mult = 0.2
        else:
            scale = max_pts(temp, c)

        # density
        if density:
            if 'hue' in kwargs:
                hue = kwargs['hue']
                if temp[hue].dtype.name == 'object':
                    pad = 0.1
                else:
                    pad = 0.0
            else:
                pad = 0.1
            # pad = pad*size_scale
            figure, tax, temp = self.density_dorito(temp, c,
                                     size_dict,
                                     density_scale=density_scale,
                                     pad=pad,
                                     density_cbar=density_cbar,
                                     density_vmax=density_vmax,
                                     **kwargs)
            scale = density_scale
            if density_cbar:
                figure.set_size_inches(size_dict['width1'],size_dict['height'])
            else:
                figure.set_size_inches(size_dict['width3'], size_dict['height'])

        # if we're jittering, adjust the points for each thing
        if jitter:
            temp, c = self.jitter_dorito(temp, c, density_scale, size_dict)

        # figure layout parameters
        fontsize = 18
        offset = 0.1
        mult = scale/5

        # if we don't already have a fig and axis from density,
        # make one
        if not density:
            figure, tax = ternary.figure(scale=scale, permutation='210')
            figure.set_facecolor('white')

        # plot gridlines below the scatterplot
        tax.gridlines(linewidth=size_dict['linewidth'],
                multiple=mult,
                color='white',
                zorder=1,
                linestyle=None)

        # scatter
        if scatter:
            figure, tax = self.scatter_dorito(temp, c,
                            figure, tax,
                            density_cbar,
                            size_dict,
                            **kwargs)

        # sectors
        if sectors:
            self.line_dorito(figure, tax, scale, size_dict, **kwargs)

        # title handler
        fdict = {}
        if not title:
            if gene:
                # title = '$\it{}$\n'.format(gene)
                title = '{}'.format(gene)
                fdict = {'fontstyle': 'italic'}
            else:
                title = ''
        else:
            if gene:
                title = '{} $\it{}$'.format(title, gene)
            else:
                title = '{}'.format(title)

        tax.set_title(title, fontsize=size_dict['title_size'],
                      fontdict=fdict, pad=size_dict['title_pad'])
        tax.boundary(linewidth=size_dict['linewidth2'], c='#e5ecf6')
        labels = ['{:.1f}'.format(n) for n in np.arange(0, 1.2, .2)]
        tax.ticks(ticks=labels,
            axis='lbr', linewidth=size_dict['tick_linewidth'], multiple=mult,
            tick_formats="%.1f", offset=size_dict['tick_offset'],
            fontsize=size_dict['tick_size'])

        tax.clear_matplotlib_ticks()
        tax.get_axes().axis('off')
        tax.set_background_color('#e5ecf6')

        if top == 'splicing_ratio':
            # top_label = 'Splicing ratio $\\beta$'
            top_label = 'Splicing ratio'
        elif top == 'intron_chain':
            # top_label = 'Intron chains $\\delta$'
            top_label = 'Intron chains'

        # tax.left_corner_label('# TSSs $\\alpha$', fontsize=fontsize)
        # tax.top_corner_label(top_label, fontsize=fontsize)
        # tax.right_corner_label('# TESs $\\gamma$', fontsize=fontsize)
        tax.left_axis_label('TSS', fontsize=size_dict['label_size'],
                            offset=size_dict['left_label_offset'])
        tax.right_axis_label(top_label, fontsize=size_dict['label_size'],
                             offset=size_dict['right_label_offset'])
        tax.bottom_axis_label('TES', fontsize=size_dict['label_size'],
                              offset=size_dict['bottom_label_offset'])

        figure.set_facecolor('white')

        # save figure
        if fname:
            # plt.savefig(fname, dpi=500, bbox_inches='tight')
            tax.savefig(fname, dpi=500, bbox_inches='tight')

        return temp

    def scatter_dorito(self,
                     counts,
                     c,
                     figure,
                     tax,
                     density,
                     size_dict,
                     order_marker_sizes=False,
                     hue=None,
                     cmap='magma',
                     marker_style=None,
                     mmap=None,
                     size=None,
                     log_size=False,
                     alpha=1,
                     legend=True,
                     **kwargs):
        """
        Parameters
          counts (pandas DataFrame): Subset of self.triplets being plotted
          c (dict of str): Dictionary of column names to plot as a, b, c
              indexed by 'a', 'b', 'c'
          figure (matplotlib figure): Figure to plot on
          tax (matplotlib ternary axes): Axes to plot on
          density (bool): Plot density
          size_dict (dict of floats): Sizes for different plotted features
          order_marker_sizes (bool): Order markers based on their sizes such
            that the smallest markers are on top
            hue (str): Column in self.triplets to color each point by
            cmap (str): Colormap to use for hue
            marker_style (str): Column in self.triplets to use to plot different
                shapes of markers
            mmap (dict): Dictionary mapping marker styles to different column values
                in self.triplets[marker_style]
            size (str): Column in self.triplets to use to control size of markers
            log_size (bool): Log marker sizes
            alpha (float): Alpha value for plotted markers
            legend (bool): Include legend

            **kwargs

        Returns:
            figure :
            tax :
        """

        def scale_col(points, counts, col,
                      min_size=10,
                      max_size=300,
                      log=False,
                      how='color'):
            if log:
                log_col = '{}_log'.format(col)
                counts[log_col] = np.log10(counts[col])
                col = log_col
            vals = counts[[col]]
            max_val = vals[col].max()
            min_val = vals[col].min()
            min_max_scaler = preprocessing.MinMaxScaler(feature_range=(min_size, max_size))
            vals = min_max_scaler.fit_transform(vals)
            max_scaled = max(vals)
            min_scaled = min(vals)

            # replace nans w/ big marker size
            # vals = [size_dict['big_marker'] if np.isnan(v) else v for v in vals]
            vals = [size_dict['big_source_marker'] if np.isnan(v) else v for v in vals]

            return vals, min_val, max_val, min_scaled, max_scaled

        # defaults
        points = [(x[0], x[1], x[2]) for x in zip_pts(counts, c)]
        labels = ['' for i in range(len(points))]
        hue_type = None
        figsize = (size_dict['width3'],size_dict['height'])
        colors = '#e78424'
        if len(points) < 60:
            sizes = [size_dict['big_marker'] for i in range(len(points))]
        else:
            sizes =  [size_dict['small_marker'] for i in range(len(points))]
        markers = 'o'
        vmin = 0
        vmax = 1
        plotted = False

        # get color
        if hue:

            # categorical
            if counts[hue].dtype.name == 'object':
                hue_type = 'cat'
                colors = counts[hue].map(cmap).tolist()
                if marker_style:
                    counts.loc[counts[hue].isnull(), hue] = counts.loc[counts[hue].isnull(), marker_style]
                    labels = counts[hue].tolist()

            # continuous
            else:
                hue_type = 'cont'
                colors, abs_min, abs_max, vmin, vmax = scale_col(points, counts, hue,
                                                                 min_size=size_dict['marker_min'],
                                                                 max_size=size_dict['marker_max'])

        # get sizes
        if size:
            sizes,_,_,_,_ = scale_col(points, counts, size, log=log_size,
                                                             min_size=size_dict['marker_min'],
                                                             max_size=size_dict['marker_max'])

        # get marker styles
        if marker_style:
            markers = counts[marker_style].map(mmap).fillna('o').tolist()
        else:
            markers = ['o' for i in range(len(counts.index))]

        if hue_type == 'cat' and density: figsize = (size_dict['width1'], size_dict['height'])
        elif hue_type == 'cat' and not density: figsize = (size_dict['width3'], size_dict['height'])
        elif hue_type == 'cont' and density: figsize = (size_dict['width2'], size_dict['height'])
        elif hue_type == 'cont' and not density: figsize = (size_dict['width1'], size_dict['height'])
        elif density: figsize = (size_dict['width1'], size_dict['height'])
        figure.set_size_inches(figsize[0], figsize[1])

        def sort_points(points, colors, sizes, labels, markers):
            # https://stackoverflow.com/questions/6618515/sorting-list-based-on-values-from-another-list
            # attr = [points, colors, labels, markers]
            # for att in attr:
            def sort_attr(attr, sizes):
                attr = [a for (a,size) in sorted(zip(attr, sizes),
                        key=lambda pair: pair[1],
                        reverse=True)]
                return attr

            points = sort_attr(points, sizes)
            colors = sort_attr(colors, sizes)
            labels = sort_attr(labels, sizes)
            markers = sort_attr(markers, sizes)
            sizes = sort_attr(sizes, sizes)
            return points, colors, sizes, labels, markers

        # actual scatter call
        if hue_type == 'cat':

            # sort points based on their sizes so that smallest will be visible
            # on top of the bigger ones
            if order_marker_sizes:
                points, colors, sizes, labels, markers = sort_points(points,
                                                                     colors,
                                                                     sizes,
                                                                     labels,
                                                                     markers)

            for point, color, size, label, marker in zip(points, colors, sizes, labels, markers):
                # if marker == 'x':
                    # import pdb; pdb.set_trace()
                if marker == 'x':
                    linewidths=None
                else:
                    linewidths=0
                tax.scatter([point], vmin=vmin, vmax=vmax,
                        s=size, c=color,
                        marker=marker,label=label,
                        alpha=alpha, zorder=3,
                        linewidths=linewidths,
                        edgecolors=None)
        else:
            tax.scatter(points, vmin=vmin, vmax=vmax,
                        s=sizes, c=colors, cmap=cmap, marker=markers,
                        alpha=alpha, zorder=3,
                        linewidths=0,
                        edgecolors=None)

        # legend handling
        if hue_type == 'cat' and legend:

            if density: x = 1.6
            else: x = 1.4
            tax.legend(bbox_to_anchor=(x, 1.05),
                    loc='upper right', prop={'size': 14})

            # fix marker size
            ax = tax.get_axes()
            lgnd = ax.get_legend()
            for handle in lgnd.legendHandles:
                handle._sizes = [size_dict['big_marker']]

        # colorbar handling
        if hue_type == 'cont':
            ax = tax.get_axes()
            norm = plt.Normalize(vmin=abs_min, vmax=abs_max)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm._A = []
            cb = plt.colorbar(sm, ax=ax, pad=0.1)
            for t in cb.ax.get_yticklabels():
                 t.set_fontsize(size_dict['cbar_size'])
            if hue == 'tss' or hue == 'tes':
                cb.set_label('# {}s'.format(hue.upper()), size=16)
            elif hue == 'intron_chain':
                cb.set_label('# {}s'.format(hue), size=16)

        return figure, tax

    def density_dorito(self,
                    counts,
                    c,
                    size_dict,
                    density_scale=20,
                    density_cmap='viridis',
                    density_vmax=None,
                    log_density=False,
                    pad=0.15,
                    density_cbar=False,
                    **kwargs):
        """
        Plot the density of a dataset on a ternary plot
        From here: https://github.com/marcharper/python-ternary/issues/81

        Parameters:
            counts
            c
            scale
            cmap

        Returns:
            fig
            tax
            counts (pandas DataFrame): Counts, scaled by factor used
        """
        hm_dict = defaultdict(int)
        log = log_density
        scale = density_scale
        cmap = density_cmap
        for i in range(0, scale+1):
            for j in range(0, scale+1):
                for k in range(0, scale+1):
                    if i+j+k == scale:
                        temp = counts.copy(deep=True)
                        if i != scale:
                            temp = temp.loc[(temp.tss_ratio*scale>=i)&(temp.tss_ratio*scale<i+1)]
                        else:
                            temp = temp.loc[(temp.tss_ratio*scale>=i)&(temp.tss_ratio*scale<=i+1)]
                        if j != scale:
                            temp = temp.loc[(temp.spl_ratio*scale>=j)&(temp.spl_ratio*scale<j+1)]
                        else:
                            temp = temp.loc[(temp.spl_ratio*scale>=j)&(temp.spl_ratio*scale<=j+1)]
                        n = len(temp.index)
                        hm_dict[i,j] += n

        # for key, item in hm_dict.items():
        #     blop = 0
        #     if item > 0:
        #         print(key)
        #         print(item)
        #         print()
        #         blop+=1
        #     if blop > 3:
        #         break
        # print(hm_dict)

        # log values if necessary
        if log:
            for key, item in hm_dict.items():
                hm_dict[key] = np.log2(item+1)
            if density_vmax:
                density_vmax = np.log2(density_vmax+1)

        # double checking stuff
        df = pd.DataFrame.from_dict(hm_dict, orient='index')
        df['i'] = [b[0] for b in df.index.tolist()]
        df['j'] = [b[1] for b in df.index.tolist()]
        df.rename({0:'val'}, axis=1, inplace=True)

        figure, tax = ternary.figure(scale=scale, permutation='210')
        tax.heatmap(hm_dict, colorbar=False, style='t',
                    adj_vlims=True, cmap=cmap, vmax=density_vmax)

        # scale according to chosen resolution
        for key in c.keys():
            counts[c[key]] = counts[c[key]]*scale

        # colorbar - hacked together by broken ternary code
        if density_cbar:
            ax = tax.get_axes()
            flat = []
            for key, item in hm_dict.items():
                flat.append(item)
            min_val = min(flat)
            max_val = max(flat)

            if density_vmax:
                max_val = density_vmax

            print(min_val)
            print(max_val)
            norm = plt.Normalize(vmin=min_val, vmax=max_val)
            sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
            sm._A = []

            def exp_format(x,pos):
                x = int(x)
                return r'$2^{{{}}}$'.format(x)

            if not log:
                cb = plt.colorbar(sm, ax=ax, pad=pad)
            else:
                cb = plt.colorbar(sm, ax=ax, pad=pad,
                    format=tck.FuncFormatter(exp_format))

            for t in cb.ax.get_yticklabels():
                t.set_fontsize(size_dict['cbar_tick_size'])
            if not log:
                cb.ax.set_yticklabels([])
                cb.ax.set_yticks([])

            cb.set_label('Density', size=size_dict['cbar_size'])

        return figure, tax, counts

    def jitter_dorito(self, counts, c, scale, size_dict):
        """
            Parameters:
                counts
                c
                scale

            Returns
                counts
                c
        """

        # figure out how much to jitter by
        sigma = size_dict['sigma']
        for d in c.keys():
            d_jitter = '{}_jitter'.format(d)
            counts[d_jitter] = counts[c[d]].apply(lambda x: np.random.normal(x, sigma))
            c[d] = d_jitter

        return counts, c

    def line_dorito(self,
                figure,
                tax,
                scale,
                size_dict,
                sect_alpha=0.5,
                sect_beta=0.5,
                sect_gamma=0.5,
                **kwargs):
        c_dict, _ = get_sector_colors()

        # scale
        alpha = sect_alpha*scale
        beta = sect_beta*scale
        gamma = sect_gamma*scale

        linewidth = size_dict['linewidth']

        # splicing line
        tax.horizontal_line(beta, linewidth=linewidth,
                         color=c_dict['splicing'],
                         linestyle='--')

        # tss
        tax.right_parallel_line(alpha, linewidth=linewidth,
                            color=c_dict['tss'],
                            linestyle='--')

        # tes
        tax.left_parallel_line(gamma, linewidth=linewidth,
                            color=c_dict['tes'],
                            linestyle='--')

# helper functions

################################################################################
############################### Plotting helpers ###############################
################################################################################
def zip_pts(df, c):
    return zip(df[c['a']], df[c['b']], df[c['c']])

def max_pts(df, c):
    return max(df[c['a']].max(), df[c['b']].max(), df[c['c']].max())

def rm_color_cats(c_dict, order, cats):
    if cats:
        keys = c_dict.keys()
        pop_list = []
        for key in keys:
            if key not in cats:
                pop_list.append(key)
        for p in pop_list:
            del c_dict[p]
        order = [o for o in order if o in cats]
    return c_dict, order

def get_sector_colors(cats=None):
    tss = '#56B4E9'
    tes = '#E69F00'
    splicing = '#CC79A7'
    simple = '#e5ecf6'
    c_dict = {'tss': tss,
              'splicing': splicing,
              'tes': tes,
              'simple': simple,
              'mixed': '#b7b7b7'}
    order = ['tss', 'splicing', 'tes', 'mixed', 'simple']

    c_dict, order = rm_color_cats(c_dict, order, cats)
    return c_dict, order


################################################################################
############################### Other helpers ###############################
################################################################################
def flatten(list_of_lists):
    """
    Flatten a list of lists into a 1D array

    Parameters:
        list_of_lists (list of lists): List of lists

    Returns:
        (list): Flattened array
    """
    if len(list_of_lists) == 0:
        return list_of_lists
    if isinstance(list_of_lists[0], list):
        return flatten(list_of_lists[0]) + flatten(list_of_lists[1:])
    return list_of_lists[:1] + flatten(list_of_lists[1:])

def count_instances(df, mode):
    count_col = 'n_{}'.format(mode)
    df = df[['gene_id', 'Name']].groupby('gene_id').count().reset_index()
    df = df.rename({'Name': count_col}, axis=1)
    df[count_col] = df[count_col].astype(int)
    return df

def subset_df_on_gene(df, gene):
    """
    Subset a dataframe based on gene name or id

    Parameters:
        df (pandas DataFrame): DF to subset
        gene (str): Gene ID or name to subset on

    Returns:
        df (pandas DataFrame): DF subset on gene
        gene (str): Gene name
    """
    temp = df.loc[df.gname == gene]
    if temp.empty:
        temp = df.loc[df.gid == gene]
        gene = df.gname.values[0]

    df = temp.copy()
    return temp, gene

def subset_df(df, subset):
    """
    Subset a dataframe based on meeting all of the conditions given

    Parameters:
        df (pandas DataFrame): DF to subset
        subset (dictionary of str): Dictionary mapping {column name :
                                    acceptable values}

    Returns:
        df (pandas DataFrame): DF subset on the dictionary values
    """
    for col, val in subset.items():
        if type(val) != list:
            val = [val]
        df = df.loc[df[col].isin(val)]
    return df

def assign_sector(df):
    """
    Assign each gene triplet a sector that reflects its behavior

    Parameters:
        df (pandas DataFrame): DataFrame with simplex coordinates and n_iso

    Returns:
        df (pandas DataFrame): DataFrame with sector listed
    """
    df['sector'] = 'simple'

    df.loc[df.tss_ratio > 0.5, 'sector'] = 'tss'
    df.loc[df.tes_ratio > 0.5, 'sector'] = 'tes'
    df.loc[df.spl_ratio > 0.5, 'sector'] = 'splicing'

    # mixed genes
    df.loc[(df.sector=='simple')&(df.n_iso>1), 'sector'] = 'mixed'
    return df

def compute_simplex_coords(df, spl_col='splicing_ratio'):
    """
    Compute simplex coordinates from gene triplets

    Parameters:
        df (pandas DataFrame): DF holding values for n tss, tes, and
            column to approximate amount of splicing
        spl_col (str): {'splicing_ratio', 'n_ic'}

    Returns:
        df (pandas DataFrame): DF w/ simplex coordinates values in columns
            ['tss_ratio', 'spl_ratio', 'tes_ratio']
    """
    df['total'] = df.n_tss+df.n_tes+df[spl_col]
    df['tss_ratio'] = df.n_tss/df.total
    df['tes_ratio'] = df.n_tes/df.total
    df['spl_ratio'] = df[spl_col]/df.total
    df.drop('total', axis=1, inplace=True)

    return df

def simplex_dist_pts(a, b, how='js'):
    """
    Compute the distance between two points on a simplex

    Parameters:
        a (np.array): Coords of pt a
        b (np.array): Coords of pt b
        how (str): How to calculate distance. {'jensenshannon'}

    Returns:
        dist (float): Distance b/w a and b using distance metric
    """
    if how == 'js':
        dist = scipy.spatial.distance.jensenshannon(a,b)
    return dist

def simplex_dist(x, suff_a=None, suff_b=None, **kwargs):
    """
    From a series, compute the distance between two points

    Parameters:
        x (pandas Series): Row-wise pandas series
        suff_a (str): Suffix for point set a (ie if column is 'tss_ratio_det'),
            provide '_det' as `suff_a`
        suff_b (str): Suffix for point set b

    Returns:
        dist (list of float): List of distances between point set a and point
            set b for each point
    """
    def get_pt(x, suff):
        tss = 'tss_ratio'
        ic = 'spl_ratio'
        tes = 'tes_ratio'
        if suff:
            tss = '{}{}'.format(tss, suff)
            ic = '{}{}'.format(ic, suff)
            tes = '{}{}'.format(tes, suff)
        pt = [x[tss], x[ic], x[tes]]
        return pt

    a = get_pt(x, suff_a)
    b = get_pt(x, suff_b)
    dist = simplex_dist_pts(a,b, **kwargs)

    return dist

def compute_splicing_ratio(df):
    """
    Compute the splicing ratio from an input df with columns
    `n_ic`, `n_tss`, and `n_tes`

    Parameters:
        df (pandas DataFrame): DF with columns for
            n tss, ic, tes

    Returns:
        df (pandas DataFrame): DF w/ splicing ratio added
    """
    df['splicing_ratio'] = df.n_ic/((df.n_tes+df.n_tss)/2)
    return df

def read(h5):
    ca = CerberusAnnotation()
    ic, tss, tes, tss_map, tes_map, m, triplets = read_h5(h5, as_pyranges=False)
    ca.set_tss(tss)
    ca.set_tes(tes)
    ca.set_ic(ic)
    ca.tss_map = tss_map
    ca.tes_map = tes_map
    ca.t_map = m
    ca.add_triplets(triplets)

    return ca

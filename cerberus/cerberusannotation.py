from cerberus.cerberus import *
import swan_vis as swan
import sparse
import pandas as pd

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
            sg (swan_vis SwanGraph): SwanGraph holding abundance information
                but also other metadata about each gene
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

        # add splicing ratio
        trip_df = compute_splicing_ratio(trip_df)

        # rename gid col
        trip_df.rename({'gene_id': 'gid'}, axis=1, inplace=True)

        trip_df = self.add_sg_info(trip_df, **kwargs)
        # if sg is not None:
        #     temp = sg.t_df[['gid', 'gname']].copy(deep=True)
        #     temp.rename({'gid':'gene_id'}, axis=1, inplace=True)
        #     temp = add_stable_gid(temp)
        #     temp.reset_index(drop=True)
        #     temp.drop_duplicates(inplace=True)
        #     temp.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        #     trip_df = trip_df.merge(temp, how='left', on='gid')

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
        t_df = swan.calc_tpm(sg.adata, obs_col='dataset')
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
            temp = t_df.melt(id_vars=['tid', 'gene_id'], var_name=gb, value_name='detected')
            temp['gene_id'] = get_stable_gid(temp, 'gene_id')
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

        trip_df = compute_splicing_ratio(trip_df)
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
        """
        df = get_feats_from_tids(tids)
        df.drop('transcript_triplet', axis=1, inplace=True)
        df = df.groupby('gene_id').nunique()
        df.rename({'tss_id': 'n_tss',
                     'ic_id': 'n_ic',
                     'tes_id': 'n_tes',
                     'transcript_id': 'n_iso'}, axis=1, inplace=True)
        df.reset_index(inplace=True)
        df = compute_splicing_ratio(df)

        df = self.add_sg_info(df, **kwargs)

        if source:
            df['source'] = source

        # # add the gene id
        # temp = self.sg.t_df[['gid', 'gname']].copy(deep=True)
        # df.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        # temp.rename({'gid':'gene_id'}, axis=1, inplace=True)
        # temp = cerberus.add_stable_gid(temp)
        # temp.reset_index(drop=True)
        # temp.drop_duplicates(inplace=True)
        # temp.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        # df = df.merge(temp, how='left', on='gid')

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

# helper functions
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

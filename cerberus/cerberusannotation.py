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

        self.triplets = None

    def get_sources(self, df):
        s = [s.split(',') for s in df.source.unique().tolist()]
        s = list(set(flatten(s)))
        return s

    def get_all_sources(self):
        """
        Get sources that are present in tss, ic, and tes
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


#     def subset_on_source(self, source):
#         self.tss = self.tss.loc[self.tss.source.str.contains(source)].copy(deep=True)
#         self.tes = self.tes.loc[self.tes.source.str.contains(source)].copy(deep=True)
#         self.ic = self.ic.loc[self.ic.source.str.contains(source)].copy(deep=True)

#         return self
    def set_sg(self, sg):
        self.sg = sg

    ## triplet computations
    def get_source_triplets(self, sg):
        """
        Compute triplets for all sources
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
            trip_df = pd.concat([trip_df, df])

        # add splicing ratio
        trip_df = compute_splicing_ratio(trip_df)

        # rename gid col
        trip_df.rename({'gene_id': 'gid'}, axis=1, inplace=True)

        # add the gene id
        temp = self.sg.t_df[['gid', 'gname']].copy(deep=True)
        temp.rename({'gid':'gene_id'}, axis=1, inplace=True)
        temp = add_stable_gid(temp)
        temp.reset_index(drop=True)
        temp.drop_duplicates(inplace=True)
        temp.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        trip_df = trip_df.merge(temp, how='left', on='gid')

        return trip_df

    def get_expressed_triplets(self, obs_col, min_tpm, subset=None):
        """
        Compute triplets for each sample for a given metadata
        column and minimum tpm
        """

        gb = obs_col

        # get sample, transcript level TPM values and format
        t_df = swan.calc_tpm(self.sg.adata, obs_col='dataset')
        t_df = t_df.sparse.to_dense()
        t_df = t_df.merge(self.sg.adata.obs[gb], how='left', left_index=True, right_index=True)
        t_df.reset_index(drop=True, inplace=True)
        t_df.set_index(gb)
        t_df = t_df.groupby(gb, observed=True).max()
        t_df = t_df.T

        # get transcript detection table and format
        t_df = (t_df >= min_tpm)

        # some formatting
        cols = ['tss_id', 'ic_id', 'tes_id', 'gid']
        dataset_cols = t_df.columns.tolist()
        t_df = t_df.merge(self.sg.t_df[cols], how='left', left_index=True, right_index=True)
        t_df.rename({'gid': 'gene_id'}, axis=1, inplace=True)
        t_df.reset_index(inplace=True)
        t_df.rename({'index': 'tid'}, axis=1, inplace=True)

        # limit only to isoforms that are in the subset set
        if isinstance(subset, pd.DataFrame):
            temp = t_df.melt(id_vars=['tid', 'gene_id'], var_name=gb, value_name='detected')
            temp['gene_id'] = get_stable_gid(temp, 'gene_id')

            subset['gene_id'] = get_stable_gid(subset, 'gid')
            subset = subset[['gene_id', 'sample', 'tid']]
            # pdb.set_trace()

            temp = temp.merge(subset, how='inner', on=['gene_id', gb, 'tid'])

            temp = temp.pivot(index=['tid', 'gene_id'], columns=[gb], values='detected').reset_index().fillna(False)
            t_df = t_df[['tid', 'tss_id', 'tes_id', 'ic_id']]
            temp = temp.merge(t_df, how='left', on='tid')
            t_df = temp.copy(deep=True)

        # count the unique tss, ic, tes, and isoforms from each expressed
        # isoform in each sample
        trip_df = pd.DataFrame()
        for c in dataset_cols:
            keep_cols = ['tid', 'tss_id', 'ic_id', 'tes_id', 'gene_id']
            temp = t_df.loc[t_df[c]].copy(deep=True)
            temp = temp[keep_cols]
            temp = temp.groupby('gene_id').nunique()
            temp.rename({'tss_id': 'n_tss',
                         'ic_id': 'n_ic',
                         'tes_id': 'n_tes',
                         'tid': 'n_iso'}, axis=1, inplace=True)
            temp[gb] = c
            temp.set_index(gb, inplace=True, append=True)
            trip_df = pd.concat([trip_df, temp])

        # add gene tpm and format
        g_df = swan.calc_tpm(self.sg.gene_adata, obs_col=gb)
        g_df = g_df.sparse.to_dense()
        g_df = g_df.T
        g_df.reset_index(inplace=True)
        g_df.rename({'index': 'gene_id'}, axis=1, inplace=True)
        g_df = add_stable_gid(g_df)
        g_df.set_index('gene_id', inplace=True)
        g_df = g_df.melt(var_name=gb, value_name='gene_tpm', ignore_index=False)
        g_df.set_index(gb, inplace=True, append=True)

        trip_df = trip_df.merge(g_df, how='left', left_index=True, right_index=True)
        trip_df = compute_splicing_ratio(trip_df)

        trip_df.reset_index(inplace=True)

        # add the gene id
        temp = self.sg.t_df[['gid', 'gname']].copy(deep=True)
        trip_df.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        temp.rename({'gid':'gene_id'}, axis=1, inplace=True)
        temp = add_stable_gid(temp)
        temp.reset_index(drop=True)
        temp.drop_duplicates(inplace=True)
        temp.rename({'gene_id': 'gid'}, axis=1, inplace=True)
        trip_df = trip_df.merge(temp, how='left', on='gid')

        return trip_df


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
    df['splicing_ratio'] = df.n_ic/((df.n_tes+df.n_tss)/2)
    return df

def read(h5):
    c_annot = CerberusAnnotation()
    ic, tss, tes, tss_map, tes_map, m = read_h5(h5, as_pyranges=False)
    c_annot.set_tss(tss)
    c_annot.set_tes(tes)
    c_annot.set_ic(ic)
    c_annot.tss_map = tss_map
    c_annot.tes_map = tes_map
    c_annot.t_map = m

    return c_annot

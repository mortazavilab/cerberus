import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def make_end_df(c,s,st,e,n, source,mode):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Start', 'End', 'Name']
    var = [c,s,st,e,n]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    # add source
    df['source'] = source

    df = format_end_df(df)

    # get end # and gene id
    if any(df.Name.isnull()):
        df['gene_id'] = np.nan
        df[mode] = np.nan
    else:
        df['gene_id'] = df.Name.str.split('_', expand=True)[0]
        df[mode] = df.Name.str.split('_', expand=True)[1]

    # get arbitrary unique ids
    df['id'] = [i for i in range(len(df.index))]

    return df#

def format_end_df(df):
    sort_cols = ['Chromosome', 'Start', 'End', 'Strand']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Start', 'End', 'Strand', 'Name', 'source']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df

def test_agg_2_ends_1(print_dfs=True):
    """
    Test agg_2_ends w/ and w/o end adding
    """

    def get_test(mode='tss'):
        # example has
        # - entries that overlap
        # - entries that don't overlap but are within a certain distance
        # - entries that are unique to either
        # - entries that overlap but aren't using the same gene_id / strand (these will be
        #      equivalent situtations b/c the gene ids for things on different strands will always
        #      differ
        # - entries in bed2 that are non continuous but overlap the same
        #      entry in bed1

        slack = 20

        n = 4
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [1, 200, 100, 300]
        e = [15, 250, 110, 340]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4']
        source = 'v1'
        bed1 = make_end_df(c,s,st,e,n, source, mode)
        bed1 = pr.PyRanges(bed1)

        n = 5
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [5, 20, 120, 500, 200]
        e = [10, 25, 140, 550, 250]
        n = ['gene1_1', 'gene1_4', 'gene1_2', 'gene1_3', 'gene2_1']
        source = 'v2'
        bed2 = make_end_df(c,s,st,e,n, source, mode)
        bed2 = pr.PyRanges(bed2)

        return bed1, bed2

    def get_ctrl(add=True):

        mode = 'tss'

        n = 6
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [1,200,100,300, 200,500]
        e = [15,250,110,340, 250,550]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1', 'gene1_5']
        source = ['v1,v2','v1','v1,v2','v1', 'v2','v2']
        df = make_end_df(c,s,st,e,n, source,mode)

        # convert a few dtypes
        df['Strand'] = df['Strand'].astype('category')
        df['Chromosome'] = df['Chromosome'].astype('category')

        # if we're not adding new ends
        if not add:
            df = df.loc[df.source != 'v2']

        # remove unnecessary columns
        df.drop(['gene_id', 'id', mode], axis=1, inplace=True)

        # fix ids
        df.reset_index(drop=True, inplace=True)

        return df

    tests = [True, False]
    for add_ends in tests:

        slack = 20
        mode = 'tss'
        sort_cols = ['Chromosome', 'Strand', 'gene_id', 'Start', 'End']
        order = ['Chromosome', 'Start', 'End', 'Strand', 'Name',
                 'gene_id', 'source', mode]
        bed1, bed2 = get_test()
        df = agg_2_ends(bed1, bed2,
                        strand=True,
                        gid=True,
                        slack=slack,
                        add_ends=add_ends,
                        mode=mode)

        test = format_end_df(df)
        ctrl = get_ctrl(add=add_ends)

        if print_dfs:
            print('test')
            print(test)
            print(test.index)
            print(test.dtypes)
            print('ctrl')
            print(ctrl)
            print(ctrl.index)
            print(ctrl.dtypes)

        pd.testing.assert_frame_equal(ctrl, test, check_like=True)

        assert len(ctrl.index) == len(test.index)

def test_agg_2_ends_2(print_dfs=True):
    """
    Test agg_2_ends w/ and w/o end adding
    """

    def get_test(mode='tss'):

        # adding a bed file that doesn't have strand or gid info
        # entries in bed2 that are duplicated based on
        # - strandedness
        # - gene id
        # entries in bed2 where both of them overlap the same region
        # in bed1

        slack = 20

        n = 6
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[-1] = '-'
        st = [1, 200, 100, 300, 260, 260]
        e = [15, 250, 110, 340, 290, 290]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1', 'gene3_1']
        source = 'v1'
        bed1 = make_end_df(c,s,st,e,n, source, mode)
        bed1 = pr.PyRanges(bed1)

        n = 5
        c = ['1' for i in range(n)]
        s = [np.nan for i in range(n)]
        st = [5, 11, 120, 500, 200]
        e = [10, 21, 140, 550, 250]
        n = [np.nan for i in range(n)]
        source = 'v2'
        bed2 = make_end_df(c,s,st,e,n, source, mode)
        bed2 = pr.PyRanges(bed2)

        return bed1, bed2

    def get_ctrl(add=True):

        mode = 'tss'

        n = 6
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[-1] = '-'
        st = [1,200,100,300, 260,260]
        e = [15,250,110,340, 290,290]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1', 'gene3_1']
        source = ['v1,v2','v1,v2','v1,v2','v1', 'v1,v2', 'v1,v2']
        df = make_end_df(c,s,st,e,n, source,mode)

        # convert a few dtypes
        df['Strand'] = df['Strand'].astype('category')
        df['Chromosome'] = df['Chromosome'].astype('category')

        # if we're not adding new ends
        if not add:
            df = df.loc[df.source != 'v2']

        # remove unnecessary columns
        df.drop(['gene_id', 'id', mode], axis=1, inplace=True)

        # fix ids
        df.reset_index(drop=True, inplace=True)

        return df

    strand = False
    gid = False
    add_ends = False

    slack = 20
    mode = 'tss'
    sort_cols = ['Chromosome', 'Strand', 'gene_id', 'Start', 'End']
    order = ['Chromosome', 'Start', 'End', 'Strand', 'Name',
             'gene_id', 'source', mode]
    bed1, bed2 = get_test()
    df = agg_2_ends(bed1, bed2,
                    strand=strand,
                    gid=gid,
                    slack=slack,
                    add_ends=add_ends,
                    mode=mode)
    test = format_end_df(df)
    ctrl = get_ctrl(add=add_ends)

    if print_dfs:
        print('test')
        print(test)
        print(test.index)
        print(test.dtypes)
        print('ctrl')
        print(ctrl)
        print(ctrl.index)
        print(ctrl.dtypes)

    pd.testing.assert_frame_equal(ctrl, test, check_like=True)
    assert len(ctrl.index) == len(test.index)

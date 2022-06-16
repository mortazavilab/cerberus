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
        #       bed1: gene1_1,
        #       bed2: gene1_1
        # - entries that don't overlap but are within a certain distance
        #       bed1: gene1_1
        #       bed2: gene1_4
        # - entries that are unique to either
        #       bed1: gene1_2, gene1_4
        #       bed2: gene1_3, gene2_1, gene1_5, gene4_1
        # - entries that overlap but aren't using the same gene_id / strand (these will be
        #      equivalent situtations b/c the gene ids for things on different strands will always
        #      differ
        #       bed1: gene1_2
        #       bed2: gene2_1
        # - entries in bed2 that are non continuous but overlap the same
        #      entry in bed1
        #       bed1: gene1_1
        #       bed2: gene1_1, gene1_4
        # - entries in bed1 that the same entry in bed2 maps to that are all the same gene
        #       bed1: gene3_1, gene3_2
        #       bed2: gene3_1
        # - entries in bed1 that the same entry in bed2 maps to that are NOT from the same gene;
        #    need to create new end
        #       bed1: gene3_1, gene3_2
        #       bed2: gene4_1
        # - entry in bed1 that multiple entries in bed2 overlap
        #       bed1: gene5_1
        #       bed2: gene5_1, gene5_2

        slack = 20

        n = 7
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [1, 200, 100, 300, 700, 760, 850]
        e = [15, 250, 110, 340, 750, 770, 870]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene3_1', 'gene3_2', 'gene5_1']
        source = 'v1'
        bed1 = make_end_df(c,s,st,e,n, source, mode)
        bed1 = pr.PyRanges(bed1)

        n = 9
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [5, 20, 120, 500, 200, 750, 750, 860, 850]
        e = [10, 25, 140, 550, 250, 775, 775, 870, 855]
        n = ['gene1_1', 'gene1_4', 'gene1_2', 'gene1_3', 'gene2_1', 'gene3_1', 'gene4_1', 'gene5_1', 'gene5_2']
        source = 'v2'
        bed2 = make_end_df(c,s,st,e,n, source, mode)
        bed2 = pr.PyRanges(bed2)

        return bed1, bed2

    def get_ctrl(add=True):

        mode = 'tss'

        n = 10
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [1,200,100,300,700,760,850, 200,500,750]
        e = [15,250,110,340,750,770,870, 250,550,775]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene3_1', 'gene3_2', 'gene5_1', 'gene2_1', 'gene1_5', 'gene4_1']
        source = ['v1,v2','v1','v1,v2','v1','v1,v2','v1,v2','v1,v2', 'v2','v2','v2']
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

    def get_ctrl_source_map(add=True):
        # n = 9
        # c = ['1' for i in range(n)]
        # s = ['+' for i in range(n)]
        # st = [1, 200, 100, 300,  5, 20, 120, 500, 200]
        # e = [15, 250, 110, 340,  10, 25, 140, 550, 250]
        # n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4',
        #      'gene1_1', 'gene1_1', 'gene1_3', 'gene1_5', 'gene2_1']
        # source = ['v1' for i in range(4)]+['v2' for i in range(5)]
        # mode = 'tss'
        # df = make_end_df(c,s,st,e,n, source,mode)
        # df = df.loc[df.source == 'v2']
        # df.drop(['gene_id', mode, 'id'], axis=1, inplace=True)

        n = 10
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        st = [5, 20, 120, 500, 200, 750, 750, 750, 860, 850]
        e = [10, 25, 140, 550, 250, 775, 775, 775, 870, 855]
        n = ['gene1_1', 'gene1_1', 'gene1_3', 'gene1_5', 'gene2_1', 'gene3_1', 'gene3_2', 'gene4_1', 'gene5_1', 'gene5_1']
        source = 'v2'
        mode = 'tss'
        df = make_end_df(c,s,st,e,n, source,mode)
        df = df.loc[df.source == 'v2']
        df.drop(['gene_id', mode, 'id'], axis=1, inplace=True)
        df['Strand'] = df['Strand'].astype('category')
        df['Chromosome'] = df['Chromosome'].astype('category')

        if not add:
            new_names = ['gene2_1', 'gene1_5', 'gene4_1']
            new_inds = df.loc[df.Name.isin(new_names)].index
            df.loc[new_inds, 'Name'] = np.nan

        return df

    tests = [True, False]
    for add_ends in tests:

        slack = 20
        mode = 'tss'
        sort_cols = ['Chromosome', 'Strand', 'gene_id', 'Start', 'End']
        order = ['Chromosome', 'Start', 'End', 'Strand', 'Name',
                 'gene_id', 'source', mode]
        bed1, bed2 = get_test()

        df, m_source = agg_2_ends(bed1, bed2,
                        strand=True,
                        gid=True,
                        slack=slack,
                        add_ends=add_ends,
                        mode=mode)

        test = format_end_df(df)
        test_m = format_end_df(m_source)
        ctrl = get_ctrl(add=add_ends)
        ctrl_m = get_ctrl_source_map(add=add_ends)

        if print_dfs:
            print('test')
            print(test)
            print(test.index)
            print(test.dtypes)
            print('ctrl')
            print(ctrl)
            print(ctrl.index)
            print(ctrl.dtypes)
            print('test source map')
            print(test_m)
            print(test_m.index)
            print(test_m.dtypes)
            print('ctrl source map')
            print(ctrl_m)
            print(ctrl_m.index)
            print(ctrl_m.dtypes)

        pd.testing.assert_frame_equal(ctrl, test, check_like=True)
        assert len(ctrl.index) == len(test.index)

        pd.testing.assert_frame_equal(ctrl_m, test_m, check_like=True)
        assert len(ctrl_m.index) == len(test_m.index)

def test_agg_2_ends_2(print_dfs=True):
    """
    Test agg_2_ends w/ and w/o end adding
    """

    def get_test(mode='tss'):

        # adding a bed file that doesn't have strand or gid info
        # entries in bed2 that are duplicated based on
        # - strandedness
        # - gene id
        #     bed1: gene2_1, gene3_1
        #     bed2: (200-250)
        # entries in bed2 where both of them overlap the same region
        # in bed1
        #     bed1: gene2_1, gene3_1
        #     bed2: (200-250)
        # - entries in bed1 that the same entry in bed2 maps to that are all the same gene
        #       bed1: gene4_1, gene4_2
        #       bed2: (750-775)


        slack = 20

        n = 8
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[5] = '-'
        st = [1, 200, 100, 300, 260, 260, 700, 760]
        e = [15, 250, 110, 340, 290, 290, 750, 770]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1', 'gene3_1', 'gene4_1', 'gene4_2']
        source = 'v1'
        bed1 = make_end_df(c,s,st,e,n, source, mode)
        bed1 = pr.PyRanges(bed1)

        n = 6
        c = ['1' for i in range(n)]
        s = [np.nan for i in range(n)]
        st = [5, 11, 120, 500, 200, 750]
        e = [10, 21, 140, 550, 250, 775]
        n = [np.nan for i in range(n)]
        source = 'v2'
        bed2 = make_end_df(c,s,st,e,n, source, mode)
        bed2 = pr.PyRanges(bed2)

        return bed1, bed2

    def get_ctrl(add=True):

        mode = 'tss'

        n = 8
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[5] = '-'
        st = [1,200,100,300, 260,260, 700,760]
        e = [15,250,110,340, 290,290, 750,770]
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1', 'gene3_1', 'gene4_1', 'gene4_2']
        source = ['v1,v2','v1,v2','v1,v2','v1', 'v1,v2', 'v1,v2', 'v1,v2', 'v1,v2']
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

    def get_ctrl_source_map(add=False):

        n = 9
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[5] = '-'
        s[6] = np.nan
        st = [5, 11, 120, 200, 200, 200, 500, 750, 750]
        e = [10, 21, 140, 250, 250, 250, 550, 775, 775]
        n = ['gene1_1', 'gene1_1', 'gene1_3', 'gene1_2', 'gene2_1', 'gene3_1', np.nan, 'gene4_1', 'gene4_2']
        source = 'v2'
        mode = 'tss'
        df = make_end_df(c,s,st,e,n, source,mode)
        df = df.loc[df.source == 'v2']
        df.drop(['gene_id', mode, 'id'], axis=1, inplace=True)
        df['Strand'] = df['Strand'].astype('category')
        df['Chromosome'] = df['Chromosome'].astype('category')

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
    df, m_source = agg_2_ends(bed1, bed2,
                    strand=strand,
                    gid=gid,
                    slack=slack,
                    add_ends=add_ends,
                    mode=mode)
    test = format_end_df(df)
    test_m = format_end_df(m_source)
    ctrl = get_ctrl(add=add_ends)
    ctrl_m = get_ctrl_source_map()

    if print_dfs:
        print('test')
        print(test)
        print(test.index)
        print(test.dtypes)
        print('ctrl')
        print(ctrl)
        print(ctrl.index)
        print(ctrl.dtypes)
        print('test source map')
        print(test_m)
        print(test_m.index)
        print(test_m.dtypes)
        print('ctrl source map')
        print(ctrl_m)
        print(ctrl_m.index)
        print(ctrl_m.dtypes)

    pd.testing.assert_frame_equal(ctrl, test, check_like=True)
    assert len(ctrl.index) == len(test.index)

    pd.testing.assert_frame_equal(ctrl_m, test_m, check_like=True)
    assert len(ctrl_m.index) == len(test_m.index)

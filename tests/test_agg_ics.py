import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def make_ics_df(c,s,co,n, source,nov):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Coordinates', 'Name']
    var = [c,s,co,n]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    df['gene_id'] = df.Name.str.split('_', expand=True)[0]
    df['ic'] = df.Name.str.split('_', expand=True)[1]
    df['source'] = source

    if nov:
        df['novelty'] = nov

    df.ic = df.ic.astype(int)
    df = format_ics_df(df)

    return df

def format_ics_df(df):
    sort_cols = ['Chromosome', 'Strand', 'Coordinates']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Strand', 'Coordinates', 'Name', 'gene_id', 'ic', 'source', 'novelty']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df

def test_agg_2_ics_1(print_dfs=True):

    def get_test():
        # new df has one matching ic
        # new df has a novel ic that needs a new number
        # new df has a novel gene

        n = 4
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3', '11-12-13-14']
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene2_1']
        source = 'v1'
        ref = 'Known'
        df1 = make_ics_df(c,s,co,n, source, ref)

        n = 3
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-4-3', '1-2-6', '11-12-13']
        n = ['gene1_1', 'gene1_2', 'gene2_1']
        source = 'v2'
        ref = 'Novel'
        df2 = make_ics_df(c,s,co,n, source, ref)

        return df1, df2

    def get_ctrl():
        n = 6
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3', '1-2-6', '11-12-13', '11-12-13-14']
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_2', 'gene2_1']
        source = ['v1', 'v1,v2', 'v1', 'v2', 'v2', 'v1']
        ref = ['Known' if 'v1' in s else 'Novel' for s in source]
        df = make_ics_df(c,s,co,n, source, ref)
        # df.drop(['gene_id', 'ic'], axis=1, inplace=True)

        return df

    df1, df2 = get_test()
    ctrl = get_ctrl()

    test = agg_2_ics(df1, df2)
    test = format_ics_df(test)

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

def test_get_novelty(print_dfs=True):

    def get_test():
        n = 11
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3', '11-12-13-14', '2-3-4-5',
              '1-2-6', '11-12-13', '1-3', '-', '3-5', '11-14'] # NNC, ISM, NNC, Monexonic, NNC, NIC
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene2_1', 'gene2_3',
             'gene1_4', 'gene2_2', 'gene1_5', 'gene1_6', 'gene2_4', 'gene2_5']
        source = ['v1', 'v1,v2', 'v1', 'v1', 'v1',
                  'v2', 'v2', 'v2', 'v2', 'v2', 'v2']
        ref = ['Known' if 'v1' in s else 'Novel' for s in source]
        df = make_ics_df(c,s,co,n, source, ref)
        return df

    def get_ctrl():
        n = 11
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3', '11-12-13-14', '2-3-4-5',
              '1-2-6', '11-12-13', '1-3', '-', '3-5', '11-14'] # NNC, ISM, NNC, Monexonic, NNC, NIC
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene2_1', 'gene2_3',
             'gene1_4', 'gene2_2', 'gene1_5', 'gene1_6', 'gene2_4', 'gene2_5']
        source = ['v1', 'v1,v2', 'v1', 'v1', 'v1',
                  'v2', 'v2', 'v2', 'v2', 'v2', 'v2']
        ref = ['Known', 'Known', 'Known', 'Known', 'Known',
               'NNC', 'ISM', 'NNC', 'Monoexonic', 'NNC', 'NIC']
        df = make_ics_df(c,s,co,n, source, ref)
        return df

    df = get_test()
    ctrl = get_ctrl()

    test = get_ic_novelty(df)
    test = format_ics_df(test)

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

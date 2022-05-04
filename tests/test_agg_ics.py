import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def make_ics_df(c,s,co,n):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Coordinates', 'Name']
    var = [c,s,co,n]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    df = format_ics_df(df)

    df['gene_id'] = df.Name.str.split('_', expand=True)[0]
    df['ic'] = df.Name.str.split('_', expand=True)[1]
    df.ic = df.ic.astype(int)

    return df

def format_ics_df(df):
    sort_cols = ['Chromosome', 'Strand', 'Coordinates']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Strand', 'Coordinates', 'Name', 'gene_id', 'ic']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df

def test_agg_2_ics_1(print_dfs=True):

    def get_test():
        # new df has one matching ic
        # new df has a novel ic that needs a new number
        # new df has a novel gene

        n = 3
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3']
        n = ['gene1_1', 'gene1_2', 'gene1_3']
        df1 = make_ics_df(c,s,co,n)

        n = 3
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-4-3', '1-2-6', '11-12-13']
        n = ['gene1_1', 'gene1_2', 'gene2_1']
        df2 = make_ics_df(c,s,co,n)

        return df1, df2

    def get_ctrl():
        n = 5
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        co = ['1-2-3', '1-4-3', '5-2-3', '1-2-6', '11-12-13']
        n = ['gene1_1', 'gene1_2', 'gene1_3', 'gene1_4', 'gene2_1']
        df = make_ics_df(c,s,co,n)
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

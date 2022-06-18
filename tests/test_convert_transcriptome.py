import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def make_end_gtf_df(c,s,st,e,g,t):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Start', 'End', 'gene_id', 'transcript_id']
    var = [c,s,st,e,g,t]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    df = format_end_df(df)

    return df

def make_end_df(c,s,st,e,n, source,mode, add_id=True):
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
    if add_id:
        df['id'] = [i for i in range(len(df.index))]

    return df

def format_end_df(df):
    sort_cols = ['Chromosome', 'Start', 'End', 'Strand']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Start', 'End', 'Strand', 'Name', 'gene_id', 'transcript_id', 'start', 'source']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df

def make_ics_gtf_df(c,s,co,g,t):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Coordinates', 'gene_id', 'transcript_id']
    var = [c,s,co,g,t]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    df = format_ics_df(df)

    return df

def make_ics_map_df(t, ic, i):
    df = pd.DataFrame()
    cols = ['transcript_id', 'ic', 'ic_id']
    var = [t, ic, i]
    for col, var in zip(cols, var):
        df[col] = var
    return df

def make_ics_df(c,s,co,n, source):
    df = pd.DataFrame()
    cols = ['Chromosome', 'Strand', 'Coordinates', 'Name']
    var = [c,s,co,n]
    for col, var in zip(cols, var):
        if type(var) == list:
            df[col] = var

    df = format_ics_df(df)

    df['gene_id'] = df.Name.str.split('_', expand=True)[0]
    df['ic'] = df.Name.str.split('_', expand=True)[1]
    df['source'] = source
    df.ic = df.ic.astype(int)

    return df

def format_ics_df(df):
    sort_cols = ['Chromosome', 'Strand', 'Coordinates']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Strand', 'Coordinates', 'Name', 'gene_id', 'transcript_id', 'ic', 'source']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df



def test_merge_ends(print_dfs=True):
    # cases where
    # - ends overlap existing ends as well as gene id
    # - ends do not overlap existing ends and we need to choose closest
    #    - and the closest is downstream but we need to choose upstream
    #    - fwd and rev strands, tss and tes (t5, t6)
    # - ends overlap existing ends but gene ids don't match
    # - ends overlap but are not on same strand
    # - equidistant matches?

    for mode in ['tss', 'tes']:

        # annot
        n = 6
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[5] = '-'
        st = [1, 200, 100, 300, 425, 583]
        e = [2, 201, 101, 300, 426, 584]
        g = ['gene1', 'gene1', 'gene2', 'gene3', 'gene5', 'gene6']
        t = ['t1', 't2', 't3', 't4', 't5', 't6']
        annot_df = make_end_gtf_df(c,s,st,e,g,t)
        annot_df = pr.PyRanges(annot_df)

        # ref
        n = 9
        c = ['1' for i in range(n)]
        s = ['+' for i in range(n)]
        s[4] = '-'
        s[8:9] = '-'
        st = [1, 150, 90, 200, 290, 410, 427, 570, 590]
        e = [40, 160, 98, 250, 340, 420, 440, 580, 600]
        n = ['gene1_1', # entry in annot intersects
             'gene1_2', # entry in annot needs to choose closest
             'gene2_1', # entry in annot matches but wrong gene
             'gene3_1', #
             'gene4_1', # entry overlaps on the wrong strand
             'gene5_1', # entries that is upstream (tss) but not closest, +
             'gene5_2',
             'gene6_1', # entries that is upstream (tss) but not closest, -
             'gene6_2']
        source = 'test'
        ref_df = make_end_df(c,s,st,e,n, source, mode, add_id=False)
        ref_df = pr.PyRanges(ref_df)

        # ctrl
        t = ['t1', 't2', 't3', 't4', 't5', 't6']
        i = ['gene1_1', 'gene1_2', 'gene2_1', 'gene3_1', 'gene5_1', 'gene6_2']
        ends = ['1', '2', '1', '1', '1', '2']
        # needs to be downstream one if tes
        if mode == 'tes':
            i[4] = 'gene5_2'
            ends[4] = '2'
            i[5] = 'gene6_1'
            ends[5] = '1'
        ctrl = pd.DataFrame()
        ctrl['transcript_id'] = t
        ctrl['{}_id'.format(mode)] = i
        ctrl[mode] = ends


        test = merge_ends(annot_df, ref_df, mode)

        def order_df(df):
            df.sort_values(by='transcript_id', inplace=True)
            df.reset_index(drop=True, inplace=True)
            return df

        ctrl = order_df(ctrl)
        test = order_df(test)

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

# merge ics
def test_merge_ics(print_dfs=True):

    # ics to annotate
    n = 3
    c = ['1' for i in range(n)]
    s = ['+' for i in range(n)]
    co = ['1-2-3', '1-2-3', '4-5-6']
    g = ['gene1', 'gene1', 'gene2']
    t = ['t1', 't2', 't3']
    annot_df = make_ics_gtf_df(c,s,co,g,t)


    # ref ics
    n = 3
    c = ['1' for i in range(n)]
    s = ['+' for i in range(n)]
    co = ['1-2-3', '4-5-6', '4-5-6']
    n = ['gene1_1', 'gene1_2', 'gene2_1']
    source = 'test'
    ref_df = make_ics_df(c,s,co,n, source)

    # ctrl
    t = ['t1', 't2', 't3']
    ic = [1, 1, 1]
    i = ['gene1_1', 'gene1_1', 'gene2_1']
    ctrl = make_ics_map_df(t, ic, i)

    test = merge_ics(annot_df, ref_df)

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

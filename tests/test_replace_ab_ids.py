import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def test_map_transcripts(print_dfs=True):
    # map file
    otid = ['1', '2', '3', '4']
    otname = ['1n', '2n', '3n', '4n']
    tid = ['1[1,1,1]', '2[1,1,1]', '3[1,1,1]', '4[1,1,1]']
    tname = ['1n[1,1,1]', '2n[1,1,1]', '3n[1,1,1]', '4n[1,1,1]']
    m_df = pd.DataFrame()
    m_df['original_transcript_id'] = otid
    m_df['original_transcript_name'] = otname
    m_df['transcript_id'] = tid
    m_df['transcript_name'] = tname

    # transcript file
    tid = ['2', '3', '4', '1']
    tname = ['2n', '3n', '4n', '1n']
    df = pd.DataFrame()
    df['annot_transcript_id'] = tid
    df['annot_transcript_name'] = tname
    df['tid_2'] = tid

    # ctrl
    otid = ['2', '3', '4', '1']
    tname = ['2n[1,1,1]', '3n[1,1,1]', '4n[1,1,1]', '1n[1,1,1]']
    tid = ['2[1,1,1]', '3[1,1,1]', '4[1,1,1]', '1[1,1,1]']
    ctrl = pd.DataFrame()
    ctrl['tid_2'] = otid
    ctrl['annot_transcript_id'] = tid
    ctrl['annot_transcript_name'] = tname

    test = map_transcripts(df, m_df, 'annot_transcript_name', 'annot_transcript_id')

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


def make_talon_ab(t,ta,n,i,e,l,ab):
    df = pd.DataFrame()
    cols = ['gene_ID', 'annot_gene_id', 'annot_gene_name',
                   'gene_novelty', 'annot_transcript_name',
                   'annot_transcript_id', 'transcript_ID']
    for c in cols:
        df[c] = t
    df['transcript_ID'] = ta
    df['transcript_novelty'] = n
    df['ISM_subtype'] = i
    df['n_exons'] = e
    df['length'] = l
    for i, a in enumerate(ab):
        c = 'ab{}'.format(i)
        df[c] = a

    return df

def test_agg_ab(print_dfs=True):
    # input
    t = ['1', '2', '2', '3', '3']
    ta = ['5', '6', '7', '8', '9']
    n = ['Known', 'Known', 'NIC', 'ISM', 'ISM']
    i = [None, None, None, 'Prefix', 'Both']
    e = [1,1,3,1,1]
    l = [1,1,5,1,1]
    ab = [[0,1,1,2,2],
          [3,0,1,1,2]]
    df = make_talon_ab(t,ta,n,i,e,l,ab)

    # ctrl
    t = ['1', '2', '3']
    ta = ['5', '6,7', '8,9']
    n = ['Known', 'Known', 'ISM']
    i = [None, None, 'Both']
    e = [1,2,1]
    l = [1,3,1]
    ab = [[0,2,4],
          [3,1,3]]
    ctrl = make_talon_ab(t,ta,n,i,e,l,ab)

    test = agg_ab(df)

    col_order = ctrl.columns.tolist()
    col_order.sort()
    test = test[col_order]
    ctrl = ctrl[col_order]
    test.sort_values(by='annot_transcript_id', inplace=True)
    ctrl.sort_values(by='annot_transcript_id', inplace=True)
    test.reset_index(inplace=True, drop=True)
    ctrl.reset_index(inplace=True, drop=True)
    ctrl.n_exons = ctrl.n_exons.astype(float)
    ctrl.length = ctrl.length.astype(float)


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

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

    # get id
    df['{}_id'.format(mode)] = df.gene_id+'_'+df[mode]

    return df

def format_end_df(df):
    sort_cols = ['Chromosome', 'Start', 'End', 'Strand']
    df = df.sort_values(by=sort_cols)
    order = ['Chromosome', 'Start', 'End', 'Strand', 'Name', 'source']
    order = [o for o in order if o in df.columns]
    df = df[order]
    df.reset_index(drop=True, inplace=True)
    return df

def make_exon_df(n,c,e,s,g,t,nt=None,ag1=None,ag2=None):
    df = pd.DataFrame()
    df['Chromosome'] = c
    df['Start'] = [i[0] for i in e]
    df['End'] = [i[1] for i in e]
    df['Strand'] = s
    df['Feature'] = 'exon'
    cols = ['gene_name', 'gene_id']
    for c in cols:
        df[c] = g
    cols = ['transcript_id', 'transcript_name']
    for c in cols:
        df[c] = t
    if nt:
        df['new_transcript_id'] = nt
    if ag1:
        df['ag1'] = ag1
    if ag2:
        df['ag2'] = ag2

    # reorder exons and starts/ stops if needed
    df['new_Start'] = df[['Start', 'End']].min(axis=1)
    df['new_End'] = df[['Start', 'End']].max(axis=1)
    df.drop(['Start', 'End'], axis=1, inplace=True)
    df.rename({'new_Start':'Start',
               'new_End':'End'}, axis=1, inplace=True)
    s = s[0]
    if s == '+':
        ascending = True
    elif s == '-':
        ascending = False
    df.sort_values(by='Start', ascending=ascending, inplace=True)
    return df

def make_hier_entry(df, how='t'):
    """
    kind {'g','t'}
    """
    agg_dict = {'min_coord': 'min', 'max_coord': 'max'}
    t_df = df.copy(deep=True)
    t_df['min_coord'] = t_df[['Start', 'End']].min(axis=1)
    t_df['max_coord'] = t_df[['Start', 'End']].max(axis=1)
    if how == 't':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id', 'transcript_id', 'transcript_name',
                   'tss_id', 'tes_id',
                   'new_transcript_id', 'original_transcript_id',
                   'original_transcript_name', 'ag1', 'ag2']
        gb_cols = list(set(gb_cols)&(set(t_df.columns)))
    elif how == 'g':
        gb_cols = ['Chromosome', 'Strand', 'gene_name',
                   'gene_id']

    cols = gb_cols + ['min_coord', 'max_coord']
    t_df = t_df[cols]
    t_df = t_df.groupby(gb_cols).agg(agg_dict).reset_index()
    t_df.rename({'min_coord': 'Start', 'max_coord': 'End'}, axis=1, inplace=True)
    if how == 't':
        t_df['Feature'] = 'transcript'
    elif how == 'g':
        t_df['Feature'] = 'gene'

    return t_df

def make_test_gtf(ts):
    df = pd.concat(ts)
        # make transcript entries
    t_df = make_hier_entry(df, how='t')
    # make gene entries
    g_df = make_hier_entry(df, how='g')

    # concat everything and sort by gene id, transcript id, feature rank (gene =0, t =1, exon=2), then start coords
    df = pd.concat([df, t_df, g_df])
    df = sort_gtf(df)
    return df

# tests for agg_gtf
def test_agg_gtf(print_dfs=True):

    # need:
    # transcripts that don't need to be aggregated
    # transcripts that do need to be aggregated

    ts = []

    # t1 - transcript that won't need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[1,10], [14,20], [25,30]]
    s = ['+' for i in range(n)]
    g = 'g1'
    nt = 'g1[1,1,1]'
    t = 'g1_t1'
    ag1='known'
    ag2='p1'
    df = make_exon_df(n,c,e,s,g,t,nt,ag1,ag2)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t2 - similar transcript from g1 gene that won't need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[1,10], [14,20], [25,50]]
    s = ['+' for i in range(n)]
    g = 'g1'
    nt = 'g1[1,1,3]'
    t = 'g1_t2'
    ag1 = 'novel'
    ag2 = 'p1'
    df = make_exon_df(n,c,e,s,g,t,nt,ag1,ag2)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_3'
    ts.append(df)

    # t3 - transcript from a different gene that does need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[90,60], [45,30], [10,8]]
    s = ['-' for i in range(n)]
    g = 'g2'
    nt = 'g2[1,1,1]'
    t = 'g2_t1'
    ag1 = 'known'
    ag2 = 'p3'
    df = make_exon_df(n,c,e,s,g,t,nt,ag1,ag2)
    df['tss_id'] = 'g2_1'
    df['tes_id'] = 'g2_1'
    ts.append(df)

    # t4 - transcript that needs to be collapsed with t3
    n = 3
    c = ['1' for i in range(n)]
    e = [[90,60], [45,30], [10,8]]
    s = ['-' for i in range(n)]
    g = 'g2'
    nt = 'g2[1,1,1]'
    t = 'g2_t2'
    ag1 = 'novel'
    ag2 ='p4'
    df = make_exon_df(n,c,e,s,g,t,nt,ag1,ag2)
    df['tss_id'] = 'g2_1'
    df['tes_id'] = 'g2_1'
    ts.append(df)

    # make test gtf
    test_df = make_test_gtf(ts)
    test_df.rename({'transcript_id': 'original_transcript_id',
                    'transcript_name': 'original_transcript_name'},
                   axis=1, inplace=True)
    test_df.rename({'new_transcript_id': 'transcript_id'},
                   axis=1, inplace=True)
    test_df['transcript_name'] = test_df['transcript_id']

    # make ctrl df
    ts = []

    # t1 - transcript that won't need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[1,10], [14,20], [25,30]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1[1,1,1]'
    ag1='known'
    ag2='p1'
    df = make_exon_df(n,c,e,s,g,t,ag1=ag1,ag2=ag2)
    df['original_transcript_id'] = 'g1_t1'
    df['original_transcript_name'] = df['original_transcript_id']
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t2 - similar transcript from g1 gene that won't need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[1,10], [14,20], [25,50]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1[1,1,3]'
    ag1 = 'novel'
    ag2 = 'p1'
    df = make_exon_df(n,c,e,s,g,t,ag1=ag1,ag2=ag2)
    df['original_transcript_id'] = 'g1_t2'
    df['original_transcript_name'] = df['original_transcript_id']
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_3'
    ts.append(df)

    # t3 /t4 collapsed- transcript from a different gene that does need to be aggregated
    n = 3
    c = ['1' for i in range(n)]
    e = [[90,60], [45,30], [10,8]]
    s = ['-' for i in range(n)]
    g = 'g2'
    t = 'g2[1,1,1]'
    ag1 = 'known,novel'
    ag2 = 'p3,p4'
    df = make_exon_df(n,c,e,s,g,t,ag1=ag1,ag2=ag2)
    df['tss_id'] = 'g2_1'
    df['tes_id'] = 'g2_1'
    df['original_transcript_id'] = 'g2_t1,g2_t2'
    df['original_transcript_name'] = df['original_transcript_id']
    ts.append(df)


    ctrl = pr.PyRanges(make_test_gtf(ts)).df
    test = pr.PyRanges(agg_gtf(test_df)).df

    ctrl.reset_index(inplace=True, drop=True)
    test.reset_index(inplace=True, drop=True)

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

# tests for update_gtf_ends
def test_update_gtf_ends(print_dfs=True):

    ts = []
    # t1 - fwd strand transcript from gene w/ >2 transcripts
    n = 3
    c = ['1' for i in range(n)]
    e = [[1,10], [14,20], [25,30]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t1.5 - fwd strand transcript from gene that doesn't
    # need boundaries updated
    n = 3
    c = ['1' for i in range(n)]
    e = [[0, 10], [14,20], [25, 35]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1_t1.5'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t2 - rev. strand transcript from gene w/ >2 transcripts
    n = 3
    c = ['1' for i in range(n)]
    e = [[90,60], [45,30], [10,8]]
    s = ['-' for i in range(n)]
    g = 'g2'
    t = 'g2_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g2_1'
    df['tes_id'] = 'g2_1'
    ts.append(df)

    # t3 - rev. strand transcript from gene w/ >2 transcripts
    n = 3
    c = ['1' for i in range(n)]
    e = [[95,60], [45,30], [10,6]]
    s = ['-' for i in range(n)]
    g = 'g2'
    t = 'g2_t2'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g2_2'
    df['tes_id'] = 'g2_2'
    ts.append(df)

    # t4 - fwd strand gene w/ monoexonic transcript
    n = 1
    c = ['1' for i in range(n)]
    e = [[20,30]]
    s = ['+']
    g = 'g3'
    t = 'g3_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g3_1'
    df['tes_id'] = 'g3_1'
    ts.append(df)

    # t5 - rev strand gene w/ monoexonic transcript
    n = 1
    c = ['1' for i in range(n)]
    e = [[50, 40]]
    s = ['-']
    g = 'g4'
    t = 'g4_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g4_1'
    df['tes_id'] = 'g4_1'
    ts.append(df)

    # tss reference
    n = 5
    mode = 'tss'
    c = ['1' for i in range(n)]
    s = ['+', '-', '-', '+', '-']
    st = [0, 85, 91, 19, 50]
    e = [3, 93, 98, 22, 55]
    n = ['g1_1', 'g2_1', 'g2_2', 'g3_1', 'g4_1']
    source = 'v1'
    tss = make_end_df(c,s,st,e,n,source,mode)
    tss = pr.PyRanges(tss)

    # tes reference
    n = 5
    mode = 'tes'
    c = ['1' for i in range(n)]
    s = ['+', '-', '-', '+', '-']
    st = [25, 7, 4, 30, 35]
    e = [35, 9, 6, 31, 40]
    n = ['g1_1', 'g2_1', 'g2_2', 'g3_1', 'g4_1']
    source = 'v1'
    tes = make_end_df(c,s,st,e,n,source,mode)
    tes = pr.PyRanges(tes)

    test_df = make_test_gtf(ts)
    # test_df = pr.PyRanges(test_df)

    test = update_gtf_ends(test_df, tss, tes)

    # ctrl for update_gtf_ends
    # tests for update_gtf_ends
    ts = []
    # t1 - fwd strand transcript
    n = 3
    c = ['1' for i in range(n)]
    e = [[0,10], [14,20], [25,35]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t1.5 - fwd strand transcript from gene that doesn't
    # need boundaries updated
    n = 3
    c = ['1' for i in range(n)]
    e = [[0, 10], [14,20], [25, 35]]
    s = ['+' for i in range(n)]
    g = 'g1'
    t = 'g1_t1.5'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g1_1'
    df['tes_id'] = 'g1_1'
    ts.append(df)

    # t2 - rev. strand transcript from gene w/ >2 transcripts
    n = 3
    c = ['1' for i in range(n)]
    e = [[93,60], [45,30], [10,7]]
    s = ['-' for i in range(n)]
    g = 'g2'
    t = 'g2_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g2_1'
    df['tes_id'] = 'g2_1'
    ts.append(df)

    # t3 - rev. strand transcript from gene w/ >2 transcripts
    n = 3
    c = ['1' for i in range(n)]
    e = [[98,60], [45,30], [10,4]]
    s = ['-' for i in range(n)]
    g = 'g2'
    t = 'g2_t2'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g2_2'
    df['tes_id'] = 'g2_2'
    ts.append(df)

    # t4 - fwd strand gene w/ monoexonic transcript
    n = 1
    c = ['1' for i in range(n)]
    e = [[19,31]]
    s = ['+']
    g = 'g3'
    t = 'g3_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g3_1'
    df['tes_id'] = 'g3_1'
    ts.append(df)

    # t4 - rev strand gene w/ monoexonic transcript
    n = 1
    c = ['1' for i in range(n)]
    e = [[55, 35]]
    s = ['-']
    g = 'g4'
    t = 'g4_t1'
    df = make_exon_df(n,c,e,s,g,t)
    df['tss_id'] = 'g4_1'
    df['tes_id'] = 'g4_1'
    ts.append(df)

    ctrl = make_test_gtf(ts)
    ctrl.reset_index(inplace=True, drop=True)
    test.reset_index(inplace=True, drop=True)

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

import pyranges as pr
import pandas as pd
import pdb
import h5py
import numpy as np
import os

##### helper functions #####

def parse_file_input(input, ext):
    """
    From either a file w/ filenames on each line or a comma-separated
    list of files, get an ordered list of filenames

    Parameters:
        input (str): Either a path to a file where each line contains
            a file path or a comma separated ordered list of filenames
        ext (str): Extension of the anticipated files

    Returns:
        input (list of str): Ordered list of filepaths
    """

    # input files given as list on cmd line
    if input.endswith('.{}'.format(ext)):
        input = input.split(',')

    # bed files given in a file
    else:
        df = pd.read_csv(input, header=None, names=['file'])
        input = df.file.tolist()

    return input

def change_all_dtypes(df, dtype):
    """
    Change all dtypes in a dataframe to a specified dtype

    Parameters:
        df (pandas DataFrame): DataFrame to alter
        dtype (type): Variable type to convert each column to

    Returns:
        df (pandas DataFrame): DataFrame with all dtypes changed
    """
    for c in df.columns:
        df[c] = df[c].astype(dtype)
    return df

def get_nov_ranks():
    """
    Get rank for each novelty category

    Returns:
        nov_rank (dict): Dictionary of novelty category: int rank
        rank_nov (dict): Dictionary of int rank: novelty category
    """
    novs = ['Known', 'NIC', 'ISM_rescue',
            'ISM', 'NNC', 'Genomic', 'Antisense',
            'Intergenic']
    ranks = [i for i in range(len(novs))]
    nov_rank = dict([(nov, rank) for nov, rank in zip(novs, ranks)])
    rank_nov = dict()
    for key, item in nov_rank.items():
        rank_nov[item] = key

    return nov_rank, rank_nov

def get_ism_ranks():
    """
    Get rank for each ISM subtype category

    Returns:
        ism_rank (dict): Dictionary of ism subtype: int rank
        rank_ism (dict): Dictionary of int rank: ism subtype
    """
    ism_rank = {'Both': 0,
                'Prefix': 1,
                'Suffix': 2,
                None: 3}
    rank_ism = dict()
    for key, item in ism_rank.items():
        rank_ism[item] = key

    return ism_rank, rank_ism

def get_non_dataset_cols():
    """
    Get a list of ordered columns from TALON abundance file that
    are the non-dataset columns

    Returns:
        non_dataset_cols (list of str): List of non-dataset columns
    """
    non_dataset_cols = ['gene_ID', 'transcript_ID', 'annot_gene_id',
                   'annot_transcript_id', 'annot_gene_name',
                   'annot_transcript_name', 'n_exons', 'length',
                   'gene_novelty', 'transcript_novelty', 'ISM_subtype']
    return non_dataset_cols

def get_dataset_cols(df):
    """
    Get a list of the dataset columns from a TALON abundance file

    Parameters:
        df (pandas DataFrame): TALON abundance file

    Returns:
        dataset_cols (list of str): List of dataset columns
    """
    non_dataset_cols = get_non_dataset_cols()
    dataset_cols = [ x for x in list(df.columns) \
                        if x not in non_dataset_cols ]
    return dataset_cols

def get_ic(gtf_pr):
    """
    Get a hyphen-separated representation of each transcript's intron chain
    from a PyRanges GTF

    Parameters:
        gtf_pr (pyranges PyRanges): GTF PyRanges object

    Returns:
        df (pandas DataFrame): DataFrame detailing intron chain, gene, strand,
            chromosome, and transcript that intron chain was seen in
    """
    df = gtf_pr.df.copy(deep=True)

    # restrict to exon entries
    df = df.loc[df.Feature == 'exon']
    cols = ['Chromosome', 'Strand', 'Start', 'End',
            'transcript_id', 'gene_id']
    df = df[cols]

    # melt to isolate individual coordinates
    df = pd.melt(df, id_vars=['Chromosome', 'Strand', 'transcript_id', 'gene_id'],
                value_vars=['Start', 'End'],
                value_name='Coord')
    df.drop('variable', axis=1, inplace=True)

    # sort to order coordinates correctly
    df.Coord = df.Coord.astype(int)
    fwd = df.loc[df.Strand == '+'].copy(deep=True)
    rev = df.loc[df.Strand == '-'].copy(deep=True)

    fwd.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, True], inplace=True)
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False], inplace=True)
    df = pd.concat([fwd, rev])

    # create intron chain strings
    df.Coord = df.Coord.astype(str)
    df = df.groupby(['Chromosome', 'Strand',
                     'transcript_id', 'gene_id'], observed=True)['Coord'].apply('-'.join).reset_index()

    # remove tss and tes from intron chain
    df['temp'] = df.Coord.str.split('-', n=1, expand=True)[1]
    df['ic'] = df.temp.str.rsplit('-', n=1, expand=True)[0]
    df.loc[~df.ic.str.contains('-'), 'ic'] = '' # for monoexonic transcripts
    df.drop(['temp', 'Coord'], axis=1, inplace=True)

    return df

def number_tss_ic_tes(df, mode):
    """
    Assign a number to each tss, intron chain, or tes per gene based on the
    status of the transcript

    Parameters:
        df (pandas DataFrame): DataFrame derived from PyRanges
            with gene id, tss/ic/tes id, and tags to order transcripts
            with
        mode (str): {'tss', 'ic', 'tes'}

    Returns:
        df (pandas DataFrame): DataFrame with an identifier for each unique
            tss/ic/tes per gene
    """
    # groupby feature but record which feature
    # each transcript id uses
    gb_cols = ['Chromosome', 'Strand', mode, 'gene_id']
    subset_cols = ['transcript_id', 'Chromosome', 'Strand', mode,
                   'basic_set', 'MANE_Select', 'appris_principal', 'gene_id']

    # for ics that are monoexonic
    df.loc[df[mode].isnull(), mode] = ''

    if mode == 'tss' or mode == 'tes':
        gb_cols += ['Start', 'End']
        subset_cols += ['Start', 'End']

    sort_cols = ['gene_id', 'MANE_Select',
                 'appris_principal', 'basic_set']

    # agg drops nan ics!
    df = df[subset_cols].groupby(gb_cols,
                           observed=True).agg({'transcript_id': ','.join,
                                     'MANE_Select': 'max',
                                     'basic_set': 'max',
                                     'appris_principal': 'min'}).reset_index()

    # compute feature number based on tags
    df['{}_num'.format(mode)] = df.sort_values(by=sort_cols,
                                 ascending=[True, False, True, False],
                                 na_position='last')\
                                 .groupby(['gene_id'])\
                                 .cumcount() + 1

    return df

def number_gtf_ends(bed, gtf, mode):
    """
    As a part of calling ends from a gtf, use the tags and gene id
    in the gtf to assign each region a name / number

    Parameters:
        bed (pyranges PyRanges): Output of get_ends_from_gtf
        gtf (str): Path to gtf file
        mode (str): {'tss', 'tes'}

    Returns:
        bed (pyranges PyRanges): Bed file with Name field
            consisting of stable gene id_end region number
    """

    # read gtf with extra tag info
    gr_gtf = get_transcript_ref(gtf)
    gr_gtf = add_stable_gid(gr_gtf)
    gr_gtf = pr.PyRanges(gr_gtf)

    # add a temporary unique identifier for each end
    bed = bed.df
    bed[mode] = [i for i in range(len(bed.index))]
    bed = pr.PyRanges(bed)

    # get ends from reference gtf and join with called ends
    cols = ['Chromosome', 'Start', 'End', 'Strand',
             'gene_id', 'transcript_id', 'MANE_Select',
             'basic_set', 'appris_principal']
    if mode == 'tss':
        temp = pr.PyRanges(gr_gtf.features.tss().df[cols])
    elif mode == 'tes':
        temp = pr.PyRanges(gr_gtf.features.tes().df[cols])
    bed = bed.join(temp,
               strandedness='same',
               how='left',
               slack=0)
    bed = bed.df
    bed = bed.loc[bed.gene_id == bed.gene_id_b]
    cols = ['Start_b', 'End_b', 'Strand_b', 'gene_id_b']
    bed.drop(cols, axis=1, inplace=True)

    # use the tags from the gtf to assign an actual end id
    bed = number_tss_ic_tes(bed, mode=mode)
    # bed['gene_id'] = bed.gene_id.str.split(pat='.', n=1, expand=True)[0]
    c = '{}_num'.format(mode)
    bed['Name'] = bed.gene_id+'_'+bed[c].astype(str)
    cols = [mode, 'transcript_id', 'MANE_Select',
            'basic_set', 'appris_principal', c, 'gene_id']
    bed.drop(cols, axis=1, inplace=True)
    bed = pr.PyRanges(bed)

    return bed

##### for aggregating ends and intron chains #####

def get_gene_feat_max(df, mode):
    """
    Get the maximum current number of a feature in an existing df

    Parameters:
        df (pandas DataFrame): DataFrame with all the entries from the old
            DataFrame
        mode (str): {ic, tss, tes}

    Returns:
        temp (pandas DataFrame): DataFrame detailing the max value for
            each gene
    """

    max_c = '{}_max'.format(mode)

    df[mode] = df[mode].astype(int)
    temp = df[['gene_id', mode]].copy(deep=True)
    temp = temp.groupby('gene_id').max().reset_index()
    temp.rename({mode: max_c}, axis=1, inplace=True)

    return temp

def renumber_new_feats(df, g_maxes, mode):
    """
    Update the id numbers of each novel feature that's
    being added

    Parameters:
        df (pandas DataFrame): DataFrame with all the novel entries
            being added to the reference
        g_maxes (pandas DataFrames): Output from `get_gene_feat_max`
        mode (str): {ic, tss, tes}

    Returns:
        df (pandas DataFrame): DataFrame with the new numbers for
            each tss, tes, or ic
    """
    max_c = '{}_max'.format(mode)
    new_c = '{}_new'.format(mode)

    # merge with new ends that are being added
    df = df.merge(g_maxes, how='left', on='gene_id')
    df[max_c].fillna(0, inplace=True)
    sort_cols = ['gene_id', mode]

    # renumber new features
    df[new_c] = df.sort_values(by=sort_cols,
                                ascending=[True, True])\
                                .groupby(['gene_id'])\
                                .cumcount()+1
    df[new_c] = df[new_c].astype(int) + df[max_c].astype(int)

    return df

def format_agg_2_ends_bed(bed, mode):
    """
    Format bed file for agg_2_ends

    Parameters:
        bed (pyranges PyRanges): DataFrame of bed file
        mode (str): {'tss', 'tes'}

    Returns:
        bed (pandas DataFrame): DataFrame with added gene id
            and strand info
        gid (bool): Whether or not the bed file contained a gene id
        strand (bool): Whether or not the bed file contained a strand
    """
    cols = ['Name', 'gene_id', mode, 'Strand']

    gid = True
    strand = True
    for col in cols:
        if col not in bed.df.columns:
            bed = bed.df
            bed[col] = [np.nan for i in range(len(bed.index))]

            # record some bools
            if col == 'gene_id':
                gid = False
            elif col == 'Strand':
                strand = False

            bed[col] = bed[col].astype(str)
            bed = pr.PyRanges(bed)

    return bed, gid, strand

def get_source_map_fname(o):
    """
    Generate file name to match output bed file from `agg_ends`.

    Parameters:
        o (str): File name

    Returns:
        o (str): File name using suffix '_source_map.bed'
    """
    return ''.join(os.path.splitext(o)[:-1])+'_source_map.bed'

def agg_2_ends(bed1, bed2,
               strand,
               gid,
               slack,
               add_ends,
               mode):
    """
    Parameters:
        bed1 (pyranges PyRanges): Bed PyRanges object for existing ends
        bed2 (pyranges PyRanges): Bed PyRanges object for new ends
        strand (bool): Whether bed2 has strand info
        gid (bool): Whether bed2 has gene id info
        slack (int): Maximum allowable distance between ends in bed1 and bed2
            to call them the same end
        add_ends (bool): Whether to initialize new regions from bed2
        mode (str): {'tss', 'tes'}
    """
    source1 = bed1.df.source.unique().tolist()[0]
    source2 = bed2.df.source.unique().tolist()[0]

    new_c = '{}_new'.format(mode)
    max_c = '{}_max'.format(mode)

    # convert into int64
    bed1 = pr.PyRanges(bed1.df, int64=True)
    bed2 = pr.PyRanges(bed2.df, int64=True)

    # depending on whether the new bed has strand information,
    # construct the join call
    if strand:
        temp_joined = bed1.join(bed2,
            strandedness='same',
            suffix='_new',
            slack=slack,
            how='left')
    elif not strand:
        bed2 = bed2.df
        bed2.drop('Strand', axis=1, inplace=True)
        bed2 = pr.PyRanges(bed2, int64=True)
        temp_joined = bed1.join(bed2,
            strandedness=False,
            suffix='_new',
            slack=slack,
            how='left')

    # format null starts as actual nans b/c of join
    temp_joined = temp_joined.df
    temp_joined.sort_values(by=['Chromosome', 'Start', 'End'], inplace=True)
    temp_joined.loc[temp_joined.Start_new == -1, 'Start_new'] = np.nan

    # df to hold final end annotations
    df = pd.DataFrame()

    ### old ends ###

    # situation 1: ends match across the datasets in coord and gene id
    if gid:
        temp = temp_joined.loc[temp_joined.gene_id == temp_joined.gene_id_new].copy(deep=True)
    else:
        temp = temp_joined.loc[~temp_joined.Start_new.isnull()].copy(deep=True)
    temp.source = temp.source+','+temp.source_new
    df = pd.concat([df, temp])

    # create source map for ends from this source
    m_source = df.loc[df.source_new == source2]
    m_source = m_source[['Chromosome', 'Start_new', 'End_new', 'Strand',
                   'source_new', 'Name', 'id_new']].copy(deep=True)
    m_source.rename({'Start_new': 'Start',
                     'End_new': 'End',
                     'source_new': 'source',
                     'id_new': 'id'},
                     axis=1, inplace=True)

    # situation 2: ends are only in the first dataset
    if gid:
        # end either didn't match something or matched wrong gid
        temp = temp_joined.loc[(temp_joined.Start_new.isnull())|\
                               (temp_joined.gene_id!=temp_joined.gene_id_new)].copy(deep=True)

        # restrict to ends that haven't been added to our main df
        temp = temp.loc[~temp.Name.isin(df.Name.tolist())]

    else:
        temp = temp_joined.loc[temp_joined.Start_new.isnull()].copy(deep=True)

    # mark these ends with a -1 id so we don't count them as already
    # processed
    temp['id_new'] = -1
    df = pd.concat([df, temp])

    # restrict to relevant columns
    cols = ['Chromosome', 'Start', 'End', 'Strand',
            'Name', 'gene_id', 'source', mode, 'id_new']
    df = df[cols]
    df.rename({'id_new': 'id'}, axis=1, inplace=True)

    ### new ends, only add if we're allowing them to be independent
    ### end support
    if add_ends and strand and gid:

        new_df = pd.DataFrame()

        drop_cols = ['Start', 'End', 'Strand', 'gene_id', 'source', 'Name', 'id', mode]
        m = {'Start_new': 'Start',
             'End_new': 'End',
             'gene_id_new': 'gene_id',
             'Strand_new': 'Strand',
             'source_new': 'source',
             'Name_new': 'Name',
             'id_new': 'id',
             new_c: mode}

        # situation 3: the ends overlapped, but the gene ids didn't match AND
        # the end hasn't been added yet
        # temp = temp_joined.loc[(temp_joined.gene_id!=temp_joined.gene_id_new)&\
        #                        (temp_joined.gene_id_new!='-1')&\
        #                        (~temp_joined.id_new.isin(df.id.tolist()))].copy(deep=True)
        temp = temp_joined.loc[(temp_joined.gene_id!=temp_joined.gene_id_new)&(temp_joined.gene_id_new!='-1')&(~temp_joined.id_new.isin(df.id.tolist()))].copy(deep=True)
        temp.drop(drop_cols, axis=1, inplace=True)
        temp.rename(m, axis=1, inplace=True)
        temp.drop_duplicates(inplace=True)

        new_df = pd.concat([new_df, temp])

        # situation 4: the ends are brand new and didn't overlap at all in the existing ends
        bed2 = bed2.df
        inds = list(set(bed2.id.tolist())-set(df.id.tolist())-set(new_df.id.tolist()))
        temp = bed2.loc[inds]
        new_df = pd.concat([new_df, temp])

        g_maxes = get_gene_feat_max(df, mode)
        new_df = renumber_new_feats(new_df, g_maxes, mode)

        # some df formatting
        new_df.drop([mode, max_c], axis=1, inplace=True)
        new_df.rename({new_c: mode}, axis=1, inplace=True)

        # finally, concatenate new and old df
        df = pd.concat([df, new_df])

    # drop duplicates that could have arisen from entries
    # on multiple strands or from multiple subregions in bed2
    df['n_sources'] = df.source.str.count(pat=',')
    df.sort_values(by='n_sources', ascending=False)
    df = df.loc[~df[['Chromosome', 'Start', 'End',
                    'Strand', 'gene_id']].duplicated(keep='last')]
    df.drop('n_sources', axis=1, inplace=True)

    # drop unnecessary columns and create new names
    # and do some extra formatting
    df.drop('id', axis=1, inplace=True)
    df['Name'] = df.gene_id+'_'+df[mode].astype(str)
    df[mode] = df[mode].astype(int)
    df['Start'] = df.Start.astype(int)
    keep_cols = ['Chromosome', 'Start', 'End', 'Strand', 'Name', 'gene_id', mode, 'source']
    df = df[keep_cols]

    df = df.sort_values(by=['Chromosome', 'Start'])
    df = df.sort_values(by='Name')

    df['id'] = [i for i in range(len(df.index))]

    # add new ends to source map or ends that weren't added
    if add_ends:
        temp = df.loc[df.source == source2].copy(deep=True)
        temp = temp[['Chromosome', 'Start', 'End', 'Strand', 'source', 'Name']]
        m_source = pd.concat([m_source, temp])
    else:
        merge_cols = ['Chromosome', 'Start', 'End', 'id']
        if strand:
            merge_cols.append('Strand')
        if gid:
            merge_cols.append('gene_id')
            # m_source = split_cerberus_id(m_source, mode)
        temp = bed2.df.copy(deep=True)
        temp = temp[merge_cols]
        temp = temp.loc[~temp.id.isin(m_source.id.tolist())]
        # m_source = m_source.merge(temp, how='outer', on=merge_cols)
        m_source = pd.concat([m_source, temp])
        m_source['source'] = source2
        # print(len(m_source[merge_cols].drop_duplicates()))

    # update dtypes and final formatting
    m_source.drop('id', axis=1, inplace=True)
    m_source['Start'] = m_source.Start.astype('int')
    m_source['End'] = m_source.End.astype('int')

    return df, m_source

def agg_2_ics(ic1, ic2):
    """
    Aggregate and reconcile intron chains from two
    intron chain dataframes

    Parameters:
        ic1, ic2 (pandas DataFrame): DataFrame with the old and
            new intron chains respecitvely

    Returns:
        df (pandas DataFrame): DataFrame with merged and
            renamed intron chains
    """

    mode = 'ic'
    max_c = '{}_max'.format(mode)
    new_c = '{}_new'.format(mode)

    # if we have more than 1 set of ics, merge on chrom, strand,
    # ic coords, and gene id
    ic1 = ic1.merge(ic2,
                  on=['Chromosome', 'Strand', 'Coordinates', 'gene_id'],
                  how='outer', suffixes=('', '_new'))
    # https://stackoverflow.com/questions/62681371/python-combining-names-with-missing-values/62681510#62681510
    ic1['source'] = ic1[['source', 'source_new']].stack().groupby(level=0).agg(','.join)

    # get new ic numbers for duplicate entries
    old = ic1.loc[~ic1.ic.isnull()].copy(deep=True)
    g_maxes = get_gene_feat_max(old, 'ic')
    new = ic1.loc[ic1.ic.isnull()].copy(deep=True)
    new = renumber_new_feats(new, g_maxes, 'ic')

    # some df formatting
    new.drop([mode, max_c, 'Name_new'], axis=1, inplace=True)
    new.rename({new_c: mode}, axis=1, inplace=True)

    # finally, concatenate new and old df
    df = pd.concat([old, new])

    # drop unnecessary columns and create new names
    # and do some extra formatting
    df['Name'] = df.gene_id+'_'+df[mode].astype(str)
    df.drop(['Name_new', new_c, 'source_new'], axis=1, inplace=True)

    return df

##### helpers for convert_transcriptome #####
def merge_ics(df, ic):
    """
    Assign each transcript in df to an intron chain in ic

    Parameters:
        df (pandas DataFrame): DataFrame w/ intron chain for each transcript
        ic (pandas DataFrame): DataFrame from cerberus reference with intron chains

    Returns:
        df (pandas DataFrame): DataFrame indicating which intron chain is used
            for each transcript
    """
    # merge on intron chain, strand, chromosome, and gene id
    df = df.merge(ic, how='left',
                  on=['Chromosome', 'Strand',
                      'Coordinates', 'gene_id'])

    # formatting
    df.rename({'Name': 'ic_id'}, axis=1, inplace=True)
    df = df[['transcript_id', 'ic', 'ic_id']]

    return df

def merge_ends(ends, ref, mode):
    """
    Merge ends from a GTF with those already annotated in a cerberus reference

    Parameters:
        ends (pyranges PyRanges): PyRanges object of ends from GTF annotation
        ref (pyranges PyRanges): PyRanges object of reference ends from cerberus
            annotation
        mode (str): {'tss', 'tes'}

    Returns:
        df (pandas DataFrame): DataFrame detailing which transcript uses which
            end from the cerberus reference
    """

    # whether we should look upstream (tss) or downstream (tes)
    if mode == 'tss':
        direction = 'upstream'
    elif mode == 'tes':
        direction = 'downstream'

    # limit to relevant columns
    ends = ends[['Chromosome', 'Start', 'End', 'Strand',
                 'gene_id', 'transcript_id']]

    # get only the relevant columns and deduplicate
    ends = ends.df
    t_ends = ends.copy(deep=True)
    ends.drop('transcript_id', axis=1, inplace=True)
    ends.drop_duplicates(inplace=True)
    ends = pr.PyRanges(ends)

    # find closest interval in ref
    ends = ends.nearest(ref,
                        strandedness=None,
                        how=direction)
    # pdb.set_trace()

    # fix the ends with mismatching gene ids - this part can be slow :(
    ends = ends.df
    fix_ends = ends.loc[ends.gene_id != ends.gene_id_b]
    fix_ends = fix_ends[['Chromosome', 'Start', 'End', 'Strand',
                         'gene_id']]
    ends = ends.loc[ends.gene_id == ends.gene_id_b]

    for i, gid in enumerate(fix_ends.gene_id.unique().tolist()):
        gene_ends = fix_ends.loc[fix_ends.gene_id == gid].copy(deep=True)
        gene_ends = pr.PyRanges(gene_ends)
        gene_refs = ref.df.loc[ref.df.gene_id == gid].copy(deep=True)
        gene_refs = pr.PyRanges(gene_refs)
        gene_ends = gene_ends.nearest(gene_refs,
                                      strandedness=None,
                                      how='direction')
        gene_ends = gene_ends.df
        ends = pd.concat([ends, gene_ends])

        if i % 100 == 0:
            print('Processed {} / {} genes'.format(i, len(fix_ends.gene_id.unique().tolist())))
    # pdb.set_trace()

    # merge back in to get transcript ids
    t_ends = t_ends.merge(ends, how='left',
                      on=['Chromosome', 'Start', 'End', 'Strand', 'gene_id'])

    # formatting
    t_ends.rename({'Name': '{}_id'.format(mode)}, axis=1, inplace=True)
    t_ends = t_ends[['transcript_id', '{}_id'.format(mode), mode]]

    return t_ends

##### helpers for replace_ab_ids #####
def map_transcripts(df, m_df, tname, tid):
    """
    Parameters:
        df (pandas DataFrame): DataFrame with transcript ids
        m_df (pandas DataFrame): Map DataFrame from cerberus reference
        tname (str): Column name for transcript name in df
        tid (str): Column name for transcript id in df

    Returns:
        df (pandas DataFrame): TALON abundance file with tids swapped out
    """

    # fix transcript ids in abundance file
    t_map = m_df[['original_transcript_id', 'original_transcript_name',
                  'transcript_name', 'transcript_id']]
    df = df.merge(t_map, left_on=tid, right_on='original_transcript_id')

    df.drop(['original_transcript_id', 'original_transcript_name',
             tid, tname],
            axis=1, inplace=True)
    df.rename({'transcript_name': tname,
               'transcript_id': tid},
              axis=1, inplace=True)

    return df

def agg_ab(df):
    """
    Aggregates expression over transcripts that now are the same based on
    assignment to cerberus reference ends

    Parameters:
        df (pandas DataFrame): TALON abundance file output from `map_transcripts`

    Returns:
        df (pandas DataFrame): Abundance file with counts and other metadata
            per transcript aggregated for duplicated transcripts
    """
    gb_cols = ['gene_ID', 'annot_gene_id', 'annot_gene_name',
               'gene_novelty', 'annot_transcript_name',
               'annot_transcript_id']

    # handle properties which won't always correspond across the transcripts
    # these are all subject to change
    nov_rank, rank_nov = get_nov_ranks()
    df['nov_rank'] = df.transcript_novelty.map(nov_rank)
    df.drop('transcript_novelty', axis=1, inplace=True)

    ism_rank, rank_ism = get_ism_ranks()
    df['ism_rank'] = df.ISM_subtype.map(ism_rank)
    df.drop('ISM_subtype', axis=1, inplace=True)

    df['transcript_ID'] = df.transcript_ID.astype(str)
    agg_dict = {'transcript_ID': ','.join,
                'n_exons': 'mean',
                'length': 'mean',
                'nov_rank': 'min',
                'ism_rank': 'min'}
    cols = gb_cols + list(agg_dict.keys())
    for c in list(set(df.columns)-set(cols)):
        agg_dict[c] = 'sum'

    df = df.groupby(gb_cols).agg(agg_dict).reset_index()
    df['transcript_novelty'] = df.nov_rank.map(rank_nov)
    df.drop(['nov_rank'], axis=1, inplace=True)
    df['ISM_subtype'] = df.ism_rank.map(rank_ism)
    df.drop(['ism_rank'], axis=1, inplace=True)

    return df

##### helpers for replace_gtf_ids #####
def get_stranded_gtf_dfs(df):
    """
    Split a GTF df into fwd and rev strands

    Parameters:
        df (pandas DataFrame): DF of gtf

    Returns:
        fwd (pandas DataFrame): DF of all forward-stranded entries from GTF
        rev (pandas DataFrame): DF of all reverse-stranded entries from GTF
    """
    rev = df.loc[df.Strand == '-'].copy(deep=True)
    fwd = df.loc[df.Strand == '+'].copy(deep=True)

    return fwd, rev

def sort_gtf(df):
    """
    Sort a GTF into its proper ordering

    Parameters:
        df (pandas DataFrame): DF of GTF

    Returns:
        df (pandas DataFrame): DF of GTF, sorted
    """
    df['feature_rank'] = df.Feature.map({'gene':0, 'transcript':1, 'exon':2})
    df.feature_rank = df.feature_rank.astype(int)

    fwd, rev = get_stranded_gtf_dfs(df)

    df = pd.DataFrame()
    for temp in [fwd, rev]:
        if len(temp.index) > 0:
            strand = temp.Strand.values.tolist()[0]
            if strand == '+':
                ascending = True
            elif strand == '-':
                ascending = False
            temp.sort_values(by=['gene_id', 'transcript_id', 'feature_rank', 'Start'],
                             ascending=[True, True, True, ascending],
                             na_position='first', inplace=True)

            df = pd.concat([df, temp], ignore_index=True)
    df.drop('feature_rank', axis=1, inplace=True)
    return df

def get_update_ends_settings(strand, mode):
    """
    Returns which columns to refer to and which min/max function
    to use depending on looking at forward / rev strand or
    tss / tes

    Parameters:
        strand (str): {'+', '-'}
        mode (str): {'tss', 'tes'}

    Returns:
        old_end (str): Name of column to modify; {'Start', 'End'}
        new_end (str): Name of column to pull new value from; {'Start_end', 'End_end'}
        gene_func (str): What function to apply to new_end; {'min', 'max'}
    """
    if mode == 'tss':
        if strand == '+':
            old_end = 'Start'
            new_end = 'Start_end'
            gene_func = 'min'
        elif strand == '-':
            old_end = 'End'
            new_end = 'End_end'
            gene_func = 'max'
    elif mode == 'tes':
        if strand == '+':
            old_end = 'End'
            new_end = 'End_end'
            gene_func = 'max'
        elif strand == '-':
            old_end = 'Start'
            new_end  = 'Start_end'
            gene_func = 'min'

    return old_end, new_end, gene_func

def update_transcript_ends(df, mode, strand):
    """
    Update the ends of transcripts and the first / last exon
    in a GTF. GTF must be sorted!

    Parameters:
        df (pandas DataFrame): Sorted DF of GTF with 'Start_end', and 'End_end'
            columns denoting the boundaries of each end region
        mode (str): {'tss', 'tes'}
        strand (str): {'+', '-'}

    Returns:
        df (pandas DataFrame): DF of GTF with transcript and
            exon ends modified
    """
    old_col, new_col, gene_func = get_update_ends_settings(strand, mode)

    temp = df[['Feature', 'gene_id', 'original_transcript_id', 'Strand', 'Start', 'End', 'Start_end', 'End_end']].copy(deep=True)
    temp = temp.loc[temp.Feature != 'gene']
    if mode == 'tss':
        inds = temp.groupby('original_transcript_id').head(2).index.tolist()
    elif mode == 'tes':
        inds = temp.groupby('original_transcript_id').head(1).index.tolist()
        inds += temp.groupby('original_transcript_id').tail(1).index.tolist()

    df.loc[inds, old_col] = df.loc[inds, new_col]

    # convert float dtypes
    df.Start = df.Start.astype(int)
    df.End = df.End.astype(int)

    return df

def update_gene_ends(df, mode, strand):
    """
    Update the ends of genes in a GTF.

    Parameters:
        df (pandas DataFrame): GTF dataframe
        mode (str): {'tss', 'tes'}
        strand (str): {'+', '-'}

    Returns:
        df (pandas DataFrame): DataFrame of GTF with gene ends updated
    """
    # determine which ends we're updating and how we're doing so
    old_col, new_col, gene_func = get_update_ends_settings(strand, mode)

    # get min or max of transcript ends depending on settings
    temp = df[['Feature', 'gene_id', old_col]].copy(deep=True)
    temp = temp.loc[temp.Feature == 'transcript']
    temp = temp.groupby(['gene_id', 'Feature'], observed=True).agg(gene_func).reset_index()
    temp.drop('Feature', axis=1, inplace=True)

    # add that coord to gene end
    df = df.merge(temp, on='gene_id', suffixes=('', '_gene'))
    inds = df.loc[df.Feature == 'gene'].index.tolist()
    df.loc[inds, old_col] = df.loc[inds, '{}_gene'.format(old_col)]
    df.drop('{}_gene'.format(old_col), axis=1, inplace=True)

    return df

def update_gtf_ends(gtf, tss, tes):
    """
    Update gene, transcript, and exon boundaries to be
    furthest upstream or downstream entry for end used

    Parameters:
        gtf (pandas DataFrame): DF of GTF
        tss (pyranges PyRanges): PyRanges object of reference TSSs
        tes (pyranges PyRanges): PyRanges object of reference TESs

    Returns:
        gtf (pandas DataFrame): DF of GTF with updated ends
            based on the TSSs and TESs used in the input beds
    """
    gtf = gtf.copy(deep=True)

    for mode, ends in zip(['tss', 'tes'], [tss, tes]):
        ends = ends.df
        ends = ends[['Start', 'End', '{}_id'.format(mode)]]
        gtf = gtf.merge(ends, how='left',
                        on='{}_id'.format(mode),
                        suffixes=('', '_end'))

        fwd, rev = get_stranded_gtf_dfs(gtf)
        df = pd.DataFrame()
        for strand, temp in zip(['+', '-'], [fwd, rev]):

            # fix exon, transcript, and gene boundaries
            temp = update_transcript_ends(temp, mode, strand)
            temp = update_gene_ends(temp, mode, strand)
            df = pd.concat([df, temp], ignore_index=True)

        gtf = df.copy(deep=True)
        gtf.drop(['Start_end', 'End_end'], axis=1, inplace=True)

    return gtf

def agg_gtf(df):
    """
    Deduplicate GTF transcripts that have the same triplet id

    Parameters:
        df (pandas DataFrame): DF of gtf from `update_ends`

    Returns:
        df (pandas DataFrame): DF of gtf with deduplicated
            transcript / exon entries based on the triplet id
    """

    def collapse_non_gb_col(x):
        x = x.fillna('')
        x = x.astype(str)
        x = x.unique().tolist()
        x = ','.join(x)
        if x == '':
            x = np.nan
        return(x)

    gb_cols = ['Chromosome',
                 'Feature',
                 'Start', 'End',
                 'Score', 'Strand', 'Frame', 'gene_id', 'gene_name',
                 'gene_status', 'gene_type', 'talon_gene',
                 'ic', 'ic_id', 'tss_id', 'tss', 'tes_id', 'tes', 'transcript_id',
                 'transcript_name']
    gb_cols = list(set(df.columns)&set(gb_cols))

    # get dictionary:function mapping for aggregation
    agg_cols = list(set(df.columns)-set(gb_cols))
    agg_dict = dict()
    for c in agg_cols:
        agg_dict[c] = collapse_non_gb_col

    # get collapsed features to add to deduplicated df
    t_df = df.loc[df.Feature == 'transcript'].copy(deep=True)
    non_gb_cols = list(set(t_df.columns.tolist())-set(gb_cols))
    t_df[non_gb_cols].fillna('', inplace=True)
    t_df = t_df.groupby(gb_cols, observed=True).agg(agg_dict).reset_index()
    t_df = t_df[['transcript_id']+non_gb_cols]
    collapsed_feats = t_df.copy(deep=True)

    # deduplicate df based only on transcript id
    temp = df[['transcript_id', 'original_transcript_id']].drop_duplicates()
    dupe_old_tids = temp.loc[temp.transcript_id.duplicated(keep='first'), 'original_transcript_id']
    df = df.loc[~df.original_transcript_id.isin(dupe_old_tids)]

    # replace the non gb columns with the ones that we already grouped
    df.drop(non_gb_cols, axis=1, inplace=True)
    df = df.merge(collapsed_feats, how='left', on='transcript_id')

    return df

##### writers #####

def write_h5(ic, tss, tes, oname,
             tss_map=None,
             tes_map=None,
             m=None):
    """
    Write a cerberus transcriptome as an h5df file

    Parameters:
        ic (pandas DataFrame): DataFrame of intron chains
        tss (pandas DataFrame): DataFrame in bed format of tsss
        tes (pandas DataFrame): DataFrame in bed format of tess
        oname (str): Output file name ending in '.h5'
        tss_map (pandas DataFrame): DataFrame in bed format of map from
            each input tss to which cerberus tss it was assigned to
        tes_map (pandas DataFrame): DataFrame in bed format of map from
            each input tes to which cerberus tes it was assigned to
        m (pandas DataFrame): DataFrame of map file
    """

    def cat_to_obj(df):
        temp = df.dtypes.to_frame()
        cols = temp.loc[temp[0] == 'category'].index.tolist()
        df[cols] = df[cols].astype('object')
        return df

    ic.to_hdf(oname, 'ic', mode='w')
    tss.to_hdf(oname, 'tss', mode='a', format='table')
    tes.to_hdf(oname, 'tes', mode='a', format='table')

    if not isinstance(tss_map, pd.DataFrame):
        tss_map = pd.DataFrame()
    tss_map = cat_to_obj(tss_map)
    tss_map.to_hdf(oname, 'tss_map', mode='a', format='table')

    if not isinstance(tes_map, pd.DataFrame):
        tes_map = pd.DataFrame()
    tes_map = cat_to_obj(tes_map)
    tes_map.to_hdf(oname, 'tes_map', mode='a', format='table')

    if not isinstance(m, pd.DataFrame):
        m = pd.DataFrame()
    pdb.set_trace()
    print('you also need to fix this')
    for c in ['tss_first_sd_issue', 'tes_last_sa_issue']:
        print('# issues w/ {} nan issue: {}'.format(c, len(m.loc[m[c].isnull()].index)))
        m.loc[m[c].isnull(), c] = True
        m[c] = m[c].astype('str')
    m.to_hdf(oname, 'map', mode='a', format='table')

def write_h5_to_tsv(h5, opref):
    """
    Write a cerberus transcriptome h5 file to a series
    of flat tsv files

    Parameters:
        h5 (str): h5 file to convert to tables
        opref (str): Output file path / prefix to give
            each of the output files
    """
    ic, tss, tes, tes_map, tss_map, m = read_h5(h5, as_pyranges=True)

    oname = '{}_ic.tsv'.format(opref)
    ic.to_csv(oname, sep='\t', index=False)

    oname = '{}_tss.bed'.format(opref)
    tss.to_bed(oname)

    oname = '{}_tes.bed'.format(opref)
    tes.to_bed(oname)

    oname = '{}_tss_map.bed'.format(opref)
    tss_map.to_bed(oname)

    oname = '{}_tes_map.bed'.format(opref)
    tes_map.to_bed(oname)

    oname = '{}_map.tsv'.format(opref)
    m.to_csv(oname, sep='\t', index=False)

##### readers #####

def get_stable_gid(df, col):
    """
    Get a list of stable gene ids from a dataframe

    Parameters:
        df (pandas DataFrame): DF w/ ENSEMBL gids in some column
        col (str): Column name of gids

    Returns:
        gids (list of str): List of stable gids
    """
    df = df.copy(deep=True)
    try:
        df[['temp', 'par_region_1', 'par_region_2']] = df[col].str.split('_', n=2, expand=True)
        df[col] = df[col].str.split('.', expand=True)[0]
        df[['par_region_1', 'par_region_2']] = df[['par_region_1',
                                                           'par_region_2']].fillna('')
        df[col] = df[col]+df.par_region_1+df.par_region_2
        df.drop(['temp', 'par_region_1', 'par_region_2'], axis=1, inplace=True)
    except:
        df[col] = df[col].str.split('.', expand=True)[0]

    return df[col].tolist()


def add_stable_gid(gtf):
    """
    Add stable gene id that accounts for PAR_X and PAR_Y
    chromosomes to gtf df

    Parameters:
        gtf (pandas DataFrame): GTF dataframe

    Returns:
        gtf (pandas DataFrame): GTF dataframe with gene id turned into its
            stable version
    """
    # try:
    #     gtf[['temp', 'par_region_1', 'par_region_2']] = gtf.gene_id.str.split('_', n=2, expand=True)
    #     gtf['gene_id'] = gtf.gene_id.str.split('.', expand=True)[0]
    #     gtf[['par_region_1', 'par_region_2']] = gtf[['par_region_1',
    #                                                        'par_region_2']].fillna('')
    #     gtf['gene_id'] = gtf.gene_id+gtf.par_region_1+gtf.par_region_2
    #     gtf.drop(['temp', 'par_region_1', 'par_region_2'], axis=1, inplace=True)
    # except:
    #     gtf['gene_id'] = gtf.gene_id.str.split('.', expand=True)[0]
    gtf['gene_id'] = get_stable_gid(gtf, 'gene_id')

    return gtf

def check_files(files, sources):
    """
    Check if a series of files exists and whether their
    associated sources are all unique. Raise an exception
    if not.

    Parameters:
        files (list of str): List of file paths
        sources (list of str): List of sources
    """
    print(files)

    for f in files:
        if not os.path.exists(f):
            raise Exception('File {} does not exist'.format(f))

    if len(set(sources)) < len(sources):
        raise Exception('Sources must be unique')

def parse_gtf_config(fname):
    """
    Parse a config file for gtf end and ic calling

    Parameters:
        fname (str): Path to gtf config file

    Returns:
        ics (list of str): List of fnames to ic files
        add_ends (list of bool): List of whether or not to add ends from
            gtf file
        sources (list of str): Source names for each ic file
    """
    df = pd.read_csv(fname, header=None, sep=',')
    df.columns = ['fname', 'add_ends', 'source']

    gtfs = df.fname.tolist()
    add_ends = df.add_ends.tolist()
    sources = df.source.tolist()

    return gtfs, add_ends, sources

def read_h5(h5, as_pyranges=True):
    """
    Read h5 representation of a transcriptome

    Parameters:
        h5 (str): .h5 file to read from
        as_pyranges (bool): Convert bed representations to PyRanges objects

    Returns:
        ic (pandas DataFrame): Table detailing intron chains
        tss (pyranges PyRanges / pandas DataFrame): Bed representation of tss regions
        tes (pyranges PyRanges / pandas DataFrame): Bed represenation of tes regions
        tss_map (pyranges PyRanges / pandas DataFrame): Bed representation of
            all input tsss and which cerberus end they were assigned to
        tes_map (pyranges PyRanges / pandas DataFrame): Bed representation of
            all input tess and which cerberus end they were assigned to
        m (pandas DataFrame): Map of transcript id to tss / ic / tes
    """

    ic = pd.read_hdf(h5, key='ic')
    tss = pd.read_hdf(h5, key='tss')
    tes = pd.read_hdf(h5, key='tes')

    # turn NaN coords into empty strings
    ic.loc[ic.Coordinates.isnull(), 'Coordinates'] = ''

    def read_empty_h5(h5, key, as_pyranges=False):
        try:
            df = pd.read_hdf(h5, key=key)
            if as_pyranges:
                df = pr.PyRanges(df)
        except:
            df = None
        return df

    m = read_empty_h5(h5, 'map', as_pyranges=False)
    tss_map = read_empty_h5(h5, 'tss_map', as_pyranges=as_pyranges)
    tes_map = read_empty_h5(h5, 'tes_map', as_pyranges=as_pyranges)

    if as_pyranges:
        tss = pr.PyRanges(tss)
        tes = pr.PyRanges(tes)

    return ic, tss, tes, tss_map, tes_map, m

def get_transcript_ref(fname):
    """
    Get transcripts with necessary tags to order based on basic set etc.

    Parameters:
        fname (str): Path to GTF file

    Returns:
        df (pandas DataFrame): DataFrame of only transcript entries
            with tags parsed out
    """
    df = pr.read_gtf(fname, duplicate_attr=True).df
    df = df.loc[df.Feature == 'transcript']

    # for some transcripts, there are no tags. replace w/ empty strings
    df.loc[df.tag.isnull(), 'tag'] = ''

    df['MANE_Select'] = df.tag.str.contains('MANE_Select')
    df['basic_set'] = df.tag.str.contains('basic')
    df['temp'] = df.tag.str.split('appris_principal_', n=1, expand=True)[1]
    df['appris_principal'] = df.temp.str.split(',', n=1, expand=True)[0]
    df['appris_principal'] = df.appris_principal.astype(float)

    # if this tag didn't exist in the annotation replace it with
    # 0s so the aggregation function won't fail
    if all(df.appris_principal.isnull()):
        df.appris_principal = 0

    return df

def split_cerberus_id(df, mode):
    """
    Splits the cerberus id for a given df into its gene id
    and number for the feature

    Parameters:
        df (pandas DataFrame): DataFrame of either intron chains
            or end regions
        mode (str): {'ic', 'tss', 'tes'}

    Returns:
        df (pandas DataFrame): DataFrame with gene id and feature
            number added
    """
    df[['gene_id', mode]] = df.Name.str.split('_', expand=True)
    df[mode] = df[mode].astype(int)

    return df


def read_ic_ref(ic_file, add_gid=True, add_num=True):
    """
    Read intron chain tsv format (output from gtf_to_ics or agg_ics)

    Parameters:
        ic_file (str): Path to ic file
        add_gid (bool): Whether to include gene id in output
        add_num (bool): Whether to add intron chain # in output

    Returns:
        df (pandas DataFrame): Dataframe with gene id and intron chain number
            added in addition to existing information
    """
    df = pd.read_csv(ic_file, sep='\t')

    # add stable gene id and intron chain #
    df = split_cerberus_id(df, 'ic')

    drop_cols = []
    if not add_gid:
        drop_cols.append('gene_id')
    elif not add_num:
        drop_cols.append('ic')
    df.drop(drop_cols, axis=1, inplace=True)

    order = ['Chromosome', 'Strand', 'Coordinates',
             'Name', 'source', 'gene_id', 'ic']
    order = [o for o in order if o in df.columns]
    df = df[order]

    # turn NaN coords into empty strings
    df.loc[df.Coordinates.isnull(), 'Coordinates'] = ''

    return df

def read_cerberus_ends(bed_file, mode,
                       add_gid=True, add_num=True):
    """
    Read end reference bed file (output from gtf_to_bed or agg_ends)

    Parameters:
        bed_file (str): Path to bed file
        mode (str): {'tss', 'tes'}
        add_gid (bool): Whether to include gene id in output
        add_num (bool): Whether to add end # in output

    Returns:
        df (pandas DataFrame): Dataframe with gene id and end number
            added in addition to existing information
    """
    df = pr.read_bed(bed_file).df

    # add stable gene id and end #
    df = split_cerberus_id(df, mode)

    drop_cols = []
    if not add_gid:
        drop_cols.append('gene_id')
    elif not add_num:
        drop_cols.append(mode)
    df.drop(drop_cols, axis=1, inplace=True)

    df.rename({'ThickStart': 'source'}, axis=1, inplace=True)

    order = ['Chromosome', 'Start', 'End', 'Strand',
             'Name', 'source', 'gene_id', mode]
    order = [o for o in order if o in df.columns]
    df = df[order]

    return df

# def read_lapa_ends(bed_file,
#                        add_gid=True):
#     """
#     Read end reference bed file (output from gtf_to_bed or agg_ends)
#
#     Parameters:
#         bed_file (str): Path to bed file
#         mode (str): {'tss', 'tes'}
#         add_gid (bool): Whether to include gene id in output
#         add_num (bool): Whether to add end chain # in output
#
#     Returns:
#         df (pandas DataFrame): Dataframe with gene id and end number
#             added in addition to existing information
#     """
#     df = pr.read_bed(bed_file).df
#     df['gene_id']
#
#     drop_cols = []
#     if not add_gid:
#         drop_cols.append('gene_id')
#     df.drop(drop_cols, axis=1, inplace=True)
#
#     df.rename({'ThickStart': 'source'}, axis=1, inplace=True)
#
#     order = ['Chromosome', 'Start', 'End', 'Strand',
#              'Name', 'source', 'gene_id']
#     order = [o for o in order if o in df.columns]
#     df = df[order]
#
#     return df

def read_cerberus_source_map(fname):
    """
    Read a source map file output from `agg_ends`

    Parameters:
        fname (str): Name of source map file

    Returns:
        df (pandas DataFrame): DF bed representation of source map
    """
    df = pr.read_bed(fname).df
    df.rename({'ThickStart': 'source'}, axis=1, inplace=True)

    order = ['Chromosome', 'Start', 'End',
             'Strand', 'source', 'Name']
    order = [o for o in order if o in df.columns]
    df = df[order]

    return df

def read_bed(bed_file, mode):
    """
    Read a bed file and determine whether it comes from
    cerberus

    Parameters:
        bed_file (str): Path to input bed file
        mode (str): {'tss', 'tes'}

    Returns:
        df (pandas DataFrame): DataFrame representation of bed file
    """

    df = pr.read_bed(bed_file).df

    # bed files without strands but with a column that the
    # strand should be
    if 'Strand' in df.columns:
        ctrl = set(['+', '-'])
        if set(df.Strand.tolist()) != ctrl:
            df.rename({'Strand': 'Other'}, axis=1, inplace=True)

    # bed files output from gtf_to_bed with gid_number names
    if 'ThickStart' in df.columns:
        if 'cerberus' in df.ThickStart.tolist():
            df = read_cerberus_ends(bed_file, mode,
                                    add_gid=True, add_num=True)
        elif 'lapa' in df.ThickStart.tolist():
            df = read_cerberus_ends(bed_file, mode,
                                    add_gid=True, add_num=True)

    df = pr.PyRanges(df)
    return df

def parse_agg_ics_config(fname):
    """
    Parse a config file for intron chain aggregation

    Parameters:
        fname (str): Path to ic config file

    Returns:
        ics (list of str): List of fnames to ic files
        sources (list of str): Source names for each ic file
    """
    df = pd.read_csv(fname, header=None, sep=',')
    df.columns = ['fname', 'source']

    ics = df.fname.tolist()
    sources = df.source.tolist()

    return ics, sources

def parse_agg_ends_config(fname):
    """
    Parse a config file for end aggregation

    Parameters:
        fname (str): Path to the tss / tes config file

    Returns:
        beds (list of str): List of fnames to beds
        add_ends (list of bool): List of whether or not
            to add new ends from each bed
        sources (list of str): Source names for each bed
    """
    if not fname:
        return [], [], []

    df = pd.read_csv(fname, header=None, sep=',')
    df.columns = ['fname', 'add_ends', 'source']

    beds = df.fname.tolist()
    add_ends = df.add_ends.tolist()
    sources = df.source.tolist()

    return beds, add_ends, sources

##### main methods #####

def get_ends_from_gtf(gtf, mode, dist, slack):
    """
    Create bed regions for each tes in a gtf and number them
    based on tags indicating priority for each gene w/i the gtf

    Parameters:
        gtf (str): Path to gtf
        dist (int): Distance by which to extend regions on either side
        slack (int): Distance allowable for merging nearby regions

    Returns:
        bed (pyranges PyRanges): PyRanges object containing extended regions
    """
    gr_gtf = pr.read_gtf(gtf)

    if 'gene_id' not in gr_gtf.columns:
        raise ValueError('No gene_id field found in {}'.format(gtf))
    else:
        gr_gtf = gr_gtf.df
        gr_gtf = add_stable_gid(gr_gtf)
        gr_gtf = pr.PyRanges(gr_gtf)

    # get and extend ends
    if mode == 'tss':
        bed = gr_gtf.features.tss().extend(dist).df[[
            'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
        ]].drop_duplicates()
    elif mode == 'tes':
        bed = gr_gtf.features.tes().extend(dist).df[[
            'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
        ]].drop_duplicates()

    bed = pr.PyRanges(bed.rename({'ThickStart': 'gene_id'}, axis=1))
    bed = bed.merge(strand=True,
                by='gene_id',
                slack=slack)

    bed = number_gtf_ends(bed, gtf, mode)

    # add cerberus as the source in thickstart column
    bed = bed.df
    bed['ThickStart'] = 'cerberus'
    bed = pr.PyRanges(bed)

    return bed

def get_ics_from_gtf(gtf):
    """
    Get a file for each intron chain in a gtf and number them
    based on tags indicating priority for each gene w/i the gtf

    Parameters:
        gtf (str): Filename for input gtf

    Returns:
        ic (pandas DataFrame): Dataframe for each unique intron chain found
            in the input gtf
    """

    # get basic status and appris_principal tag for each transcript
    t_df = get_transcript_ref(gtf)
    t_df = t_df[['transcript_id', 'MANE_Select', 'basic_set', 'appris_principal']]

    # get unique intron chains from gtf
    df = pr.read_gtf(gtf)
    df = df.df
    df = add_stable_gid(df)
    df = pr.PyRanges(df)
    df = get_ic(df)

    # add basic annotation, appris principal number, and gene id
    df = df.merge(t_df, on='transcript_id', how='left')

    # add number for each unique intron chain
    df = number_tss_ic_tes(df, mode='ic')

    # make coords into tuple and perform additional
    # formatting for this table
    # df['ic'] = df.ic.str.split('-')
    # df['ic'] = [tuple(c) for c in df.ic.tolist()]
    ic = df.copy(deep=True)
    ic.rename({'ic': 'Coordinates'},
               axis=1, inplace=True)
    # ic['gene_id'] = ic.gene_id.str.split('.', n=1, expand=True)[0]
    ic['Name'] = ic['gene_id']+'_'+ic.ic_num.astype(str)
    # print(ic.loc[ic.transcript_id == 'ENCODEHT000206942'])

    cols = ['transcript_id', 'MANE_Select', 'basic_set', 'appris_principal',
            'gene_id', 'ic_num']
    ic.drop(cols, axis=1, inplace=True)


    return ic

def aggregate_ends(beds, sources, add_ends, slack, mode):
    """
    Aggregate ends from more than one bed source.

    Parameters:
        beds (list of str): List of bed file names
        sources (list of str): List of source names for each bed
        add_ends (list of bool): List of booleans indicating whether
            to add novel ends for each bed file
        slack (int): Allowable distance to an existing end for new
            ends to be called the same end
        mode (str): {'tss', 'tes'}

    Returns:
        df (pandas DataFrame): DataFrame of regions from
            aggregated bed files
        source_map (pandas DataFrame): Map ends from each source and which
            cereberus end it was mapped to
    """

    df = pd.DataFrame()
    i = 0
    for bed_fname, source, add in zip(beds, sources, add_ends):

        # read in bed file and do some formatting
        bed = read_bed(bed_fname, mode)
        bed = bed.df
        bed['source'] = source
        bed['id'] = [i for i in range(len(bed.index))]
        bed = pr.PyRanges(bed)

        # first bed; just accept all these ends
        if len(df.index) == 0:

            if not add:
                raise Exception('Must add ends from first bed file')

            if 'gene_id' not in bed.df.columns and 'Strand' not in bed.df.columns:
                raise Exception('First bed must contain Strand and gene_id columns')

            df = bed.df

            # source map
            m_source = df[['Chromosome', 'Start', 'End', 'Strand',
                           'source', 'Name']].copy(deep=True)

        # more than one bed; merge and reconcile ends
        else:

            if 'gene_id' not in bed.df.columns or 'Strand' not in bed.df.columns:
                if add:
                    print('Cannot add new ends from {} because '+\
                          'it does not contain gene_id information.'.format(source))

            # add missing columns but keep track of what information we'll be
            # able to merge on
            bed, gid, strand = format_agg_2_ends_bed(bed, mode)
            df = pr.PyRanges(df)
            df, temp = agg_2_ends(df, bed,
                            strand, gid,
                            slack, add, mode)

            # update source map
            m_source = pd.concat([m_source, temp])
            i += 1

    drop_cols = ['id', mode, 'gene_id']
    df.drop(drop_cols, axis=1, inplace=True)
    return df, m_source

def aggregate_ics(ics, sources):
    """
    Aggregate intron chains from multiple gtf_to_ics calls into one
    ics tsv table. Preferentially number intron chains based on their
    numbers computed in the order that the ic files are passed in.

    Parameters:
        ics (list of str): List of ic.tsv files in priority ordering to
            number each unique ic
        sources (list of str): List of strings describing the source of each
            ic file

    Returns:
        df (pandas DataFrame): Dataframe detailing the chromosome, strand,
            coordinates, and new names computed for each unique intron
            chain
    """
    df = pd.DataFrame()
    for ic, source in zip(ics, sources):
        temp = read_ic_ref(ic)
        temp['source'] = source

        # start with priority 1
        if len(df.index) == 0:
            df = temp.copy(deep=True)
        else:
            df = agg_2_ics(df, temp)

    # drop gene id and ic number as they are captured in name
    df.drop(['gene_id', 'ic'], axis=1, inplace=True)

    return df

def assign_triplets(gtf_df, tss, ic, tes):
    """
    Determines which tss, intron chain, and tes are used from a cerberus
    reference are used in a specific gtf.

    Parameters:
        gtf (pyranges PyRanges): PyRanges GTF object for transcriptome to assign triplets to
        ic_file (pandas DataFrame): df of intron chains from cerberus ref.
        tss_bed (str): PyRanges obj of tsss from cerberus ref
        tes_bed (str): PyRanges obj of tess from cerberus ref

    Returns:
        df (pandas DataFrame): File that maps each transcript from gtf to
            a tss, intron chain, and tes given by the cerberus reference.
            Also includes new transcript ids to refer to each transcript
    """

    ### intron chain ###

    # get intron chains from input transcriptome
    df = gtf_df.copy()
    df = get_ic(df)

    # record locations of first splice site and last splice site
    s_df = df.copy(deep=True)
    # s_df.loc[s_df.ic.isnull(), 'ic'] = ''
    s_df['coords'] = s_df.ic.str.split('-')
    first_sd = [coords[0] for coords in s_df.coords.values.tolist()]
    s_df['first_sd'] = first_sd
    last_sa = [coords[-1] for coords in s_df.coords.values.tolist()]
    s_df['last_sa'] = last_sa
    s_df.loc[s_df.first_sd == '', 'first_sd'] = None
    s_df.loc[s_df.last_sa == '', 'last_sa'] = None
    s_df[['first_sd', 'last_sa']] = s_df[['first_sd', 'last_sa']].astype(float)

    # merge ics with annotated ics
    df.rename({'ic': 'Coordinates'}, axis=1, inplace=True)

    df = merge_ics(df, ic)
    print('merged ic')

    ### ends ###
    for mode, ref in zip(['tss', 'tes'], [tss, tes]):
        if mode == 'tss':
            ends = gtf_df.features.tss()
        elif mode == 'tes':
            ends = gtf_df.features.tes()

        t_ends = merge_ends(ends, ref, mode)
        print('merge ends')

        # merge with ic ids
        df = df.merge(t_ends, how='left', on='transcript_id')

    ### creating map file ###

    # record whether or not this transcript has the bug
    # add tss / tes coords
    s_df = s_df.merge(df[['transcript_id', 'tss_id', 'tes_id']],
        on='transcript_id', how='left')
    for mode, ref in zip(['tss', 'tes'], [tss, tes]):
        temp = ref.df[['Start', 'End', 'Name']].copy(deep=True)
        temp.rename({'Start': 'Start_{}'.format(mode),
                     'End': 'End_{}'.format(mode)},
                     axis=1, inplace=True)
        s_df = s_df.merge(temp, left_on='{}_id'.format(mode), right_on='Name')

    fwd, rev = get_stranded_gtf_dfs(s_df)
    fwd['new_tss'] = fwd[['Start_tss', 'End_tss']].min(axis=1)
    fwd['new_tes'] = fwd[['Start_tes', 'End_tes']].max(axis=1)
    fwd['tss_first_sd_issue'] = fwd.new_tss > fwd.first_sd
    fwd['tes_last_sa_issue'] = fwd.new_tes <  fwd.last_sa
    rev['new_tss'] = rev[['Start_tss', 'End_tss']].max(axis=1)
    rev['new_tes'] = rev[['Start_tes', 'End_tes']].min(axis=1)
    rev['tss_first_sd_issue'] = rev.new_tss < rev.first_sd
    rev['tes_last_sa_issue'] = rev.new_tes >  rev.last_sa
    s_df = pd.concat([fwd, rev])
    s_df = s_df[['transcript_id', 'tss_first_sd_issue', 'tes_last_sa_issue']]
    s_df.rename({'transcript_id': 'original_transcript_id'}, axis=1, inplace=True)

    # get gene id / name and transcript name from original gtf
    gtf_df = gtf_df.df
    gtf_df = gtf_df.loc[gtf_df.Feature == 'transcript']
    if 'gene_name' not in gtf_df.columns:
        gtf_df['gene_name'] = gtf_df.transcript_name.str.split('-', n=1, expand=True)[0]
    gtf_df = gtf_df[['gene_id', 'gene_name',
                 'transcript_id', 'transcript_name']]
    df = df.merge(gtf_df, how='left', on='transcript_id')

    # create triplets and rename old ids
    df.rename({'transcript_id': 'original_transcript_id',
               'transcript_name': 'original_transcript_name'},
               axis=1, inplace=True)
    # pdb.set_trace()

    print('fixing issue that you need to debug later...')
    for beep in ['tss', 'tes', 'ic']:
        print('# affected transcripts w/ null {}: {}'.format(beep,len(df.loc[df[beep].isnull()].index)))
        df[beep] = df[beep].fillna(1)

    df['transcript_triplet'] = '['+df.tss.astype(int).astype(str)+','+\
                                   df.ic.astype(int).astype(str)+','+\
                                   df.tes.astype(int).astype(str)+']'
    df['transcript_id'] = df['gene_id']+df.transcript_triplet
    df['transcript_name'] = df['gene_name']+df.transcript_triplet

    df = df.merge(s_df, how='left', on='original_transcript_id')

    return df

###### routines called from main #####
####### actual routines #######
def gtf_to_bed(gtf, mode, o, dist=50, slack=50):
    bed = get_ends_from_gtf(gtf, mode, dist, slack)
    bed.to_bed(o)

def gtf_to_ics(gtf, o):
    df = get_ics_from_gtf(gtf)
    df.to_csv(o, index=False, sep='\t')

def agg_ends(input, mode, slack, o):
    beds, add_ends, sources = parse_agg_ends_config(input)
    bed, source_map = aggregate_ends(beds, sources, add_ends, slack, mode)
    bed = pr.PyRanges(bed)
    source_map = pr.PyRanges(source_map)
    bed.to_bed(o)
    o2 = get_source_map_fname(o)
    source_map.to_bed(o2)

def agg_ics(input, o):
    ics, sources = parse_agg_ics_config(input)
    ic = aggregate_ics(ics, sources)
    ic.to_csv(o, sep='\t', index=False)

def write_reference(tss_fname, tes_fname, ic_fname, o):
    def open_source_map(fname):
        try:
            map_fname = get_source_map_fname(fname)
            df = read_cerberus_source_map(map_fname)
        except:
            raise Warning('Could not find matching source map file for {} '.format(fname),
                          'Resulting entry in reference will be empty.')
            return None
        return df

    ic = read_ic_ref(ic_fname)
    tss = read_cerberus_ends(tss_fname, mode='tss')
    tes = read_cerberus_ends(tes_fname, mode='tes')

    tss_map = open_source_map(tss_fname)
    tes_map = open_source_map(tes_fname)

    write_h5(ic, tss, tes, o,
             tss_map=tss_map,
             tes_map=tes_map)

def gen_reference(ref_gtf, o, ref_tss, ref_tes,
                  gtf_tss_dist, gtf_tss_slack,
                  gtf_tes_dist, gtf_tes_slack,
                  tss_slack, tes_slack,
                  verbosity, tmp_dir,
                  keep_tmp):

    # parse config files and check if the input files exist
    gtfs, gtf_add_ends, gtf_sources = parse_gtf_config(ref_gtf)
    tss_beds, tss_add_ends, tss_sources = parse_agg_ends_config(ref_tss)
    tes_beds, tes_add_ends, tes_sources = parse_agg_ends_config(ref_tes)

    check_files(gtfs+tss_beds, gtf_sources+tss_sources)
    check_files(gtfs+tes_beds, gtf_sources+tss_sources)

    # make tmp_dir if it doesn't already exist
    if not os.path.exists(tmp_dir):
        os.makedirs(tmp_dir)

    # get ends / ics from gtfs
    gtf_tss_beds = []
    gtf_tes_beds = []
    ics = []
    for gtf, source in zip(gtfs, gtf_sources):

        if verbosity > 0:
            print('Finding TSSs from {}'.format(source))

        tss_fname = '{}/{}_tss.bed'.format(tmp_dir, source)
        gtf_to_bed(gtf, 'tss', tss_fname, gtf_tss_dist, gtf_tss_slack)

        if verbosity > 0:
            print('Finding TESs from {}'.format(source))

        tes_fname = '{}/{}_tes.bed'.format(tmp_dir, source)
        gtf_to_bed(gtf, 'tes', tes_fname, gtf_tes_dist, gtf_tes_slack)

        if verbosity > 0:
            print('Finding intron chains from {}'.format(source))

        ic_fname = '{}/{}_ic.tsv'.format(tmp_dir, source)
        ic = gtf_to_ics(gtf, ic_fname)

        gtf_tss_beds.append(tss_fname)
        gtf_tes_beds.append(tes_fname)
        ics.append(ic_fname)

    # get list of tmp files to (maybe) remove later
    tmp_files = gtf_tss_beds + gtf_tes_beds + ics

    # aggregate ends
    tss_beds = gtf_tss_beds + tss_beds
    tss_add_ends = gtf_add_ends + tss_add_ends
    tss_sources = gtf_sources + tss_sources

    if verbosity > 0:
        print('Aggregating TSSs')

    tss, tss_map = aggregate_ends(tss_beds,
                                  tss_sources,
                                  tss_add_ends,
                                  tss_slack,
                                  'tss')

    tes_beds = gtf_tes_beds + tes_beds
    tes_add_ends = gtf_add_ends + tes_add_ends
    tes_sources = gtf_sources + tes_sources

    if verbosity > 0:
        print('Aggregating TESs')

    tes, tes_map = aggregate_ends(tes_beds,
                                  tes_sources,
                                  tes_add_ends,
                                  tes_slack,
                                  'tes')

    # aggregate ics
    if verbosity > 0:
        print('Aggregating intron chains')

    ics = aggregate_ics(ics, gtf_sources)

    # write h5df
    write_h5(ics, tss, tes, o,
             tss_map=tss_map,
             tes_map=tes_map)

    # delete tmp files
    if not keep_tmp:
        for f in tmp_files:
            os.remove(f)

def convert_transcriptome(gtf, h5, o):
    """
    Create an h5 cerberus transcriptome from a GTF using an existing
    cerberus annotation

    Parameters:
        gtf (str): Path to GTF to annotate
        h5 (str): Path to cerberus annotation in h5 format
        o (str): Output .h5 file path / name
    """

    # read in / format existing reference
    ic, tss, tes, tss_map, tes_map, _ = read_h5(h5, as_pyranges=False)
    ic = split_cerberus_id(ic, 'ic')
    tss = split_cerberus_id(tss, 'tss')
    tes = split_cerberus_id(tes, 'tes')
    tss = pr.PyRanges(tss)
    tes = pr.PyRanges(tes)
    print('wow look at that im actually runing')

    # read in transcriptome to convert to cerberus
    gtf_df = pr.read_gtf(gtf).df
    print('read gtf')

    # format gid
    gtf_df = add_stable_gid(gtf_df)

    gtf_df = pr.PyRanges(gtf_df)

    df = assign_triplets(gtf_df, tss, ic, tes)

    # write h5 file
    tss = tss.df
    tes = tes.df
    write_h5(ic, tss, tes, o,
        tss_map=tss_map,
        tes_map=tes_map,
        m=df)

def replace_ab_ids(ab, h5, agg, o):
    """
    Replace the transcript ids and transcript names in a TALON abundance file
    with the new transcript ids that contain the triplet

    Parameters:
        ab (str): Path to TALON abundance file
        h5 (str): Path to h5 annotation (output from assign)
        agg (bool): Aggregate / collapse transcripts with the same triplets
            and sum up their count values

    Returns:
        df (pandas DataFrame): TALON abundance file with updated
            transcript ids / names
    """
    df = pd.read_csv(ab, sep='\t')
    _, _, _, _, _, m_df = read_h5(h5)

    # temporary fix for problematic transcripts
    rm_tids = m_df.loc[(m_df.tss_first_sd_issue)|(m_df.tes_last_sa_issue), 'original_transcript_id'].tolist()
    df = df.loc[~df.annot_transcript_id.isin(rm_tids)]

    df = map_transcripts(df, m_df, 'annot_transcript_name', 'annot_transcript_id')

    # aggregate counts if requested
    if agg:
        df = agg_ab(df)

    # reorder columns
    c1 = get_non_dataset_cols()
    c2 = get_dataset_cols(df)
    df = df[c1+c2]

    # write to file
    df.to_csv(o, sep='\t', index=False)

def replace_gtf_ids(h5, gtf, update_ends, agg, o):
    """
    Replace transcript ids in a gtf with the new cerberus ends.
    Optionally update the end coordinates of each transcript
    based on tss / tes regions and aggregate transcripts
    that use the same triplets.

    Parameters:
        h5 (str): Path to cerberus h5 file output from `convert_transcriptome`
        gtf (str): Path to GTF file to update
        update_ends (bool): Update the ends of each transcript based on regions
            in h5 tss / tes file
        agg (bool): Collapse transcripts with the same triplet and only report one
        o (str): Path to output file
    """

    if agg:
        if not update_ends:
            raise ValueError('Must update ends to aggregate transcripts')

    df = pr.read_gtf(gtf).df
    entry_types = ['gene', 'transcript', 'exon']
    df = df.loc[df.Feature.isin(entry_types)]
    df = sort_gtf(df)

    if not update_ends:
        _, _, _, _, _, m_df = read_h5(h5)
    else:
        _, tss, tes, _, _, m_df = read_h5(h5)
        tss = tss.df
        tes = tes.df
        tss['tss_id'] = tss.gene_id+'_'+tss.tss.astype(str)
        tes['tes_id'] = tes.gene_id+'_'+tes.tes.astype(str)
        tss = pr.PyRanges(tss)
        tes = pr.PyRanges(tes)

    # temporary fix for problematic transcripts
    rm_tids = m_df.loc[(m_df.tss_first_sd_issue)|(m_df.tes_last_sa_issue), 'original_transcript_id'].tolist()
    df = df.loc[~df.transcript_id.isin(rm_tids)]

    print('Adding cerberus transcript ids...')
    m_df.drop(['transcript_triplet',
               'gene_name', 'gene_id'], axis=1, inplace=True)
    df = df.merge(m_df, how='left',
                    left_on=['transcript_name', 'transcript_id'],
                    right_on=['original_transcript_name', 'original_transcript_id'],
                    suffixes=('', '_cerberus'))
    df.drop(['transcript_id', 'transcript_name'], axis=1, inplace=True)
    df.rename({'transcript_id_cerberus': 'transcript_id',
               'transcript_name_cerberus': 'transcript_name'},
               axis=1, inplace=True)

    # update the ends of each transcript based on the end it was assigned to
    if update_ends:
        print('Updating ends of transcripts...')
        df = update_gtf_ends(df, tss, tes)

    # deduplicate transcripts with the same triplets
    if agg:
        print('Deduplicating transcripts...')
        df = agg_gtf(df)

    # df.drop(['transcript_id', 'transcript_name'], axis=1, inplace=True)
    # df.rename({'transcript_id_cerberus': 'transcript_id',
    #             'transcript_name_cerberus': 'transcript_name'},
    #            axis=1, inplace=True)

    # write gtf
    df = pr.PyRanges(df)
    df.to_gtf(o)

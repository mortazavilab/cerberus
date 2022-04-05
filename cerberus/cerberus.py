import pyranges as pr
import pandas as pd
import pdb

def get_nov_ranks():
    """
    Get rank for each novelty category

    Returns:
        nov_rank (dict): Dictionary of novelty category: int rank
        rank_nov (dict): Dictionary of int rank: novelty category
    """
    nov_rank = {'Known': 0,
           'NIC': 1,
           'ISM': 2,
           'NNC': 3,
           'Genomic': 4,
           'Antisense': 5,
           'Intergenic': 6}

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
    df['basic_set'] = df.tag.str.contains('basic')
    df['temp'] = df.tag.str.split('appris_principal_', n=1, expand=True)[1]
    df['appris_principal'] = df.temp.str.split(',', n=1, expand=True)[0]

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
    cols = ['transcript_id', mode, 'basic_set', 'appris_principal', 'gene_id']
    df = df[cols].groupby([mode, 'gene_id']).agg({'transcript_id': ','.join,
                                     'basic_set': 'max',
                                     'appris_principal': 'min'}).reset_index()

    # compute intron chain number
    df['{}_num'.format(mode)] = df.sort_values(by=['gene_id', 'appris_principal', 'basic_set'],
                                 ascending=[True, False, True],
                                 na_position='last')\
                                 .groupby(['gene_id'])\
                                 .cumcount() + 1

    return df

def get_ends_from_gtf(gtf, mode, dist, slack):
    """
    Create bed regions for each tes in a gtf

    Parameters:
        gtf (str): Path to gtf
        dist (int): Distance by which to extend regions on either side
        slack (int): Distance allowable for merging nearby regions

    Returns:
        bed (pyranges PyRanges): PyRanges object containing extended regions
    """
    gr_gtf = pr.read_gtf(gtf)
    if 'gene_id' not in gr_gtf.df.columns:
        raise ValueError('No gene_id field found in {}'.format(gtf))

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

    return bed

def aggregate_ends(beds, mode):
    """
    Aggregate the end regions within the same gene
    from multiple bed files.

    Parameters:
        beds (str): Comma-separated list of paths to bed files
            or path of text file with path of each bed file on each line
            Bed files will be aggregated in the order that they are given
        mode (str): 'tss' or 'tes'
        slack (int): Allowable distance b/w ends to cluster them.
            Default: 50

    Returns:
        bed (pyranges PyRanges): PyRanges object containing aggregated ends
    """

    # bed files given as list on cmd line
    if beds.endswith('.bed'):
        beds = beds.split(',')

    # bed files given in a file
    else:
        df = pd.read_csv(beds, header=None, names=['bed_file'])
        beds = df.bed_file.tolist()

    if len(beds) > 1:
        raise ValueError('Currently more than one bed file is not supported')

    for bed_fname in beds:
        bed = pr.read_bed(bed_fname)

    return bed

def add_triplets(gtf, tss_bed, tes_bed):
    """
    Merges end regions with ends from GTF and assigns triplets to each detected
    intron chain, tss, and tes for each gene.
    Preferentially orders based on transcript support level (appris_principal)
    and which transcripts are in the basic annotation

    Parameters:
        gtf (str): File path to GTF
        tss_bed (str): File path to TSS bed; output from `aggregate_ends` or
            `get_ends_from_gtf`
        tes_bed (str): File path to TES bed: output from `aggregate_ends` or
            `get_ends_from_gtf`

    Returns:
        gtf (pyranges PyRanges): GTF PyRanges object with added field
            "transcript_triplet"
        tss_bed (pyranges PyRanges): BED PyRanges object with TSS named
            gene_id_tss_#
        tes_bed (pyranges PyRanges): BED PyRanges object with TES named
            gene_id_tes_#
        df (pandas DataFrame): DataFrame including new and old transcript IDs
    """

    ### intron chain annotation ###

    # get basic status and appris_principal tag for each transcript
    t_df = get_transcript_ref(gtf)
    t_df = t_df[['transcript_id', 'gene_id', 'basic_set', 'appris_principal']]

    # get unique intron chains from gtf
    df = pr.read_gtf(gtf).df

    # restrict to exon entries
    df = df.loc[df.Feature == 'exon']
    cols = ['Chromosome', 'Strand', 'Start', 'End', 'transcript_id', 'gene_id']
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
                    ascending=[True, True, True])
    rev.sort_values(by=['Chromosome', 'transcript_id', 'Coord'],
                    ascending=[True, True, False])
    df = pd.concat([fwd, rev])

    # create intron chain strings
    df.Coord = df.Coord.astype(str)
    df = df.groupby(['Chromosome', 'Strand', 'transcript_id'])['Coord'].apply('-'.join).reset_index()

    # remove tss and tes from intron chain
    df['temp'] = df.Coord.str.split('-', n=1, expand=True)[1]
    df['ic'] = df.temp.str.rsplit('-', n=1, expand=True)[0]

    # add basic annotation, appris principal number, and gene id
    df = df.merge(t_df, on='transcript_id', how='left')

    # add number for each unique intron chain
    df = number_tss_ic_tes(df, mode='ic')

    # keep track of which transcripts use each intron chain
    df.transcript_id = df.transcript_id.str.split(',')
    ttrip_df = df[['transcript_id', 'ic_num']].explode(column='transcript_id')

    ### tss / tes numbering ###
    beds = {}
    for mode, bed_file in zip(['tss', 'tes'], [tss_bed, tes_bed]):

        ends = pr.read_bed(bed_file).df
        ends.rename({'ThickStart': 'gene_id'}, axis=1, inplace=True)
        ends[mode] = [i for i in range(len(ends.index))]
        ends = pr.PyRanges(ends)

        # merge gtf ends with bed ends
        df = pr.read_gtf(gtf)
        if mode == 'tss':
            df = df.features.tss().df
        elif mode == 'tes':
            df = df.features.tes().df

        # merge in information about transcript status
        t_df = get_transcript_ref(gtf)
        t_df = t_df[['transcript_id', 'appris_principal', 'basic_set']]
        df = df.merge(t_df, how='left', on='transcript_id')
        df = pr.PyRanges(df)

        # join with bed ends
        df = df.join(ends,
                     strandedness='same',
                     how='left').df
        df = df.loc[df.gene_id == df.gene_id_b]

        # add number for each unique end region
        df = number_tss_ic_tes(df, mode=mode)

        # keep track of which transcripts use each intron chain
        df.transcript_id = df.transcript_id.str.split(',')
        temp = df[['transcript_id', '{}_num'.format(mode)]].explode(column='transcript_id')
        ttrip_df = ttrip_df.merge(temp, on='transcript_id')

        # add numbers back into bed file
        ends = ends.df
        ends = ends.merge(df[[mode, '{}_num'.format(mode)]], on=mode)

        ends['Name'] = ends.gene_id+'_'+ends['{}_num'.format(mode)].astype(str)
        ends.drop([mode, 'gene_id', '{}_num'.format(mode)], axis=1, inplace=True)
        beds[mode] = pr.PyRanges(ends)

    # create formatted string to represent each transcript
    ttrip_df['transcript_triplet'] = '['+ttrip_df['tss_num'].astype(str)+', '+\
                                         ttrip_df['ic_num'].astype(str)+', '+\
                                         ttrip_df['tes_num'].astype(str)+']'

    # final gtf file formatting
    df = pr.read_gtf(gtf, duplicate_attr=True).df
    df = df.merge(ttrip_df[['transcript_id', 'transcript_triplet']],
                  how='left',
                  on='transcript_id')
    gtf = pr.PyRanges(df)

    # transcript to triplet reference file

    # add gene name from transcript name
    df = df.loc[df.Feature == 'transcript']
    if 'gene_name' not in df.columns:
        df['gene_name'] = df.transcript_name.str.split('-', n=1, expand=True)[0]
    df = df[['transcript_id', 'transcript_name',
             'gene_name', 'gene_id',
             'transcript_triplet']]
    df.rename({'transcript_id': 'original_transcript_id',
               'transcript_name': 'original_transcript_name'},
               axis=1, inplace=True)
    df['transcript_id'] = df.gene_id+' '+df.transcript_triplet
    df['transcript_name'] = df.gene_name+' '+df.transcript_triplet

    return gtf, beds['tss'], beds['tes'], df

def replace_gtf_ids(gtf, m, agg):
    """
    Replace transcript ids and names in a gtf with the triplets
    calculated from assign_triplets

    Parameters:
        gtf (str): Path to gtf file
        m (str): Path to map file output from assign_triplets
        agg (bool): Whether or not to collapse transcripts with
            duplicate triplets

    Returns:
        df (pyranges PyRanges): PyRanges gtf table with updated ids
    """

    df = pr.read_gtf(gtf).df
    m_df = pd.read_csv(m, sep='\t')

    # groupby transcripts that are the same
    gb_cols = ['gene_name', 'gene_id', 'transcript_triplet',
               'transcript_id', 'transcript_name']
    temp = m_df[['transcript_id',
                 'original_transcript_id',
                 'original_transcript_name']].copy(deep=True)
    m_df = m_df.groupby(gb_cols).agg({'original_transcript_id': ','.join,
                                      'original_transcript_name': ','.join}).reset_index()
    m_df = m_df.merge(temp, on='transcript_id', suffixes=('','_merge'))
    m_df.drop(['gene_name', 'gene_id', 'transcript_triplet'],
              axis=1, inplace=True)

    # add new transcript ids
    df = df.merge(m_df, left_on=['transcript_id', 'transcript_name'],
                  right_on=['original_transcript_id_merge',
                            'original_transcript_name_merge'],
                 suffixes=('_x', ''))

    # drop old tids
    df.drop(['transcript_id_x', 'transcript_name_x',
             'original_transcript_id_merge',
             'original_transcript_name_merge'],
            axis=1, inplace=True)

    # if we're only reporting one transcript per triplet
    if agg:
        df = df.loc[~((df.duplicated(subset=['transcript_id',
                                             'transcript_name',
                                             'exon_id']))&(df.Feature.isin(['transcript', 'exon'])))]

    df = pr.PyRanges(df)

    return df

def replace_ab_ids(ab, m, agg):
    """
    Replace the transcript ids and transcript names in a TALON abundance file
    with the new transcript ids that contain the triplet

    Parameters:
        ab (str): Path to TALON abundance file
        m (str): Path to map file (output from assign_triplets)
        agg (bool): Aggregate / collapse transcripts with the same triplets
            and sum up their count values

    Returns:
        df (pandas DataFrame): TALON abundance file with updated
            transcript ids / names
    """
    df = pd.read_csv(ab, sep='\t')
    m_df = pd.read_csv(m, sep='\t')

    # fix transcript ids in abundance file
    ab_map = m_df[['original_transcript_id', 'original_transcript_name',
                   'transcript_name', 'transcript_id']]
    df = df.merge(ab_map, left_on='annot_transcript_id', right_on='original_transcript_id')

    df.drop(['original_transcript_id', 'original_transcript_name',
             'annot_transcript_id', 'annot_transcript_name'],
            axis=1, inplace=True)
    df.rename({'transcript_name': 'annot_transcript_name',
               'transcript_id': 'annot_transcript_id'},
              axis=1, inplace=True)

    # aggregate counts if requested
    if agg:
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

    # reorder columns
    c1 = get_non_dataset_cols()
    c2 = get_dataset_cols(df)
    df = df[c1+c2]

    return df

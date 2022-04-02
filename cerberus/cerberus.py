import pyranges as pr
import pandas as pd
import pdb

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
        df (pyranges PyRanges): GTF PyRanges object with added field
            "transcript_triplet"
        tss_bed (pyranges PyRanges): BED PyRanges object with TSS named
            gene_id_tss_#
        tes_bed (pyranges PyRanges): BED PyRanges object with TES named
            gene_id_tes_#
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
    df = pr.PyRanges(df)

    return df, beds['tss'], beds['tes']

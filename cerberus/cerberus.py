import pyranges as pr
import pandas as pd
import pdb

def get_tes_from_gtf(gtf, dist):
    """
    Create bed regions for each tes in a gtf

    Parameters:
        gtf (str): Path to gtf
        dist (int): Distance by which to extend regions on either side

    Returns:
        bed (pyranges PyRanges): PyRanges object containing extended regions
    """
    gr_gtf = pr.read_gtf(gtf)

    df_tes = gr_gtf.features.tes().extend(dist).df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
    ]].drop_duplicates()

    return df_tes

def get_tss_from_gtf(gtf, dist):
    """
    Create bed regions for each tss in a gtf

    Parameters:
        gtf (str): Path to gtf
        dist (int): Distance by which to extend regions on either side

    Returns:
        bed (pyranges PyRanges): PyRanges object containing extended regions
    """
    gr_gtf = pr.read_gtf(gtf)

    df_tss = gr_gtf.features.tss().extend(dist).df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
    ]].drop_duplicates()

    return df_tss

def aggregate_ends(beds, mode, slack):
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
        if 'ThickStart' not in bed.df.columns:
            raise ValueError('No gene_id column found in {}'.format(bed_fname))
        bed = pr.PyRanges(bed.df.rename({'ThickStart': 'gene_id'}, axis=1))
        bed = bed.merge(strand=True,
                        by='gene_id',
                        slack=slack)

    return bed

import pyranges as pr
import pandas as pd
import pdb

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

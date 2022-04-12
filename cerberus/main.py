import click
from cerberus.cerberus import *

@click.group()
def cli():
    pass

@cli.command()
@click.option('--mode',
              help='Choose tss or tes',
              required=True)
@click.option('--gtf',
              help='GTF file',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
@click.option('--dist',
              help='Distance (bp) to extend regions on either side',
              default=50)
@click.option('--slack',
              help='Distance allowable for merging regions',
              default=50)
def gtf_to_bed(mode, gtf, o, dist=50, slack=50):
    bed = get_ends_from_gtf(gtf, mode, dist, slack)
    bed.to_bed(o)

@cli.command()
@click.option('--gtf',
              help='GTF file',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def gtf_to_ics(gtf, o):
    df = get_ics_from_gtf(gtf)
    df.to_csv(o, index=False, sep='\t')

@cli.command()
@click.option('--mode',
              help='Choose tss or tes',
              required=True)
@click.option('--input',
              help='Path to file w/ path to BED '+\
                'files on each line or comma-separated '+\
                ' list of file paths; ordered by priority',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def agg_ends(mode, input, o):
    bed = aggregate_ends(input, mode)
    bed.to_bed(o)

@cli.command()
@click.option('--gtf',
              help='GTF of isoforms to assign triplets to',
              required=True)
@click.option('--ic',
              help='Intron chain file',
              required=True)
@click.option('--tss_bed',
              help='Bed file of TSS regions',
              required=True)
@click.option('--tes_bed',
              help='Bed file of TES regions',
              required=True)
@click.option('-o',
              help='Output file name',
              required=True)
def assign_triplets(gtf, ic, tss_bed, tes_bed, o):
    df = add_triplets(gtf, ic, tss_bed, tes_bed)

    # read in references that we'll just dump back to a new
    # h5 file
    tss_bed = pr.read_bed(tss_bed).df
    tes_bed = pr.read_bed(tes_bed).df
    ic = pd.read_csv(ic, sep='\t')

    ic = change_all_dtypes(ic, str)
    tss_bed = change_all_dtypes(tss_bed, str)
    tes_bed = change_all_dtypes(tes_bed, str)

    ic.to_hdf(o, 'ic', mode='w')
    tss_bed.to_hdf(o, 'tss', mode='a', format='table')
    tes_bed.to_hdf(o, 'tes', mode='a', format='table')
    df.to_hdf(o, 'map', mode='a')

@cli.command()
@click.option('--map',
              help='transcript ID map from assign_triplets',
              required=True)
@click.option('--gtf',
              help='GTF of isoforms',
              required=False,
              default=None)
@click.option('--ab',
              help='TALON abundance file',
              required=False,
              default=None)
@click.option('--collapse',
              help='collapse transcripts with the same triplets',
              is_flag=True,
              required=False,
              default=False)
@click.option('--opref',
              help='Output file prefix to save updated gtf / ab',
              required=True)
def replace_ids(map, gtf, ab, collapse, opref):
    if gtf:
        df = replace_gtf_ids(gtf, map, collapse)
        oname = '{}.gtf'.format(opref)
        df.to_gtf(oname)
    if ab:
        df = replace_ab_ids(ab, map, collapse)
        oname = '{}.tsv'.format(opref)
        df.to_csv(oname, index=False, sep='\t')

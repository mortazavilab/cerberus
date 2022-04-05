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
              help='GTF of isoforms',
              required=True)
@click.option('--tss_bed',
              help='Bed file of TSS regions',
              required=True)
@click.option('--tes_bed',
              help='Bed file of TES regions',
              required=True)
@click.option('--opref',
              help='Output file prefix to save beds / gtf w/ triplets',
              required=True)
def assign_triplets(gtf, tss_bed, tes_bed, opref):
    gtf, tss_bed, tes_bed, df = add_triplets(gtf, tss_bed, tes_bed)

    oname = '{}.gtf'.format(opref)
    gtf.to_gtf(oname)

    oname = '{}_tss.bed'.format(opref)
    tss_bed.to_bed(oname)

    oname = '{}_tes.bed'.format(opref)
    tes_bed.to_bed(oname)

    oname = '{}_tid_map.tsv'.format(opref)
    df.to_csv(oname, index=False, sep='\t')

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

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
                'files on each line; ordered by priority',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def agg_ends(mode, input, o):
    bed = aggregate_ends(input, mode)
    bed.to_bed(o)

def compute_triplets(mode, gtf, ofile, dist=50):
    click.echo('Syncing')

def update_talon_abundance(mode, gtf, ofile, dist=50):
    click.echo('Syncing')

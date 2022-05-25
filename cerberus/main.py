import click
from cerberus.cerberus import *

@click.group()
def cli():
    pass

####### create a reference #######
@cli.command(name='gtf_to_bed')
@click.option('--gtf',
              help='GTF file',
              required=True)
@click.option('--mode',
              help='Choose tss or tes',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
@click.option('--dist',
              help='Distance (bp) to extend regions on either side',
              default=50,
              show_default=True)
@click.option('--slack',
              help='Distance allowable for merging regions',
              default=50,
              show_default=True)
def gtf_to_bed_command(gtf, mode, o, dist=50, slack=50):
    return gtf_to_bed(gtf, mode, o, dist, slack)

@cli.command(name='gtf_to_ics')
@click.option('--gtf',
              help='GTF file',
              required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def gtf_to_ics_command(gtf, o):
    return gtf_to_ics(gtf, o)

@cli.command(name='agg_ends')
@click.option('--input',
              help='Path to config file. Each line contains'+\
                   'file path,whether to add ends (True / False),source name',
              required=True)
@click.option('--mode',
            help='Choose tss or tes',
            required=True)
@click.option('--slack',
              help='Distance (bp) allowable for merging regions',
              default=20,
              show_default=True)
@click.option('-o',
            help='Output file name',
            required=True)
def agg_ends_command(input, mode, slack, o):
    return agg_ends(input, mode, slack, o)

@cli.command(name='agg_ics')
@click.option('--input',
              help='Path to config file. Each line contains'+\
                   'file path,source name',
              required=True)
@click.option('-o',
              help='Output file name',
              required=True)
def agg_ics_command(input, o):
    return agg_ics(input, o)

# wrapper to create a reference
@cli.command(name='gen_reference')
@click.option('--ref_gtf',
              help='Path to config file for each GTF, in priority order. '+\
                   'Each line contains file path, whether to add ends, source name',
              required=True)
@click.option('-o',
              help='Output h5 file name',
              required=True)
@click.option('--ref_tss',
              help='Path to config file for each TSS bed file, '+\
                   'in priority order. Each line contains file path, '+\
                   'whether to add ends, source. Note: these beds will be ordered '+\
                   'after the beds from the GTF file. ',
              required=False)
@click.option('--ref_tes',
            help='Path to config file for each TES bed file, '+\
                 'in priority order. Each line contains file path, '+\
                 'whether to add ends, source. Note: these beds will be ordered '+\
                 'after the beds from the GTF file. ',
            required=False)
@click.option('--gtf_tss_dist',
              help='Distance (bp) to extend TSS regions on either side of '+\
                   'each end from the GTFs',
              default=50,
              show_default=True)
@click.option('--gtf_tss_slack',
              help='Distance allowable for merging TSS regions from each GTF',
              default=50,
              show_default=True)
@click.option('--gtf_tes_dist',
            help='Distance (bp) to extend TES regions on either side of '+\
                 'each end from the GTFs',
            default=50,
            show_default=True)
@click.option('--gtf_tes_slack',
            help='Distance allowable for merging TSS regions from each GTF',
            default=50,
            show_default=True)
@click.option('--tss_slack',
              help='Distance (bp) allowable for merging TSS regions',
              default=20,
              show_default=True)
@click.option('--tes_slack',
            help='Distance (bp) allowable for merging TES regions',
            default=20)
@click.option('--verbosity',
              help='Verbosity. Higher numbers mean more output.',
              default=1,
              show_default=True)
@click.option('--tmp_dir',
              help='Prefix / file path to save temporary files.',
              default='temp',
              show_default=True)
@click.option('--keep_tmp',
              help='Keep intermediate bed and ic files instead of deleting them',
              default=False,
              is_flag=True)
def gen_reference_command(ref_gtf, o,
                          ref_tss, ref_tes,
                          gtf_tss_dist, gtf_tss_slack,
                          gtf_tes_dist, gtf_tes_slack,
                          tss_slack, tes_slack,
                          verbosity, tmp_dir,
                          keep_tmp):
    gen_reference(ref_gtf, o, ref_tss, ref_tes,
                  gtf_tss_dist, gtf_tss_slack,
                  gtf_tes_dist, gtf_tes_slack,
                  tss_slack, tes_slack,
                  verbosity, tmp_dir,
                  keep_tmp)

####### annotate a transcriptome ######
@cli.command(name='convert_transcriptome')
@click.option('--gtf',
              help='GTF file',
              required=True)
@click.option('--h5',
            help='cerberus reference from gen_reference',
            required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def convert_transcriptome_command(gtf, h5, o):
    convert_transcriptome(gtf, h5, o)

###### replace transcript ids of an annotated transcriptome
###### with cerberus triplet ids
@cli.command(name='replace_ab_ids')
@click.option('--h5',
              help='cerberus reference from gen_reference',
              required=True)
@click.option('--ab',
              help='TALON abundance file to replace ids in',
              required=True)
@click.option('--collapse',
              help='collapse transcripts with the same triplets',
              is_flag=True,
              required=False,
              default=False)
@click.option('-o',
              help='Output file name',
              required=True)
def replace_ab_ids_command(h5, ab, collapse, o):
    replace_ab_ids(ab, h5, collapse, o)


# @cli.command()
# @click.option('--h5',
#               help='h5 file output from assign-triplets')
# @click.option('--gtf',
#               help='GTF of isoforms',
#               required=False,
#               default=None)
# @click.option('--ab',
#               help='TALON abundance file',
#               required=False,
#               default=None)
# @click.option('--collapse',
#               help='collapse transcripts with the same triplets',
#               is_flag=True,
#               required=False,
#               default=False)
# @click.option('--opref',
#               help='Output file prefix to save updated gtf / ab',
#               required=True)
# def replace_ids(h5, gtf, ab, collapse, opref):
#     if gtf:
#         df = replace_gtf_ids(gtf, h5, collapse)
#         oname = '{}.gtf'.format(opref)
#         df.to_gtf(oname)
#     if ab:
#         df = replace_ab_ids(ab, h5, collapse)
#         oname = '{}.tsv'.format(opref)
#         df.to_csv(oname, index=False, sep='\t')
#
# @cli.command(name='h5_to_tsv')
# @click.option('--h5',
#               help='h5 transcriptome file output from cerberus assign-triplets',
#               required=True)
# @click.option('--opref',
#               help='output file prefix',
#               required=True)
# def h5_to_tsv(h5, opref):
#     write_h5_to_tsv(h5, opref)

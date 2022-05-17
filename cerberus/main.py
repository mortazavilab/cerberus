import click
from cerberus.cerberus import *

@click.group()
def cli():
    pass

####### actual routines #######
def gtf_to_bed(gtf, mode, o, dist=50, slack=50):
    bed = get_ends_from_gtf(gtf, mode, dist, slack)
    bed.to_bed(o)

def gtf_to_ics(gtf, o):
    df = get_ics_from_gtf(gtf)
    df.to_csv(o, index=False, sep='\t')

def agg_ends(input, mode, slack, o):
    beds, add_ends, sources = parse_agg_ends_config(input)
    bed = aggregate_ends(beds, sources, add_ends, slack, mode)
    bed = pr.PyRanges(bed)
    bed.to_bed(o)

def agg_ics(input, o):
    ics, sources = parse_agg_ics_config(input)
    ic = aggregate_ics(ics, sources)
    ic.to_csv(o, sep='\t', index=False)

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
            print('Finding TSSs from {}'.format(source))

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

    tss = aggregate_ends(tss_beds,
                         tss_sources,
                         tss_add_ends,
                         tss_slack,
                         'tss')

    tes_beds = gtf_tes_beds + tes_beds
    tes_add_ends = gtf_add_ends + tes_add_ends
    tes_sources = gtf_sources + tes_sources

    if verbosity > 0:
        print('Aggregating TESs')

    tes = aggregate_ends(tes_beds,
                         tes_sources,
                         tes_add_ends,
                         tes_slack,
                         'tes')

    # aggregate ics
    if verbosity > 0:
        print('Aggregating intron chains')

    ics = aggregate_ics(ics, gtf_sources)

    # write h5df
    write_h5(ics, tss, tes, o)

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
    ic, tss, tes, _ = read_h5(h5, as_pyranges=False)
    ic = split_cerberus_id(ic, 'ic')
    tss = split_cerberus_id(tss, 'tss')
    tes = split_cerberus_id(tes, 'tes')
    tss = pr.PyRanges(tss)
    tes = pr.PyRanges(tes)

    # read in transcriptome to convert to cerberus
    gtf_df = pr.read_gtf(gtf).df
    gtf_df['gene_id'] = gtf_df.gene_id.str.split('.', expand=True)[0]
    gtf_df = pr.PyRanges(gtf_df)

    df = assign_triplets(gtf_df, tss, ic, tes)

    # write h5 file
    tss = tss.df
    tes = tes.df
    write_h5(ic, tss, tes, o, m=df)

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
    return gtf_to_bed(gtf, mode, dist, slack)

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
            help='cerberus reference from gen-reference',
            required=True)
@click.option('-o',
            help='Output file name',
            required=True)
def convert_transcriptome_command(gtf, h5, o):
    convert_transcriptome(gtf, h5, o)


@cli.command()
@click.option('--h5',
              help='h5 file output from assign-triplets')
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
def replace_ids(h5, gtf, ab, collapse, opref):
    if gtf:
        df = replace_gtf_ids(gtf, h5, collapse)
        oname = '{}.gtf'.format(opref)
        df.to_gtf(oname)
    if ab:
        df = replace_ab_ids(ab, h5, collapse)
        oname = '{}.tsv'.format(opref)
        df.to_csv(oname, index=False, sep='\t')

@cli.command(name='h5_to_tsv')
@click.option('--h5',
              help='h5 transcriptome file output from cerberus assign-triplets',
              required=True)
@click.option('--opref',
              help='output file prefix',
              required=True)
def h5_to_tsv(h5, opref):
    write_h5_to_tsv(h5, opref)

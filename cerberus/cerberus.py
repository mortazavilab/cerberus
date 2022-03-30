import pyranges as pr


def get_tes_from_gtf(gtf, dist):
    gr_gtf = pr.read_gtf(gtf)

    df_tes = gr_gtf.features.tes().extend(dist).df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
    ]].drop_duplicates()

    return df_tes

def get_tss_from_gtf(gtf, dist):
    gr_gtf = pr.read_gtf(gtf)

    df_tss = gr_gtf.features.tss().extend(dist).df[[
        'Chromosome', 'Start', 'End', 'Strand', 'gene_id',
    ]].drop_duplicates()

    return df_tss

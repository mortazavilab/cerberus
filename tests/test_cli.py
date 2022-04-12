import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import gtf_to_bed, gtf_to_ics, agg_ends, assign_triplets, replace_ids

def test_gtf_to_bed(tmp_path):
    bed_path = str(tmp_path/'test.bed')
    runner = CliRunner()
    cmd = '--mode tss --gtf {} -o {} --dist 50 --slack 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed, cmd)
    assert result.exit_code == 0

    bed = pr.read_bed(bed_path)
    # import pdb
    # pdb.set_trace()
    # bed.head()
    # assert all(bed.lengths() == 101)

    cmd = '--mode tes --gtf {} -o {} --dist 50 --slack 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed, cmd)
    assert result.exit_code == 0

    bed = pr.read_bed(bed_path)
    # import pdb
    # pdb.set_trace()
    # bed.head()
    # assert all(bed.lengths() == 101)

def test_gtf_to_ics(tmp_path):
    out_path = str(tmp_path/'test.tsv')
    runner = CliRunner()
    cmd = '--gtf {} -o {}'.format(canx_gtf, out_path)
    print(cmd)
    result = runner.invoke(gtf_to_ics, cmd)
    assert result.exit_code == 0


def test_agg_ends(tmp_path):
    ofile = str(tmp_path/'test.bed')
    runner = CliRunner()

    # 1 file, passed from cmd
    cmd = '--mode tss --input tests/files/Canx_tss.bed -o {}'.format(ofile)
    print(cmd)
    result = runner.invoke(agg_ends, cmd)
    assert result.exit_code == 0

    # 1 file, passed from config file
    cmd = '--mode tss --input tests/files/Canx_tss_beds.txt -o {}'.format(ofile)
    print(cmd)
    result = runner.invoke(agg_ends, cmd)
    assert result.exit_code == 0

    # >1 file, passed from cmd

def test_assign_triplets(tmp_path):
    opref = str(tmp_path/'test')
    runner = CliRunner()

    # saving as series of 4 files
    cmd = '--gtf {} --tss_bed {} --tes_bed {} --opref {}'.format(canx_gtf, canx_tss_bed, canx_tes_bed, opref)
    print(cmd)
    result = runner.invoke(assign_triplets, cmd)
    assert result.exit_code == 0

    # saving as h5 file
    cmd = '--gtf {} --tss_bed {} --tes_bed {} --opref {} --h5'.format(canx_gtf, canx_tss_bed, canx_tes_bed, opref)
    print(cmd)
    result = runner.invoke(assign_triplets, cmd)
    assert result.exit_code == 0

def test_replace_ids(tmp_path):
    opref = str(tmp_path/'test')
    runner = CliRunner()

    # ab+gtf
    cmd = '--map {} --gtf {} --ab {} --collapse --opref {}'.format(canx_tid_map,
                                                             canx_gtf,
                                                             canx_ab,
                                                             opref)
    print(cmd)
    result = runner.invoke(replace_ids, cmd)
    assert result.exit_code == 0

    # ab
    cmd = '--map {} --ab {} --collapse --opref {}'.format(canx_tid_map,
                                                             canx_ab,
                                                             opref)
    print(cmd)
    result = runner.invoke(replace_ids, cmd)
    assert result.exit_code == 0

    # gtf
    cmd = '--map {} --gtf {} --collapse --opref {}'.format(canx_tid_map,
                                                             canx_gtf,
                                                             opref)
    print(cmd)
    result = runner.invoke(replace_ids, cmd)
    assert result.exit_code == 0

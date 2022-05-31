import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import *
from cerberus.main import *

def test_gtf_to_bed(tmp_path):
    bed_path = str(tmp_path/'test.bed')
    runner = CliRunner()
    cmd = '--mode tss --gtf {} -o {} --dist 50 --slack 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed_command, cmd)
    assert result.exit_code == 0

    bed = pr.read_bed(bed_path)
    # import pdb
    # pdb.set_trace()
    # bed.head()
    # assert all(bed.lengths() == 101)

    cmd = '--mode tes --gtf {} -o {} --dist 50 --slack 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed_command, cmd)
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
    result = runner.invoke(gtf_to_ics_command, cmd)
    assert result.exit_code == 0


def test_agg_ends(tmp_path):
    ofile = str(tmp_path/'test.bed')
    runner = CliRunner()

    # 1 file, passed from config file
    cmd = '--input tests/files/agg_tss_config.csv --mode tss --slack 20 -o {}'.format(ofile)
    print(cmd)
    result = runner.invoke(agg_ends_command, cmd)
    assert result.exit_code == 0


# def test_agg_ics(tmp_path):
#
#     ofile = str(tmp_path/'test_ic.tsv')
#     runner = CliRunner()
#
#     # 1 file, passed from cmd
#     cmd = '--input tests/files/Canx_ics.tsv -o {}'.format(ofile)
#     print(cmd)
#     result = runner.invoke(agg_ics_command, cmd)
#     assert result.exit_code == 0
#
#     # 1 file, passed from config file
#     cmd = '--input tests/files/Canx_ic_files.txt -o {}'.format(ofile)
#     print(cmd)
#     result = runner.invoke(agg_ics_command, cmd)
#     assert result.exit_code == 0
#
#     # >1 file, passed from cmd
#     cmd = '--input tests/files/Canx_1_ics.tsv,tests/files/Canx_2_ics.tsv -o {}'.format(ofile)
#     print(cmd)
#     result = runner.invoke(agg_ics_command, cmd)
#     assert result.exit_code == 0
#
#     # >1 file, passed from config file
#     cmd = '--input tests/files/Canx_ic_files_2.txt -o {}'.format(ofile)
#     print(cmd)
#     result = runner.invoke(agg_ics_command, cmd)
#     assert result.exit_code == 0

# def test_convert_transcriptome(tmp_path):
#     o = str(tmp_path/'test.h5')
#     runner = CliRunner()
#
#     # saving as h5 file
#     cmd = '--gtf {} --h5 {} -o {}'.format(canx_gtf, , o)
#     print(cmd)
#     result = runner.invoke(convert_transcriptome, cmd)
#     assert result.exit_code == 0
#
# def test_replace_ids(tmp_path):
#     opref = str(tmp_path/'test')
#     runner = CliRunner()
#
#     # ab+gtf
#     cmd = '--h5 {} --gtf {} --ab {} --collapse --opref {}'.format(canx_h5,
#                                                              canx_gtf,
#                                                              canx_ab,
#                                                              opref)
#     print(cmd)
#     result = runner.invoke(replace_ids, cmd)
#     assert result.exit_code == 0
#
#     # ab
#     cmd = '--h5 {} --ab {} --collapse --opref {}'.format(canx_h5,
#                                                              canx_ab,
#                                                              opref)
#     print(cmd)
#     result = runner.invoke(replace_ids, cmd)
#     assert result.exit_code == 0
#
#     # gtf
#     cmd = '--h5 {} --gtf {} --collapse --opref {}'.format(canx_h5,
#                                                              canx_gtf,
#                                                              opref)
#     print(cmd)
#     result = runner.invoke(replace_ids, cmd)
#     assert result.exit_code == 0
#
# def test_h5_to_tsv(tmp_path):
#     opref = str(tmp_path/'test')
#     runner = CliRunner()
#
#     cmd = '--h5 {} --opref {}'.format(canx_h5, opref)
#     print(cmd)
#     result = runner.invoke(h5_to_tsv, cmd)
#     assert result.exit_code == 0

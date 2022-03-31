import pyranges as pr
import pandas as pd
import pdb
from click.testing import CliRunner
from conftest import canx_gtf
from cerberus.main import gtf_to_bed, agg_ends

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

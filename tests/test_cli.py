import pyranges as pr
import pandas as pd
from click.testing import CliRunner
from conftest import canx_gtf
from cerberus.main import gtf_to_bed

def test_gtf_to_bed(tmp_path):
    bed_path = str(tmp_path/'test.bed')
    print(bed_path)
    runner = CliRunner()
    cmd = '--mode tss --gtf {} -o {} --dist 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed, cmd)
    assert result.exit_code == 0

    bed = pr.read_bed(bed_path)
    # import pdb
    # pdb.set_trace()
    # bed.head()
    assert all(bed.lengths() == 101)

    cmd = '--mode tes --gtf {} -o {} --dist 50'.format(canx_gtf, bed_path)
    print(cmd)
    result = runner.invoke(gtf_to_bed, cmd)
    assert result.exit_code == 0

    bed = pr.read_bed(bed_path)
    import pdb
    pdb.set_trace()
    # bed.head()
    assert all(bed.lengths() == 101)

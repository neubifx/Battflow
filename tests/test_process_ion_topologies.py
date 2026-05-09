# tests/test_process_ion_topologies.py

from pathlib import Path
from battflow.md_setup import process_ion_topologies

def test_process_ion_topologies(tmp_path):

    pack_path = tmp_path / "pack"
    md_em_path = tmp_path / "em"
    md_eq_path = tmp_path / "eq"
    md_prod_path = tmp_path / "prod"

    pack_path.mkdir()
    md_em_path.mkdir()
    md_eq_path.mkdir()
    md_prod_path.mkdir()

    BASE_DIR = Path(__file__).resolve().parents[1] / "battflow"

    ions_itp_files, topol_main_file, ions_pdb = process_ion_topologies(
        BASE_DIR=BASE_DIR,
        config={},
        ions=["Li"],
        pack_path=pack_path,
        md_em_path=md_em_path,
        md_eq_path=md_eq_path,
        md_prod_path=md_prod_path
    )

    assert len(ions_pdb) == 1
    assert (pack_path / "li.pdb").exists()
    assert (md_em_path / "li.itp").exists()

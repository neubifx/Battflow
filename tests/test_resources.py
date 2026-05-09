from pathlib import Path
import battflow

def test_default_yaml_exists():

    base = Path(battflow.__file__).parent

    assert (base / "config" / "default.yaml").exists()
# tests/test_user_config.py

from pathlib import Path
from battflow.database import config_path

def test_user_config(tmp_path):

    config_file = tmp_path / "test.yaml"

    config_file.write_text("""
mongodb:
    host: localhost
""")

    BASE_DIR, config = config_path(config_file)

    assert config["mongodb"]["host"] == "localhost"

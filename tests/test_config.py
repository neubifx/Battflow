from battflow.database import config_path

def test_config_loading():
    BASE_DIR, config = config_path()

    assert "mongodb" in config
    assert "md_run_env" in config
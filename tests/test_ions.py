import subprocess

def test_cli_help():

    result = subprocess.run(
        ["battflow", "--help"],
        capture_output=True,
        text=True
    )

    assert result.returncode == 0
    assert "Run Battflow workflow" in result.stdout
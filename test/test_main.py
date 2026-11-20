import os
import shutil
import subprocess

import pytest

import sequana_pipelines.superin.main as m  # noqa: F401
from . import test_dir


def test_version():
    cmd = ["sequana_superin", "--version"]
    subprocess.check_call(cmd)


def test_help():
    cmd = ["sequana_superin", "--help"]
    subprocess.check_call(cmd)


def test_standalone_subprocess(tmp_path):
    fq = tmp_path / "fq"
    plasmids = tmp_path / "plasmids"
    for d in (fq, plasmids):
        d.mkdir()
    (fq / "s1.fastq.gz").touch()
    (plasmids / "p1.fa").touch()
    (tmp_path / "SI.fa").touch()
    (tmp_path / "ss.csv").write_text("sample,plasmid,reference\ns1,p1,p1.fa\n")

    cmd = [
        "sequana_superin",
        "--input-directory", str(fq),
        "--samplesheet", str(tmp_path / "ss.csv"),
        "--plasmids-directory", str(plasmids),
        "--SI-reference-file", str(tmp_path / "SI.fa"),
        "--working-directory", str(tmp_path / "wd"),
        "--force",
    ]
    subprocess.check_call(cmd)
    assert (tmp_path / "wd" / "superin.rules").exists()


def _run_pipeline(tmp_path, extra_options=()):
    data = os.path.join(test_dir, "data")
    wkdir = tmp_path / "wd"
    cmd = [
        "sequana_superin",
        "--input-directory", os.path.join(data, "fastq"),
        "--input-bam-pattern", os.path.join(data, "bam"),
        "--samplesheet", "samplesheet.csv",  # relative path on purpose (rulegraph runs in a subdirectory)
        "--plasmids-directory", os.path.join(data, "references"),
        "--SI-reference-file", os.path.join(data, "references", "SI.fasta"),
        "--working-directory", str(wkdir),
        "--force",
        *extra_options,
    ]
    shutil.copy(os.path.join(data, "samplesheet.csv"), tmp_path)
    subprocess.check_call(cmd, cwd=tmp_path)
    subprocess.check_call(["bash", "superin.sh"], cwd=wkdir)

    assert (wkdir / "nanocap_filtered_fastqs" / "stats.csv").exists()
    assert (wkdir / "plasmids" / "deletions_insertions_stats.csv").exists()
    assert (wkdir / "SI" / "position_deletions.png").exists()
    assert (wkdir / ".sequana" / "rulegraph.svg").exists()
    # references are copied in the working directory for book keeping
    for name in ("SI.fasta", "pV827.fa", "pV830.fa"):
        assert (wkdir / "references" / name).exists()


def test_full(tmp_path):
    """Run the whole pipeline with tools found in the PATH (no containers)"""
    _run_pipeline(tmp_path)


@pytest.mark.skipif(
    shutil.which("apptainer") is None or not os.environ.get("SUPERIN_TEST_APPTAINER"),
    reason="needs apptainer and SUPERIN_TEST_APPTAINER=1 (images are downloaded)",
)
def test_full_apptainer(tmp_path):
    """Run the whole pipeline with the containers (images downloaded on first use)

    The image cache can be set with SUPERIN_APPTAINER_PREFIX to avoid downloading
    the images each time.
    """
    prefix = os.environ.get("SUPERIN_APPTAINER_PREFIX", str(tmp_path / "images"))
    _run_pipeline(tmp_path, ["--apptainer-prefix", prefix])

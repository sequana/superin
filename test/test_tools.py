import pysam

from sequana_pipelines.superin import tools


def test_get_INDELs_uses_contig_of_the_read(tmp_path):
    """Deleted sequences must come from the contig the read is mapped on, not from the first contig."""
    contig1 = "A" * 2000
    # contig2 has a deleted region (positions 500-700) with a recognisable, non-A content
    block = "ACGTTGCAAGCT"
    deleted = (block * 20)[:200]
    contig2 = "C" * 500 + deleted + "G" * 500

    fasta = tmp_path / "ref.fa"
    fasta.write_text(f">chr1\n{contig1}\n>plasmid\n{contig2}\n")

    # read made of the 500 first bases and the last 500 bases of contig2 (200 bp deletion)
    query = "C" * 500 + "G" * 500
    header = {"HD": {"VN": "1.0"}, "SQ": [{"SN": "chr1", "LN": len(contig1)}, {"SN": "plasmid", "LN": len(contig2)}]}
    bamfile = tmp_path / "test.bam"
    with pysam.AlignmentFile(str(bamfile), "wb", header=header) as fout:
        a = pysam.AlignedSegment()
        a.query_name = "read1"
        a.query_sequence = query
        a.flag = 0
        a.reference_id = 1
        a.reference_start = 0
        a.mapping_quality = 60
        a.cigarstring = "500M200D500M"
        a.query_qualities = pysam.qualitystring_to_array("I" * len(query))
        fout.write(a)

    df, stats = tools.get_INDELs(str(bamfile), str(fasta), minimal_length_indels=100)

    assert stats["deletions"] == 1
    assert df.iloc[0]["typ"] == "D"
    assert df.iloc[0]["sequence"] == deleted

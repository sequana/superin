# Regenerates test/data/{bam,fastq} from the full B6060-3 data set (paths are local to the author).
import gzip, re, subprocess, sys
import pysam

src = "/home/cokelaer/Work/BiomicsProjects/B6060-3/data"
out = "/home/cokelaer/Work/github/forked/sequana_pipelines/superin/test/data"
samples = {"Y305": "pV827", "Y308": "pV830", "Y306": "pV827"}
N_PLASMID, N_SI_INDEL, N_SI_PLAIN = 60, 60, 30
indel = re.compile(r"\d{3,}[ID]")

for sample, plasmid in samples.items():
    inbam = pysam.AlignmentFile(f"{src}/bam/{sample}.sorted.bam")
    keep = []
    n = 0
    for r in inbam.fetch(plasmid):
        if r.flag in (0, 16) and r.cigarstring and indel.search(r.cigarstring) and r.query_length < 6000:
            keep.append(r); n += 1
            if n >= N_PLASMID: break
    n_i = n_p = 0
    SI0 = 302056
    cig = re.compile(r"(\d+)([MIDNSHP=X])")

    def si_indels(r):
        """positions (SI coordinates) of indels longer than 100"""
        pos = r.reference_start - SI0
        found = []
        for l, op in cig.findall(r.cigarstring):
            l = int(l)
            if op in "M=XN":
                pos += l
            elif op == "D":
                if l > 100: found.append(pos)
                pos += l
            elif op == "I" and l > 100:
                found.append(pos)
        return found

    for r in inbam.fetch("NC_002506.1", SI0 + 10600, 448942):
        if r.flag not in (0, 16) or r.query_length > 6000 or r.reference_end > 448942 or r.reference_start < SI0 + 10000:
            continue
        pos = si_indels(r)
        if pos and min(pos) > 10600 and n_i < N_SI_INDEL:
            keep.append(r); n_i += 1
        elif not pos and n_p < N_SI_PLAIN:
            keep.append(r); n_p += 1
        if n_i >= N_SI_INDEL and n_p >= N_SI_PLAIN: break
    print(sample, len(keep), n_i, n_p, file=sys.stderr)
    keep.sort(key=lambda r: (r.reference_id, r.reference_start))
    with pysam.AlignmentFile(f"{out}/bam/{sample}.sorted.bam", "wb", header=inbam.header) as o:
        for r in keep: o.write(r)
    pysam.index(f"{out}/bam/{sample}.sorted.bam")
    names = {r.query_name for r in keep}
    found = 0
    with gzip.open(f"{src}/fastq/{sample}.fastq.gz", "rt") as fin, gzip.open(f"{out}/fastq/{sample}.fastq.gz", "wt") as fout:
        while True:
            h = fin.readline()
            if not h: break
            rest = [fin.readline() for _ in range(3)]
            if h[1:].split()[0] in names:
                fout.write(h + "".join(rest)); found += 1
                if found == len(names): break
    print(sample, "fastq found", found, "of", len(names), file=sys.stderr)

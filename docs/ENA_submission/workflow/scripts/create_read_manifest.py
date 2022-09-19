"""https://ena-docs.readthedocs.io/en/latest/submit/reads/webin-cli.html#manifest-file"""
import sys

sys.stderr = open(snakemake.log[0], "w")

from pathlib import Path
import pandas as pd

MADAGASCAR = "madagascar"
INDIA = "india"


def get_country(sample: str) -> str:
    if sample.startswith("P"):
        return MADAGASCAR
    elif sample.startswith("Test"):
        return INDIA
    else:
        raise ValueError(f"Got unknown sample {sample}")


def sample_name(s: str) -> str:
    for sfx in ["_rep1", "_rep2", "-singleplex", "-multiplex"]:
        s = s.replace(sfx, "")
    return s


def main():
    sample_acc_df = pd.read_csv(snakemake.input.samples, index_col="alias")
    project = snakemake.params.project
    alias = snakemake.wildcards.sample
    sample = sample_name(alias)
    biosample = sample_acc_df.at[sample, "biosample"]
    tech = snakemake.wildcards.tech
    source = snakemake.wildcards.source
    reads_dir = Path(snakemake.params.reads_dir)

    if tech == "nanopore":
        instrument = "MinION"
        fastq = [reads_dir / f"{tech}/filtered/{alias}/{alias}.filtered.fq.gz"]
    elif tech == "illumina":
        site = get_country(alias)
        if site == MADAGASCAR:
            instrument = "Illumina HiSeq 2500"
        elif site == INDIA:
            instrument = "NextSeq 500"
        else:
            raise KeyError(f"Unrecognised site {site}")
        fastq = [
            reads_dir / f"{tech}/filtered/{source}/{alias}/{alias}_R{i}.filtered.fq.gz"
            for i in [1, 2]
        ]
    else:
        raise ValueError(f"Unknown technology {tech}")

    for f in map(Path, fastq):
        assert f.exists(), f

    library_selection = "unspecified"
    library_source = "GENOMIC"
    library_strategy = "WGS" if source == "culture" else "METAGENOMIC"

    with open(snakemake.output.manifest, "w") as fp:
        print(f"STUDY\t{project}", file=fp)
        print(f"SAMPLE\t{biosample}", file=fp)
        alias = f"{alias}-{tech}-{source}"
        print(f"NAME\t{alias}", file=fp)
        print(f"INSTRUMENT\t{instrument}", file=fp)
        print(f"LIBRARY_SOURCE\t{library_source}", file=fp)
        print(f"LIBRARY_SELECTION\t{library_selection}", file=fp)
        print(f"LIBRARY_STRATEGY\t{library_strategy}", file=fp)
        for p in fastq:
            print(f"FASTQ\t{p}", file=fp)


main()

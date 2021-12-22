"""This script generates a HTML file with a table containing information about the
composition and lineage of each sample.
"""
import re
import sys
from itertools import takewhile

sys.stderr = open(snakemake.log[0], "w")
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import TextIO, Tuple, Optional, List, Union
import json
import jinja2
import pandas as pd

NORD0 = "#2e3440"
NORD11 = "#bf616a"
NORD12 = "#d08770"


class RipgrepError(Exception):
    pass


def ripgrep_search(file: Path) -> Tuple[int, int, int]:
    pattern = "\[INFO\]:\s(?P<num>\d+)\sread"
    extra_params = ["--replace", "$num", "--only-matching", "--no-line-number"]
    process = subprocess.Popen(
        ["rg", *extra_params, pattern, str(file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise RipgrepError(
            f"Failed to execute the rg command on the file {file} due to the "
            f"following error:\n{process.stderr.read()}"
        )
    values = map(int, map(str.rstrip, process.stdout.readlines()))
    to_keep, contam, unmapped = list(values)
    return to_keep, contam, unmapped


def ripgrep_extract_covg(file: Path) -> float:
    pattern = r"Actual cov.*\s(?P<covg>\d+?\.?\d+)x"
    extra_params = ["--replace", "$covg", "--only-matching", "--no-line-number"]
    process = subprocess.Popen(
        ["rg", *extra_params, pattern, str(file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise RipgrepError(
            f"Failed to execute the rg command on the file {file} due to the "
            f"following error:\n{process.stderr.read()}"
        )
    return float(process.stdout.read().strip())


LINEAGE_DELIM = "."
LINEAGE_REGEX = re.compile(
    rf"^(lineage)?(?P<major>[\w]+)\{LINEAGE_DELIM}?(?P<minor>.*)?$"
)


class InvalidLineageString(Exception):
    pass


class Lineage:
    """A class representing a lineage and any of its sublineage information."""

    def __init__(
        self,
        major: str,
        minor: Optional[Union[str, List[str]]] = None,
        minor_delim: str = ".",
    ):
        self.major = major
        if minor is None:
            minor = tuple()
        if isinstance(minor, str):
            minor = tuple(minor.split(minor_delim))
        self.minor = tuple(minor)
        self._minor_delim = minor_delim

    def __str__(self):
        if not self.minor:
            return self.major
        return self.major + self._minor_delim + self._minor_delim.join(self.minor)

    @staticmethod
    def from_str(s: str) -> "Lineage":
        """Create a Lineage object from a str."""
        match = LINEAGE_REGEX.search(s)
        if not match:
            return Lineage.from_name(s)
        major = match.group("major")
        minor = match.group("minor") or None
        return Lineage(major=major, minor=minor)

    @staticmethod
    def from_name(s: str) -> "Lineage":
        if "beijing" in s.lower():
            lin = 2
        elif "european" in s.lower():
            lin = 4
        elif "east_africa" in s.lower():
            lin = 1
        elif "central_asia" in s.lower():
            lin = 3
        else:
            lin = s
        return Lineage(lin)

    def __eq__(self, other: "Lineage") -> bool:
        return self.major == other.major and self.minor == other.minor

    def __lt__(self, other: "Lineage") -> bool:
        if (not self.minor and not other.minor) or not self.minor:
            return False
        if not other.minor:
            return True
        # we now know that both have minors
        return len(self.minor) > len(other.minor)

    def mrca(self, other: "Lineage") -> Optional["Lineage"]:
        """Determine the most recent common ancestor between two Lineages.
        Returns None if there is no MRCA.
        """
        if self.major != other.major:
            return Lineage("Mixed")
        if not self.minor or not other.minor:
            return Lineage(major=self.major)

        try:
            common_minors, _ = zip(
                *takewhile(lambda xy: xy[0] == xy[1], zip(self.minor, other.minor))
            )
            minor_str = LINEAGE_DELIM.join(common_minors)
        except ValueError:
            minor_str = None

        return Lineage(self.major, minor_str)

    @staticmethod
    def call(lineages: List["Lineage"]) -> Optional["Lineage"]:
        """Returns the Lineage with the most specific minor.
        if the majors are different, returns None. If minors are the same, then takes
        MRCA of the lineages with the same minor length.
        """
        if not lineages:
            return Lineage("Unknown")
        if len(lineages) == 1:
            return lineages[0]
        lineages.sort()
        minors_of_same_len = filter(
            lambda l: len(l.minor) == len(lineages[0].minor), lineages
        )
        lineage = lineages[0]
        for lin in minors_of_same_len:
            lineage = lin.mrca(lineage)

        return lineage


contam_warning = snakemake.params.get("contam_warning", 5.0) / 100
unmapped_warning = snakemake.params.get("unmapped_warning", 5.0) / 100
covg_warning = snakemake.params.get("covg_warning", 30)


def highlight_high_contam(col: pd.Series):
    """Highlights cells in a column if their contamination level is over the
    threshold
    """
    return [
        f"background-color: {NORD12}" if val > contam_warning else "" for val in col
    ]


def highlight_high_unmapped(col: pd.Series):
    """Highlights cells if their unmapped level is over the threshold"""
    return [
        f"background-color: {NORD12}" if val > unmapped_warning else "" for val in col
    ]


def highlight_low_coverage(col: pd.Series):
    """Highlights cells if their coverage is less than the threshold"""
    return [f"background-color: {NORD12}" if val < covg_warning else "" for val in col]


def highlight_abnormal_lineages(col: pd.Series):
    """Highlights cells if their lineage is not one of the numbered majors."""
    return [f"background-color: {NORD12}" if v.major.isalpha() else "" for v in col]


def highlight_abnormal_species(col: pd.Series):
    """Highlights cells if their lineage is not one of the numbered majors."""
    return [
        f"background-color: {NORD12}" if "tuberculosis" not in v else "" for v in col
    ]


data = defaultdict(dict)
logfiles = snakemake.input.filter_logs
for file in map(Path, logfiles):
    sample = file.name.split(".")[0]
    num_keep, num_contam, num_unmapped = ripgrep_search(file)
    total = sum([num_keep, num_contam, num_unmapped])
    data[sample].update(
        {
            f"keep": num_keep,
            f"keep%": num_keep / total,
            f"contam": num_contam,
            f"contam%": num_contam / total,
            f"unmapped": num_unmapped,
            f"unmapped%": num_unmapped / total,
            f"total": total,
        }
    )

for file in map(Path, snakemake.input.subsample_logs):
    covg = ripgrep_extract_covg(file)
    sample = file.name.split(".")[0]
    data[sample]["coverage"] = round(covg, 1)

assignment_files = snakemake.input.lineages
for file in map(Path, assignment_files):
    sample = file.name.split(".")[0]
    with open(file) as fp:
        phylo = json.load(fp)[sample]["phylogenetics"]
        lineage_calls = phylo["lineage"]
        if "lineage" in lineage_calls:
            lineages = [Lineage.from_str(l) for l in lineage_calls["lineage"]]
        else:
            lineages = [Lineage.from_str(l) for l in lineage_calls.keys()]

        data[sample]["lineage"] = Lineage.call(lineages)
        species = list(phylo["species"].keys())
        if len(species) != 1:
            raise ValueError(f"Unexpected number of species {species}")
        species = species[0]
        data[sample]["species"] = species

df = pd.DataFrame(data).T
df.index.name = "sample"

percent_format_cols = [s for s in data[list(data.keys())[0]].keys() if s.endswith("%")]
df_styled = (
    df.style.apply(highlight_high_contam, subset=["contam%"])
    .apply(highlight_abnormal_lineages, subset=["lineage"])
    .apply(highlight_abnormal_species, subset=["species"])
    .apply(highlight_low_coverage, subset=["coverage"])
    .apply(highlight_high_unmapped, subset=["unmapped%"])
    .format("{:.2%}", subset=percent_format_cols)
)

table_html = df_styled.render()
template = snakemake.params.template
template_content = Path(template).read_text()
html = jinja2.Template(template_content).render(
    table=table_html,
    contam_warning=contam_warning,
    unmapped_warning=unmapped_warning,
    covg_warning=covg_warning,
)

with open(snakemake.output.html, "w") as outfile:
    outfile.write(html)

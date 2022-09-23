import sys


sys.stderr = open(snakemake.log[0], "w")

import xml.etree.ElementTree as ET
from pathlib import Path
from typing import Tuple, List
import re

ALIAS_REGEX = re.compile(r"(?P<alias>mada_[12]-?\d+|R\d+|1[6-9]_\d+)")


def extract_accessions_from_receipt(file: Path) -> Tuple[str, str, str, str, str]:
    tree = ET.parse(file)
    root = tree.getroot()

    if root.attrib["success"] != "true":
        raise ValueError("Failed response XML given")

    experiment = None
    run_acc = None

    for child in root:
        if child.tag == "EXPERIMENT":
            experiment = child.attrib["accession"]
        if child.tag == "RUN":
            run_acc = child.attrib["accession"]
            alias = child.attrib["alias"]

    if experiment is None:
        raise KeyError(f"Got no experiment accession for successful submission {file}")
    if run_acc is None:
        raise KeyError(f"Got no run accession for successful submission {file}")

    tech = alias.split("-")[-2]
    source = alias.split("-")[-1]
    name = alias.replace("webin-reads-", "").rsplit("-", maxsplit=2)[0]
    return name, tech, source, experiment, run_acc


failed: List[Path] = []
passed: List[Tuple[str, str, str, str, str]] = []

for d in map(Path, snakemake.input.receipts):
    xmls = list(d.rglob("*receipt.xml"))
    if len(xmls) != 1:
        raise ValueError(f"{d} has too many receipts")
    receipt = xmls[0]
    try:
        data = extract_accessions_from_receipt(receipt)
        passed.append(data)
    except ValueError:
        failed.append(receipt)

if failed:
    for f in failed:
        print(f"[ERROR]: {f} failed submission", file=sys.stderr)
    sys.exit(1)

with open(snakemake.output.sheet, "w") as fp:
    print("alias,tech,source,experiment,run", file=fp)
    for fields in passed:
        print(",".join(fields), file=fp)

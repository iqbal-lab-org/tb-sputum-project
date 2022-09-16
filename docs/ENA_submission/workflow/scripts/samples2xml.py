"""
This scripts takes a samplesheet and creates an ENA-compliant Sample XML file.
See https://ena-docs.readthedocs.io/en/latest/submit/samples/programmatic.html and
https://ena-docs.readthedocs.io/en/latest/submit/samples.html
USAGE: python samples2xml.py samplesheet.csv
"""
# Good example for Sample entry https://www.ebi.ac.uk/ena/browser/view/SAMEA95444668
import sys

sys.stderr = open(snakemake.log[0], "w")

import xml.etree.ElementTree as ET

# https://www.ebi.ac.uk/ena/browser/view/ERC000028
CHECKLIST = "ERC000028"
TAXON_ID = ET.Element("TAXON_ID")
TAXON_ID.text = "1773"
SCI_NAME = ET.Element("SCIENTIFIC_NAME")
SCI_NAME.text = "Mycobacterium tuberculosis"
TITLE = ET.Element("TITLE")
TITLE.text = "Mycobacterium tuberculosis WGS"


def attribute(tag: str, value) -> ET.Element:
    el = ET.Element("SAMPLE_ATTRIBUTE")
    tag_el = ET.Element("TAG")
    tag_el.text = tag
    el.append(tag_el)
    value_el = ET.Element("VALUE")
    value_el.text = value
    el.append(value_el)
    return el


def attributes(loc: str) -> ET.Element:
    attrs = ET.Element("SAMPLE_ATTRIBUTES")
    attrs.append(attribute("collection date", snakemake.params.date))
    attrs.append(attribute("host health state", "diseased"))
    attrs.append(attribute("host scientific name", "homo sapien"))
    attrs.append(attribute("isolate", "Mycobacterium tuberculosis"))
    attrs.append(attribute("isolation_source", "sputum"))
    attrs.append(attribute("geographic location (country and/or sea)", loc))
    attrs.append(attribute("ENA-CHECKLIST", CHECKLIST))

    return attrs


def add_sample_name(el: ET.Element):
    name = ET.Element("SAMPLE_NAME")
    name.append(TAXON_ID)
    name.append(SCI_NAME)
    el.append(name)


def location(sample: str) -> str:
    if sample.startswith("P"):
        return "Madagascar"
    elif sample.startswith("Test"):
        return "India"
    else:
        raise ValueError(f"Can't infer location for {sample}")


def main():
    samples = snakemake.params.samples
    sys.stdout = open(snakemake.output.xml, "w")
    print('<?xml version="1.0" encoding="UTF-8"?>')
    root = ET.Element("SAMPLE_SET")
    for sample in samples:
        sample_el = ET.Element("SAMPLE", attrib=dict(alias=sample))
        sample_el.append(TITLE)
        add_sample_name(sample_el)
        loc = location(sample)
        attrs = attributes(loc=loc)
        sample_el.append(attrs)
        root.append(sample_el)

    ET.indent(root)
    ET.dump(root)
    sys.stdout.close()


if __name__ == "__main__":
    main()

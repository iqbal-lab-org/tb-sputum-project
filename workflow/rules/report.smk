rule plot_nanopore_sample_composition:
    input:
        tsv=rules.generate_nanopore_krona_input.output.krona_input,
    output:
        chart=report(
            ont_results / "plots/krona/{sample}.krona.html",
            category="Krona",
            subcategory="Nanopore",
            caption=report_dir / "krona.rst",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * GB,
    conda:
        str(env_dir / "krona.yaml")
    log:
        rule_log_dir / "plot_nanopore_sample_composition/{sample}.log",
    shell:
        "ktImportText {input.tsv} -o {output.chart} &> {log}"


rule plot_illumina_sample_composition:
    input:
        tsv=rules.generate_illumina_krona_input.output.krona_input,
    output:
        chart=report(
            illumina_results / "plots/krona/{source}/{sample}.krona.html",
            category="Krona",
            subcategory="Illumina {source}",
            caption=report_dir / "krona.rst",
        ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * GB,
    conda:
        str(env_dir / "krona.yaml")
    log:
        rule_log_dir / "plot_illumina_sample_composition/{source}/{sample}.log",
    shell:
        "ktImportText {input.tsv} -o {output.chart} &> {log}"


rule nanopore_composition_report:
    input:
        filter_logs=expand(
            rule_log_dir / "filter_nanopore_contamination/{sample}.log",
            sample=ont_samplesheet.index,
        ),
        lineages=expand(
            ont_results / "amr_predictions/{sample}.mykrobe.json",
            sample=ont_samplesheet.index,
        ),
        subsample_logs=expand(
            rule_log_dir / "subsample_nanopore_reads/{sample}.log",
            sample=ont_samplesheet.index,
        ),
    output:
        html=report(
            report_dir / "nanopore_composition.html",
            category="Composition",
            subcategory="Nanopore",
            caption=report_dir / "composition.rst",
        ),
    threads: 1
    resources:
        mem_mb=GB,
    conda:
        str(env_dir / "composition_report.yaml")
    params:
        template=report_dir / "composition.html.jinja",
        contam_warning=5.0,
        unmapped_warning=5.0,
        covg_warning=30,
    log:
        rule_log_dir / "nanopore_composition_report.log",
    script:
        str(scripts_dir / "composition_report.py")


def infer_illumina_filter_logs(wildcards):
    src = wildcards.source
    subdf = illumina_samplesheet.query("source==@src")
    return [
        rule_log_dir / f"filter_illumina_contamination/{src}/{sample}.log"
        for sample in subdf["sample"]
    ]


def infer_illumina_lineages_input(wildcards):
    src = wildcards.source
    subdf = illumina_samplesheet.query("source==@src")
    return [
        illumina_results / f"amr_predictions/{src}/{sample}.mykrobe.json"
        for sample in subdf["sample"]
    ]


def infer_illumina_subsample_logs(wildcards):
    src = wildcards.source
    subdf = illumina_samplesheet.query("source==@src")
    return [
        rule_log_dir / f"subsample_illumina_reads/{src}/{sample}.log"
        for sample in subdf["sample"]
    ]


rule illumina_composition_report:
    input:
        filter_logs=infer_illumina_filter_logs,
        lineages=infer_illumina_lineages_input,
        subsample_logs=infer_illumina_subsample_logs,
    output:
        html=report(
            report_dir / "illumina_{source}_composition.html",
            category="Composition",
            subcategory="Illumina",
            caption=report_dir / "composition.rst",
        ),
    threads: 1
    resources:
        mem_mb=GB,
    conda:
        str(env_dir / "composition_report.yaml")
    params:
        template=report_dir / "composition.html.jinja",
        contam_warning=5.0,
        unmapped_warning=5.0,
        covg_warning=20,
    log:
        rule_log_dir / "illumina_composition_report/{source}.log",
    script:
        str(scripts_dir / "composition_report.py")


rule nanopore_amr_report:
    input:
        predictions=expand(
            ont_results / "amr_predictions/{sample}.mykrobe.json",
            sample=ont_samplesheet.index,
        ),
        template=report_dir / "amr.html.jinja",
    output:
        report=report(
            report_dir / "nanopore_amr.html",
            category="AMR predictions",
            subcategory="Nanopore",
            caption=report_dir / "amr.rst",
        ),
    resources:
        mem_mb=int(4 * GB),
    params:
        script=scripts_dir / "amr_report.py",
    log:
        rule_log_dir / "nanopore_amr_report.log",
    conda:
        str(env_dir / "amr_report.yaml")
    shell:
        "python {params.script} {input.template} {input.predictions} 2> {log} > {output.report}"


rule illumina_amr_report:
    input:
        predictions=infer_illumina_lineages_input,
        template=report_dir / "amr.html.jinja",
    output:
        report=report(
            report_dir / "illumina_{source}_amr.html",
            category="AMR predictions",
            caption=report_dir / "amr.rst",
        ),
    resources:
        mem_mb=int(4 * GB),
    params:
        script=scripts_dir / "amr_report.py",
    log:
        rule_log_dir / "illumina_amr_report/{source}.log",
    conda:
        str(env_dir / "amr_report.yaml")
    shell:
        "python {params.script} {input.template} {input.predictions} 2> {log} > {output.report}"

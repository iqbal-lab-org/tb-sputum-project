rule mykrobe_nanopore:
    input:
        reads=rules.subsample_nanopore_reads.output.reads,
    output:
        report=ont_results / "amr_predictions/{sample}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(4 * GB),
    container:
        containers["mykrobe"]
    log:
        rule_log_dir / "mykrobe_nanopore/{sample}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "--debug",
                "--format json",
                "-e 0.08",
                "--min_proportion_expected_depth 0.20",
                "--ploidy haploid",
                "--species tb",
            ]
        ),
    threads: 4
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.reads} \
          --sample {wildcards.sample} \
          -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """

rule mykrobe_illumina:
    input:
        r1=rules.subsample_illumina_reads.output.r1,
        r2=rules.subsample_illumina_reads.output.r2,
    output:
        report=illumina_results / "amr_predictions/{source}/{sample}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(4 * GB),
    container:
        containers["mykrobe"]
    log:
        rule_log_dir / "mykrobe_illumina/{source}/{sample}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "--debug",
                "--format json",
                "--min_proportion_expected_depth 0.20",
                "--ploidy haploid",
                "--species tb",
            ]
        ),
    threads: 4
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.r1} {input.r2} \
          --sample {wildcards.sample} \
          -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """
rule build_decontamination_db:
    output:
        fasta=decontam_db / "remove_contam.fa.gz",
        metadata=decontam_db / "remove_contam.tsv",
    threads: 1
    resources:
        mem_mb=GB,
    params:
        script=scripts_dir / "download_tb_reference_files.pl",
        outdir=lambda wildcards, output: Path(output.fasta).parent,
    conda:
        str(env_dir / "decontam_db.yaml")
    log:
        rule_log_dir / "build_decontamination_db.log",
    shell:
        """
        perl {params.script} {params.outdir} &> {log}
        tmpfile=$(mktemp)
        sed 's/NTM\t0/NTM\t1/g' {output.metadata} > "$tmpfile"
        mv "$tmpfile" {output.metadata}
        """


rule index_decontam_db_minimap2:
    input:
        fasta=rules.build_decontamination_db.output.fasta,
    output:
        mm2_index=decontam_db / "remove_contam.fa.gz.map-ont.mmi",
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
    params:
        extras="-I 12G -x map-ont",
    log:
        rule_log_dir / "index_decontam_db_minimap2.log",
    conda:
        str(env_dir / "aln_tools.yaml")
    shell:
        """
        minimap2 {params.extras} \
            -t {threads} \
            -d {output.mm2_index} \
            {input.fasta} 2> {log}
        """


rule index_decontam_db_bwa:
    input:
        fasta=rules.build_decontamination_db.output.fasta,
    output:
        bwa_index=multiext(str(decontam_db / "remove_contam.fa.gz"), *BWA_EXTNS),
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(32 * GB),
    log:
        rule_log_dir / "index_decontam_db_bwa.log",
    conda:
        str(env_dir / "aln_tools.yaml")
    shell:
        "bwa index {input.fasta} 2> {log}"


rule illumina_preprocessing:
    input:
        r1=illumina_dir / "{source}/{sample}_R1.fq.gz",
        r2=illumina_dir / "{source}/{sample}_R2.fq.gz",
    output:
        r1=illumina_results / "preprocessing/{source}/{sample}_R1.fq.gz",
        r2=illumina_results / "preprocessing/{source}/{sample}_R2.fq.gz",
        report=report(
            illumina_results / "preprocessing/{source}/{sample}.fastp.html",
            category="QC",
            caption=report_dir / "illumina_preprocessing.rst",
        ),
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(4 * GB),
    log:
        rule_log_dir / "illumina_preprocessing/{source}/{sample}.log",
    container:
        containers["fastp"]
    params:
        opts="-z 6 -l 30 --cut_tail",
    shadow:
        "shallow"
    shell:
        """
        fastp --in1 {input.r1} --in2 {input.r2} --out1 {output.r1} --out2 {output.r2} \
          -w {threads} -h {output.report} 2> {log}
        """


rule map_nanopore_to_decontam_db:
    input:
        index=rules.index_decontam_db_minimap2.output.mm2_index,
        query=rules.combine_fastqs.output.fastq,
    output:
        bam=ont_results / "mapped/{sample}.sorted.bam",
        index=ont_results / "mapped/{sample}.sorted.bam.bai",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(16 * GB),
    params:
        map_extras="-aL2 -x map-ont",
    conda:
        str(env_dir / "aln_tools.yaml")
    log:
        rule_log_dir / "map_to_decontam_db/{sample}.log",
    shell:
        """
        (minimap2 {params.map_extras} -t {threads} {input.index} {input.query} | \
            samtools sort -@ {threads} -o {output.bam}) 2> {log}
        samtools index -@ {threads} {output.bam} &>> {log}
        """


rule map_illumina_to_decontam_db:
    input:
        index=rules.index_decontam_db_bwa.output.bwa_index,
        ref=rules.build_decontamination_db.output.fasta,
        r1=rules.illumina_preprocessing.output.r1,
        r2=rules.illumina_preprocessing.output.r2,
    output:
        bam=illumina_results / "mapped/{source}/{sample}.sorted.bam",
        index=illumina_results / "mapped/{source}/{sample}.sorted.bam.bai",
    threads: 8
    resources:
        mem_mb=lambda wildcards, attempt: attempt * int(12 * GB),
    params:
        map_extras="-M",
    conda:
        str(env_dir / "aln_tools.yaml")
    log:
        rule_log_dir / "map_illumina_to_decontam_db/{source}/{sample}.log",
    shell:
        """
        (bwa mem {params.map_extras} -t {threads} {input.ref} {input.r1} {input.r2} | \
            samtools sort -n -@ {threads} | \
            samtools fixmate -m -@ {threads} - - | \
            samtools sort -@ {threads} | \
            samtools markdup -r -S -O bam - {output.bam}) 2> {log} 
        samtools index -b -@ {threads} {output.bam} 2>> {log}
        """


rule filter_nanopore_contamination:
    input:
        bam=rules.map_nanopore_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        keep_ids=ont_results / "filtered/{sample}/keep.reads",
        contam_ids=ont_results / "filtered/{sample}/contaminant.reads",
        unmapped_ids=ont_results / "filtered/{sample}/unmapped.reads",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    conda:
        str(env_dir / "filter.yaml")
    params:
        script=scripts_dir / "filter_contamination.py",
        extra="--verbose --ignore-secondary",
        outdir=lambda wildcards, output: Path(output.keep_ids).parent,
    log:
        rule_log_dir / "filter_nanopore_contamination/{sample}.log",
    shell:
        """
        python {params.script} {params.extra} \
            -i {input.bam} \
            -m {input.metadata} \
            -o {params.outdir} 2> {log}
        """


rule filter_illumina_contamination:
    input:
        bam=rules.map_illumina_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        keep_ids=illumina_results / "filtered/{source}/{sample}/keep.reads",
        contam_ids=illumina_results / "filtered/{source}/{sample}/contaminant.reads",
        unmapped_ids=illumina_results / "filtered/{source}/{sample}/unmapped.reads",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * 4 * GB,
    conda:
        str(env_dir / "filter.yaml")
    params:
        script=scripts_dir / "filter_contamination.py",
        extra="--verbose --ignore-secondary",
        outdir=lambda wildcards, output: Path(output.keep_ids).parent,
    log:
        rule_log_dir / "filter_illumina_contamination/{source}/{sample}.log",
    shell:
        """
        python {params.script} {params.extra} \
            -i {input.bam} \
            -m {input.metadata} \
            -o {params.outdir} 2> {log}
        """


rule extract_decontaminated_nanopore_reads:
    input:
        reads=rules.map_nanopore_to_decontam_db.input.query,
        read_ids=rules.filter_nanopore_contamination.output.keep_ids,
    output:
        reads=ont_results / "filtered/{sample}/{sample}.filtered.fq.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(4 * GB) * attempt,
    log:
        rule_log_dir / "extract_decontaminated_nanopore_reads/{sample}.log",
    container:
        containers["seqkit"]
    shell:
        "seqkit grep -o {output.reads} -f {input.read_ids} {input.reads} 2> {log}"


rule extract_decontaminated_illumina_reads:
    input:
        r1=rules.map_illumina_to_decontam_db.input.r1,
        r2=rules.map_illumina_to_decontam_db.input.r2,
        read_ids=rules.filter_illumina_contamination.output.keep_ids,
    output:
        r1=illumina_results / "filtered/{source}/{sample}/{sample}_R1.filtered.fq.gz",
        r2=illumina_results / "filtered/{source}/{sample}/{sample}_R2.filtered.fq.gz",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(4 * GB) * attempt,
    log:
        rule_log_dir / "extract_decontaminated_illumina_reads/{source}/{sample}.log",
    container:
        containers["seqkit"]
    shell:
        """
        seqkit grep -o {output.r1} -f {input.read_ids} {input.r1} 2> {log}
        seqkit grep -o {output.r2} -f {input.read_ids} {input.r2} 2>> {log}
        """


rule subsample_nanopore_reads:
    input:
        reads=rules.extract_decontaminated_nanopore_reads.output.reads,
    output:
        reads=ont_results / "subsampled/{sample}/{sample}.subsampled.fq.gz",
    threads: 1
    resources:
        mem_mb=int(0.5 * GB),
    container:
        containers["rasusa"]
    params:
        covg=config["max_covg"]["nanopore"],
        genome_size=config["genome_size"],
        seed=88,
    log:
        rule_log_dir / "subsample_nanopore_reads/{sample}.log",
    shell:
        """
        rasusa -c {params.covg} \
            -g {params.genome_size} \
            -i {input.reads} \
            -o {output.reads} \
            -s {params.seed} 2> {log}
        """


rule subsample_illumina_reads:
    input:
        reads=[
            rules.extract_decontaminated_illumina_reads.output.r1,
            rules.extract_decontaminated_illumina_reads.output.r2,
        ],
    output:
        r1=illumina_results
        / "subsampled/{source}/{sample}/{sample}_R1.subsampled.fq.gz",
        r2=illumina_results
        / "subsampled/{source}/{sample}/{sample}_R2.subsampled.fq.gz",
    threads: 1
    resources:
        mem_mb=int(0.5 * GB),
    container:
        containers["rasusa"]
    params:
        covg=config["max_covg"]["illumina"],
        genome_size=config["genome_size"],
        seed=88,
    log:
        rule_log_dir / "subsample_illumina_reads/{source}/{sample}.log",
    shell:
        """
        rasusa -c {params.covg} \
            -g {params.genome_size} \
            -i {input.reads} \
            -o {output.r1} {output.r2} \
            -s {params.seed} 2> {log}
        """


rule generate_nanopore_krona_input:
    input:
        bam=rules.map_nanopore_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        krona_input=ont_results / "plots/krona/{sample}.krona.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(1 * GB) * attempt,
    conda:
        str(env_dir / "krona.yaml")
    params:
        script=scripts_dir / "generate_krona_input.py",
        extras="--ignore-secondary",
    log:
        rule_log_dir / "generate_nanopore_krona_input/{sample}.log",
    shell:
        """
        python {params.script} {params.extras} \
            -i {input.bam} -m {input.metadata} -o {output.krona_input} 2> {log}
        """


rule generate_illumina_krona_input:
    input:
        bam=rules.map_illumina_to_decontam_db.output.bam,
        metadata=rules.build_decontamination_db.output.metadata,
    output:
        krona_input=illumina_results / "plots/krona/{source}/{sample}.krona.tsv",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: int(1 * GB) * attempt,
    conda:
        str(env_dir / "krona.yaml")
    params:
        script=scripts_dir / "generate_krona_input.py",
        extras="--ignore-secondary",
    log:
        rule_log_dir / "generate_illumina_krona_input/{source}/{sample}.log",
    shell:
        """
        python {params.script} {params.extras} \
            -i {input.bam} -m {input.metadata} -o {output.krona_input} 2> {log}
        """


rule nanopore_summary_stats:
    input:
        reads=rules.subsample_nanopore_reads.input.reads,
    output:
        table=ont_results / "summary/{sample}.nanoq.tsv",
    container:
        containers["nanoq"]
    shell:
        "nanoq -i {input.reads} -Hs 2> {output}"

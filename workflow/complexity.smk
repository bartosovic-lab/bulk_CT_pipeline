rule downscale_bam:
    input:
        bam   = "{out}/{sample}/{sample}_merged.bam",
        index = "{out}/{sample}/{sample}_merged.bam.bai",
    output:
        bam   = temp("{out}/{sample}/complexity/{sample}_down_{fraction}.bam"),
        index = temp("{out}/{sample}/complexity/{sample}_down_{fraction}.bam.bai"),
    conda: "../envs/bulkCT.yaml"
    resources:
        mem_mb = 8000
    shell:
        "samtools view -h -s {wildcards.fraction} {input.bam} | samtools view -bS > {output.bam}; "
        "samtools index {output.bam} -o {output.index}"

rule  picard_duplicates_downscaled:
    input:
        bam   = "{out}/{sample}/complexity/{sample}_down_{fraction}.bam",
    output:
        bam      = temp("{out}/{sample}/complexity/{sample}_down_{fraction}_picard_dups.bam"),
        report   = temp("{out}/{sample}/complexity/{sample}_down_{fraction}_report.txt"),
    resources:
        mem_mb = 32000
    conda: "../envs/bulkCT.yaml"
    shell:
        "picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.report}"

rule report_aggregate:
    input:
        downsampled = lambda wildcards: expand("{out}/{sample}/complexity/{sample}_down_{fraction}_report.txt", out=wildcards.out, sample = wildcards.sample, fraction = subsample_fraction),
        full_bam    = "{out}/{sample}/{sample}_marked_dup_metrics.txt",
    output:
        "{out}/{sample}/complexity/duplicates_table.tsv",
        "{out}/{sample}/complexity/sequencing_saturation.png",
        "{out}/{sample}/complexity/library_size_plot.png",
    params:
        script = workflow.basedir + "/../scripts/get_percent_duplication.R",
        outdir = "{out}/{sample}/complexity/"
    conda: "../envs/bulkCT.yaml"
    shell:
        "Rscript {params.script} -i {input.downsampled} {input.full_bam} -o {params.outdir} "

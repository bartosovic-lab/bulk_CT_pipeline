rule downscale_bam:
    input:
        bam   = temp("{out}/{sample}/{sample}_merged.bam"),
        index = temp("{out}/{sample}/{sample}_merged.bam.bai"),
    output:
        bam   = "{out}/{sample}/{sample}_down_{fraction}.bam",
        index = "{out}/{sample}/{sample}_down_{fraction}.bam.bai",
    conda: "../envs/bulkCT.yaml"
    shell:
        "samtools view -h -s {wildcards.fraction} {input.bam} | samtools view -bS > {output.bam}; "
        "samtools index {output.bam} -o {output.index}"

rule  picard_duplicates_downscaled:
    input:
        bam   = "{out}/{sample}/{sample}_down_{fraction}.bam",
    output:
        bam      = temp("{out}/{sample}/{sample}_down_{fraction}_picard_dups.bam"),
        report   = temp("{out}/{sample}/{sample}_down_{fraction}_report.txt")   ,
    conda: "../envs/bulkCT.yaml"
    shell:
        "picard MarkDuplicates -I {input.bam} -O {output.bam} -M {output.report}"

rule report_aggregate:
    input:
        downsampled = lambda wildcards: expand("{out}/{sample}/{sample}_down_{fraction}_report.txt", out=wildcards.out, sample = wildcards.sample, fraction = subsample_fraction),
        full_bam    = "{out}/{sample}/{sample}_marked_dup_metrics.txt",
    output:
        "{out}/{sample}/downsample_report.txt"
    params:
        script = workflow.basedir + "/../scripts/get_percent_duplication.R"
    shell:
        "Rscript {params.script} {input.downsampled} {input.full_bam} > {output}"

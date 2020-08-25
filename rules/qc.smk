# =================================================================================================
#     Basic QC Stats
# =================================================================================================

rule fastqc:
    input:
        unpack(get_fastq)
    output:
        html="qc/fastqc/{sample}-{unit}.html",
        zip="qc/fastqc/{sample}-{unit}.zip"
    log:
        "logs/fastqc/{sample}-{unit}.log"
    benchmark:
        "benchmarks/fastqc/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.27.1/bio/fastqc"

rule samtools_stats:
    input:
        get_mapping_result()
    output:
        "qc/samtools-stats/{sample}-{unit}.txt"
    log:
        "logs/samtools-stats/{sample}-{unit}.log"
    benchmark:
        "benchmarks/samtools-stats/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.27.1/bio/samtools/stats"

rule samtools_flagstat:
    input:
        get_mapping_result()
    output:
        "qc/samtools-flagstats/{sample}-{unit}.txt"
    log:
        "logs/samtools-flagstats/{sample}-{unit}.log"
    benchmark:
        "benchmarks/samtools-flagstats/{sample}-{unit}.bench.log"
    group:
        "qc"
    wrapper:
        "0.64.0/bio/samtools/flagstat"

rule qualimap:
    input:
        get_mapping_result()
    output:
        "qc/qualimap/{sample}-{unit}/genome_results.txt",
        "qc/qualimap/{sample}-{unit}/qualimapReport.html",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/coverage_histogram.txt",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/genome_fraction_coverage.txt",
        "qc/qualimap/{sample}-{unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt"

        # Alternatively, specify just the out dir.
        # But that gives us less control over whether the rule executed succeeded.
        # outdir=directory("qc/qualimap/{sample}-{unit}")
    params:
        extra=config["params"]["qualimap"]["extra"],
        outdir="qc/qualimap/{sample}-{unit}"
    threads:
        config["params"]["qualimap"]["threads"]
    log:
        stderr="logs/qualimap/{sample}-{unit}_qualimap.stderr",
        stdout="logs/qualimap/{sample}-{unit}_qualimap.stdout"
    group:
        "qc"
    conda:
        "../envs/qualimap.yaml"
    shell:
        "unset DISPLAY; qualimap bamqc -bam {input} -nt {threads} "
        "-outdir {params.outdir} -outformat HTML "
        "{params.extra} > {log.stdout} 2> {log.stderr}"

# =================================================================================================
#     MultiQC
# =================================================================================================

# Different trimming tools produce different summary files. This function simply returns these file
# names as strings, without replacing the wildcards. Then, when the function is called below,
# these are expanded.
def get_trimming_report():
    result=[]
    for smp in samples.itertuples():
        # Get the corret suffix for single and paired end fastq files.
        # We need this, as otherwise the trimming rules have colliding output.
        suffix = "se" if is_single_end(smp.sample, smp.unit) else "pe"

        # Now append the file for the sample to the result list
        if config["settings"]["trimming-tool"] == "adapterremoval":
            result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + ".settings" )
        elif config["settings"]["trimming-tool"] == "cutadapt":
            result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".qc-" + suffix + ".txt" )
        elif config["settings"]["trimming-tool"] == "skewer":
            result.append( "trimmed/" + smp.sample + "-" + smp.unit + "-" + suffix + "-trimmed.log" )
        elif config["settings"]["trimming-tool"] == "trimmomatic":
            result.append( "trimmed/" + smp.sample + "-" + smp.unit + ".trimlog.log" )
        else:
            raise Exception("Unknown trimming-tool: " + config["settings"]["trimming-tool"])
    return result

# Different dedup tools produce different summary files. See above for details.
def get_dedup_report():
    # Switch to the chosen duplicate marker tool
    if config["settings"]["duplicates-tool"] == "picard":
        return "qc/dedup/{u.sample}-{u.unit}.metrics.txt"
    elif config["settings"]["duplicates-tool"] == "dedup":
        return "dedup/{u.sample}-{u.unit}.sorted.dedup.json"
    else:
        raise Exception("Unknown duplicates-tool: " + config["settings"]["duplicates-tool"])

# Unfortunately, in some environments, multiqc does not work due to char encoding issues, see
# https://github.com/ewels/MultiQC/issues/484 ... If you run into this issue, try running it locally.
rule multiqc:
    input:
        # Trimming
        get_trimming_report(),

        # Dedup
        expand(get_dedup_report(), u=samples.itertuples()),

        # QC tools
        expand("qc/fastqc/{u.sample}-{u.unit}.zip", u=samples.itertuples()),
        expand("qc/samtools-stats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),
        expand("qc/samtools-flagstats/{u.sample}-{u.unit}.txt", u=samples.itertuples()),

        # Qualimap
        # expand("qc/qualimap/{u.sample}-{u.unit}", u=samples.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/genome_results.txt", u=samples.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/qualimapReport.html", u=samples.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/coverage_histogram.txt", u=samples.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/genome_fraction_coverage.txt", u=samples.itertuples()),
        expand("qc/qualimap/{u.sample}-{u.unit}/raw_data_qualimapReport/mapped_reads_gc-content_distribution.txt", u=samples.itertuples()),

        # Annotation
        "snpeff/all.csv" if config["settings"]["snpeff"] else []
    output:
        report("qc/multiqc.html", caption="../reports/multiqc.rst", category="Quality control")
    params:
        config["params"]["multiqc"]["extra"],
    log:
        "logs/multiqc.log"
    # conda:
    #     "../envs/multiqc.yaml"
    wrapper:
        "0.64.0/bio/multiqc"

# Rule is not submitted as a job to the cluster.
localrules: multiqc

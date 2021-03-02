# =================================================================================================
#     mapDamage
# =================================================================================================

rule mapdamage:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads
        # "mapped/{sample}-{unit}.sorted.bam"
    output:
        "mapdamage/{sample}-{unit}/Runtime_log.txt"
    params:
        index=config["data"]["reference"]["genome"],
        extra=config["params"]["mapdamage"]["extra"],
        outdir="mapdamage/{sample}-{unit}"
    log:
        "logs/mapdamage/{sample}-{unit}.log"
    conda:
        "../envs/mapdamage.yaml"
    shell:
        "mapDamage -i {input} -r {params.index} -d {params.outdir} {params.extra} > {log} 2>&1"

# =================================================================================================
#     DamageProfiler
# =================================================================================================

rule damageprofiler:
    input:
        # Get either the normal mapping output, or, if additional filtering via `samtools view`
        # is set in the config settings: filter-mapped-reads, use the filtered output instead.
        get_mapped_reads
        # "mapped/{sample}-{unit}.sorted.bam"
    output:
        "damageprofiler/{sample}-{unit}/DamageProfiler.log"
    params:
        index=config["data"]["reference"]["genome"],
        extra=config["params"]["damageprofiler"]["extra"],
        outdir="damageprofiler/{sample}-{unit}"
    log:
        "logs/damageprofiler/{sample}-{unit}.log"
    conda:
        "../envs/damageprofiler.yaml"
    shell:
        "damageprofiler -i {input} -r {params.index} -o {params.outdir} {params.extra} > {log} 2>&1"
        # "java -jar DamageProfiler-0.5.0.jar "

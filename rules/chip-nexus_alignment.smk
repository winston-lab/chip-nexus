#!/usr/bin/env python

localrules: bowtie2_build,
    index_bam

rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    log: "logs/bowtie2_build-{basename}.log"
    shell: """
        (bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.bt2", num=[1,2]),
        fastq = f"fastq/cleaned/{{sample}}_{FACTOR}-chipnexus-clean.fastq.gz"
    output:
        bam = temp(f"alignment/{{sample}}_{FACTOR}-chipnexus-uniquemappers.bam"),
        unaligned_fastq = f"fastq/{{sample}}_{FACTOR}-chipnexus-unaligned.fastq.gz",
        log = "logs/align/align_{sample}.log"
    params:
        outbase =  config["bowtie2"]["index-path"] + "/" +  config["combinedgenome"]["name"],
        minmapq = config["bowtie2"]["minmapq"]
    threads : config["threads"]
    shell: """
        (bowtie2 -x {params.outbase} -U {input.fastq} --un-gz {output.unaligned_fastq} -p {threads} | samtools view -buh -q {params.minmapq} - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        """

rule remove_PCR_duplicates:
    input:
        f"alignment/{{sample}}_{FACTOR}-chipnexus-uniquemappers.bam"
    output:
        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam"
    log: "logs/remove_PCR_duplicates/remove_PCR_duplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        """

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam"
    output:
        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input}) &> {log}
        """

#peakcalling programs (MACS2 and Qnexus) require BAM input
rule bam_separate_species:
    input:
        bam = f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam",
        bai = f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        f"alignment/{{sample}}_{FACTOR}-chipnexus-noPCRduplicates-{{species}}.bam",
    params:
        filterprefix = lambda wc: config["combinedgenome"]["spikein_prefix"] if wc.species=="experimental" else config["combinedgenome"]["experimental_prefix"],
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.species=="experimental" else config["combinedgenome"]["spikein_prefix"]
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {params.prefix} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{params.prefix}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """


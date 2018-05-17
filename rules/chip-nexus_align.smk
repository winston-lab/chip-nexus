#!/usr/bin/env python

rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.bt2", num=[1,2]),
        fastq = expand("fastq/cleaned/{{sample}}_{factor}-chipnexus-clean.fastq.gz", factor=FACTOR)
    output:
        bam = temp(expand("alignment/{{sample}}_{factor}-chipnexus-uniquemappers.bam", factor=FACTOR)),
        unaligned_fq = expand("fastq/unaligned/{{sample}}_{factor}-chipnexus-unaligned.fastq.gz", factor=FACTOR),
        log = "logs/align/align_{sample}.log"
    params:
        outbase =  config["bowtie2"]["index-path"] + "/" +  config["combinedgenome"]["name"],
        minmapq = config["bowtie2"]["minmapq"]
    threads : config["threads"]
    shell: """
        (bowtie2 -x {params.outbase} -U {input.fastq} --un-gz {output.unaligned_fq} -p {threads} | samtools view -buh -q {params.minmapq} - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
        """

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}_{factor}-chipnexus-uniquemappers.bam"
    output:
        "alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam"
    log: "logs/remove_PCR_duplicates/remove_PCR_duplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        """

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        "alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam"
    output:
        "alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input.bam}) &> {log}
        """

#peakcalling programs (MACS2 and Qnexus) require BAM input
rule bam_separate_species:
    input:
        bam = "alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam",
        bai = "alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    params:
        filterprefix = lambda wc: config["combinedgenome"]["spikein_prefix"] if wc.species==config["combinedgenome"]["experimental_prefix"] else config["combinedgenome"]["experimental_prefix"],
    output:
        "alignment/{sample}-{species}only.bam"
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{wildcards.species}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """


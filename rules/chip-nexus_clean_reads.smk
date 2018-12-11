#!/usr/bin/env python

rule clean_reads:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    output:
        fastq = temp(f"fastq/cleaned/{{sample}}_{FACTOR}-chipnexus-noadapter.fastq.gz"),
        log =  "logs/remove_adapter/remove_adapter-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    conda:
        "../envs/cutadapt.yaml"
    threads:
        config["threads"]
    shell: """
        cutadapt --front=^NNNNNCTGA --error-rate=0 --no-trim --discard-untrimmed --cores={threads} {input} | cutadapt --adapter={params.adapter} --error-rate=0.1 --nextseq-trim={params.trim_qual} --minimum-length=11 --cores={threads} -o {output.fastq} - &> {output.log}
       """

rule extract_molecular_barcode:
    input:
        expand("fastq/cleaned/{{sample}}_{factor}-chipnexus-noadapter.fastq.gz", factor=FACTOR)
    output:
        fastq = expand("fastq/cleaned/{{sample}}_{factor}-chipnexus-clean.fastq.gz", factor=FACTOR),
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    log: "logs/extract_molecular_barcode/extract_molecular_barcode-{sample}.log"
    shell: """
        (python scripts/extract_nexus_molecular_barcode.py {input} {output.fastq} {output.barcodes} {output.ligation}) &> {log}
        """


#!/usr/bin/env python

localrules: combine_peaks

rule callpeaks_macs2:
    input:
        bam = lambda wc: expand("alignment/{sample}_{factor}-chipnexus-noPCRduplicates-{species}.bam", sample=[k for k,v in PASSING.items() if v["group"]==wc.group], factor=FACTOR, species=wc.species),
        fasta = lambda wc: os.path.abspath(build_annotations(config["genome"]["fasta"])) if wc.species=="experimental" else os.path.abspath(config["spike_in"]["fasta"]),
    output:
        tsv = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.xls",
        peaks = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.narrowPeak",
        summits = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_summits.bed",
        script = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_model.r",
        pdf = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_model.pdf",
        treat_bg = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_treat_pileup.bdg",
        cntrl_bg = "peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_control_lambda.bdg"
    params:
        bw = config["macs2"]["bw"],
        llocal = config["macs2"]["llocal"],
        qscore = config["macs2"]["fdr"]
    conda:
        "../envs/macs2.yaml"
    log: "logs/macs2/macs2_{group}-{species}-{factor}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(faidx {input.fasta} -i chromsizes | awk '{{sum += $2}} END {{print sum}}') --keep-dup all --bdg -n peakcalling/macs/{wildcards.group}/{wildcards.group}_{wildcards.species}-{wildcards.factor}-chipnexus --SPMR --bw {params.bw} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript {output.script}) &>> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.peaks}) &>> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.summits}) &>> {log}
        """

rule combine_peaks:
    input:
        cond = "peakcalling/macs/{condition}/{condition}_{type}-{factor}-chipnexus_peaks.narrowPeak",
        ctrl = "peakcalling/macs/{control}/{control}_{type}-{factor}-chipnexus_peaks.narrowPeak",
    output:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}_{type}-{factor}-peaks.bed"
    log: "logs/combine_peaks/combine_peaks-{condition}-v-{control}-{type}-{factor}.log"
    shell: """
        (bedtools multiinter -i {input} | bedtools merge -i stdin | sort -k1,1 -k2,2n > {output}) &> {log}
        """


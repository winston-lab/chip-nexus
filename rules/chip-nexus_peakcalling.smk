#!/usr/bin/env python

rule callpeaks_macs2:
    input:
        bam = lambda wc: expand("alignment/{sample}_{factor}-chipnexus-noPCRduplicates-{species}.bam", sample=[k for k,v in PASSING.items() if v["group"]==wc.group], factor=FACTOR, species=wc.species),
        chrsizes = lambda wc: config["genome"]["chrsizes"] if wc.species==config["combinedgenome"]["experimental_prefix"] else config["genome"]["si-chrsizes"],
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
    log: "logs/macs2/macs2-{group}-{species}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --keep-dup all --bdg -n peakcalling/macs/{wildcards.group}/{wildcards.group}_{wildcards.species}-{wildcards.factor}-chipnexus --SPMR --bw {params.bw} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript {output.script}) &>> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.peaks}) &>> {log}
        (sed -i -e 's/peakcalling\/macs\/{wildcards.group}\///g' {output.summits}) &>> {log}
        """

#!/usr/bin/env python

localrules:
    combine_peaks

rule callpeaks_macs2:
    input:
        bam = "alignment/{sample}_{factor}-chipnexus-noPCRduplicates-{species}.bam",
        fasta = lambda wc: os.path.abspath(build_annotations(config["genome"]["fasta"])) if wc.species=="experimental" else os.path.abspath(config["spike_in"]["fasta"]),
    output:
        tsv = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_peaks.xls",
        peaks = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_peaks.narrowPeak",
        summits = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_summits.bed",
        treat_bg = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_treat_pileup.bdg",
        cntrl_bg = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_control_lambda.bdg"
    params:
        script = "peakcalling/sample_peaks/{sample}_{species}-{factor}-chipnexus_model.r",
        bw = config["peakcalling"]["bandwidth"],
        llocal = config["peakcalling"]["l_local"],
        build_model = "" if config["peakcalling"]["build_model"] else "--nomodel",
        mfold_low = config["peakcalling"]["mfold_low"],
        mfold_high = config["peakcalling"]["mfold_high"],
        shift = config["peakcalling"]["shift"] if not config["peakcalling"]["build_model"] else 0,
        extsize = config["peakcalling"]["extsize"] if not config["peakcalling"]["build_model"] else 0
    conda:
        "../envs/macs2.yaml"
    log:
        "logs/macs2/macs2_{sample}-{species}-{factor}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(faidx {input.fasta} -i chromsizes | \
                awk '{{sum += $2}} END {{print sum}}') \
         --keep-dup all --bdg -n peakcalling/sample_peaks/{wildcards.sample}_{wildcards.species}-{wildcards.factor}-chipnexus --SPMR --bw {params.bw} {params.build_model} --extsize {params.extsize} --shift {params.shift} --mfold {params.mfold_low} {params.mfold_high} --llocal {params.llocal} --call-summits -q 1) &> {log}
        if [ {params.build_model} = "" ]
        then
          (Rscript {params.script}) &>> {log}
        fi
        (sed -i -e 's/peakcalling\/sample_peaks\///g' {output.peaks}) &>> {log}
        (sed -i -e 's/peakcalling\/sample_peaks\///g' {output.summits}) &>> {log}
        """

rule idr:
    input:
        #NOTE: for now we take the first two samples since the IDR script only takes two
        #change this if we find a better way to aggregate results
        lambda wc: ["peakcalling/sample_peaks/" + x + "_{species}-{factor}-chipnexus_peaks.narrowPeak".format(**wc) for x in PASSING if PASSING[x]['group']==wc.group][0:2]
    output:
        allpeaks = "peakcalling/{group}/{group}_{species}-{factor}-chipnexus-idrpeaks-all.tsv",
        filtered = "peakcalling/{group}/{group}_{species}-{factor}-chipnexus-idrpeaks-filtered.tsv",
        narrowpeak = "peakcalling/{group}/{group}_{species}-{factor}-chipnexus-idrpeaks.narrowPeak",
        summits = "peakcalling/{group}/{group}_{species}-{factor}-chipnexus-idrpeaks-summits.bed",
    params:
        idr = int(-125*log2(config["peakcalling"]["idr"]))
    conda:
        "../envs/idr.yaml"
    log:
        "logs/idr/idr-{group}-{species}-{factor}.log"
    shell: """
        (idr -s {input} --input-file-type narrowPeak --rank q.value -o {output.allpeaks} -l {log} --plot --peak-merge-method max) &> {log}
        (awk '$5>{params.idr} || $9=="inf"' {output.allpeaks} | \
         LC_COLLATE=C sort -k1,1 -k2,2n | \
         tee {output.filtered} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4, $5, $6, $7, $11, $12, $10}}' | \
         tee {output.narrowpeak} | \
         awk 'BEGIN{{FS=OFS="\t"}}{{start=$2+$10; print $1, start, start+1, $4, $5, $6}}' > {output.summits}) &>> {log}
        """

rule combine_peaks:
    input:
        cond = "peakcalling/{condition}/{condition}_{species}-{factor}-chipnexus-idrpeaks.narrowPeak",
        ctrl = "peakcalling/{control}/{control}_{species}-{factor}-chipnexus-idrpeaks.narrowPeak",
    output:
        "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed"
    log:
        "logs/combine_peaks/combine_peaks-{condition}-v-{control}-{species}-{factor}.log"
    shell: """
        (bedtools multiinter -i {input} | \
         bedtools merge -i stdin | \
         awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, ".", 0, "."}}' | \
         sort -k1,1 -k2,2n > {output}) &> {log}
        """


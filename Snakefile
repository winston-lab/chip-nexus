#!/usr/bin/env python
import os
from math import log2

configfile: "config.yaml"

SAMPLES = config["samples"]
# name_string = " ".join(SAMPLES)
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
# pass_string = " ".join(PASSING)
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]
controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

# CATEGORIES = ["genic", "intragenic", "intergenic", "antisense", "convergent", "divergent"]

localrules: all,
            make_stranded_annotations,
            join_window_counts,
            bedgraph_to_bigwig,
            gzip_deeptools_table,
            # melt_matrix,
            cat_matrices

rule all:
    input:
        #FastQC
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        #alignment
        expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES),
        #coverage
        expand("coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-minus.bedgraph", sample=SAMPLES, factor = config["factor"]),
        # expand("coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-STRANDED.bedgraph",sample=SAMPLES, factor = config["factor"], norm=["counts"])
        # expand("coverage/counts/{sample}-{factor}-windowcounts.bedgraph", sample=SAMPLES, factor = config["factor"]),
        # "coverage/counts/union-bedgraph-tfiib-chipnexus-allwindowcounts.tsv",
        #initial QC
        # "qual_ctrl/all/all-pca-scree-libsizenorm.png",
        #datavis
        # expand("datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-heatmap-bygroup.png", annotation = config["annotations"], norm = ["libsizenorm", "spikenorm"], factor=config["factor"]),
        # expand("datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bygroup.png", annotation = config["annotations"], norm = ["libsizenorm", "spikenorm"], factor=config["factor"]),
        #qnexus
        expand("peakcalling/qnexus/{sample}-{factor}-Q-treatment.bedgraph", sample=SAMPLES, factor=config["factor"]),
        #macs2
        expand("peakcalling/macs/{group}_peaks.narrowPeak", group = GROUPS),
        expand("peakcalling/macs/{group}_peaks-intragenic.narrowPeak", group = GROUPS),
        "peakcalling/macs/intragenic-peaks-overlap.tsv"

rule fastqc_raw:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        "qual_ctrl/fastqc/raw/{sample}"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        mkdir -p {output}
        (fastqc -o {output} --noextract -t {threads} {input}) &> {log}
        """

rule remove_adapter:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        temp("fastq/cleaned/{sample}-noadapter.fastq")
    params:
        adapter = config["cutadapt"]["adapter"]
    log: "logs/remove_adapter/remove_adapter-{sample}.log"
    shell: """
       (less {input} | grep -B 1 -A 2 --no-group-separator '^.....CTGA' | cutadapt -e 0.1 -a {params.adapter} -m 11 -o {output} -) &> {log}
       """

rule remove_molec_barcode:
    input:
        "fastq/cleaned/{sample}-noadapter.fastq"
    output:
        fq = "fastq/cleaned/{sample}-clean.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/removeMBC-{sample}.log"
    shell: """
        (python scripts/extractNexusMolecularBarcode.py {input} fastq/cleaned/{wildcards.sample}-clean.fastq {output.barcodes} {output.ligation}) &> {log}
        pigz -f fastq/cleaned/{wildcards.sample}-clean.fastq
        """

rule fastqc_cleaned:
    input:
        "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        html = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.html",
        folder  = "qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip"
    threads : config["threads"]
    log: "logs/fastqc/cleaned/fastqc-cleaned-{sample}.log"
    shell: """
        mkdir -p qual_ctrl/fastqc/cleaned/{wildcards.sample}
        (fastqc -o qual_ctrl/fastqc/cleaned/{wildcards.sample} --noextract -t {threads} {input}) &> {log}
        """

rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2])
    params:
        name = config["combinedgenome"]["name"]
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} ../genome/bowtie2_indexes/{params.name}) &> {log}
        """

rule align:
    input:
        expand("../genome/bowtie2_indexes/{basename}.{num}.bt2", basename=config["combinedgenome"]["name"], num = [1,2,3,4]),
        expand("../genome/bowtie2_indexes/{basename}.rev.{num}.bt2", basename=config["combinedgenome"]["name"], num=[1,2]),
        fastq = "fastq/cleaned/{sample}-clean.fastq.gz"
    output:
        bam = temp("alignment/{sample}-unique.bam"),
        aligned_fq = "fastq/aligned/{sample}-aligned.fastq.gz",
        unaligned_fq = "fastq/unaligned/{sample}-unaligned.fastq.gz"
    params:
        outbase = "../genome/bowtie2_indexes/" + config["combinedgenome"]["name"],
        minmapq = config["bowtie2"]["minmapq"]
    threads : config["threads"]
    log: "logs/align/align-{sample}.log"
    shell: """
        (bowtie2 -x {params.outbase} -U {input.fastq} --al-gz {output.aligned_fq} --un-gz {output.unaligned_fq} -p {threads} | samtools view -buh -q {params.minmapq} - | samtools sort -T {wildcards.sample} -@ {threads} -o {output.bam} -) &> {log}
        """

rule remove_PCR_duplicates:
    input:
        "alignment/{sample}-unique.bam"
    output:
        "alignment/{sample}-noPCRdup.bam"
    log: "logs/remove_PCR_duplicates/removePCRduplicates-{sample}.log"
    shell: """
        (python scripts/removePCRdupsFromBAM.py {input} {output}) &> {log}
        """

rule qnexus:
    input:
        "alignment/{sample}-noPCRdup.bam"
    output:
        "peakcalling/qnexus/{sample}-{factor}-Q-narrowPeak.bed",
        "peakcalling/qnexus/{sample}-{factor}-Q-qfrag-binding-characteristics.R",
        "peakcalling/qnexus/{sample}-{factor}-Q-qfrag-binding-characteristics.pdf",
        "peakcalling/qnexus/{sample}-{factor}-Q-quality-statistics.tab",
        "peakcalling/qnexus/{sample}-{factor}-Q-runinfo.txt",
        "peakcalling/qnexus/{sample}-{factor}-Q-summit-info.tab",
        "peakcalling/qnexus/{sample}-{factor}-Q-treatment.bedgraph",
        coverage = "coverage/counts/{sample}-{factor}-chipnexus-qfragcounts.bedgraph"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
    log: "logs/qnexus/qnexus-{sample}-{factor}.log"
    shell: """
        (../programs/Q/bin/Q --nexus-mode -t {input} -o peakcalling/qnexus/{wildcards.sample}-{wildcards.factor} -v -wbt) &> {log}
        # (sed -i 's/peakcalling\///g' peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-qfrag-binding-characteristics.R) &>> {log}
        (Rscript peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-qfrag-binding-characteristics.R) &>> {log}
        (grep {params.exp_prefix} peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-treatment.bedgraph | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.coverage}) &>> {log}
        """

# rule get_qfrag_counts:
    # input:
        # "peakcalling/qnexus/{sample}-Q-treatment.bedgraph"
    # output:
    #     "coverage/counts/{sample}-{factor}-chipnexus-qfragcounts.bedgraph"
    # params:
    #     exp_prefix = config["combinedgenome"]["experimental_prefix"],
    # log: "get_qfrag_counts/get_qfrag_counts-{sample}-{factor}.log"
    # shell: """
        # grep {params.exp_prefix} {input} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output}
        # """

rule build_macs2_input:
    input:
        bam = "alignment/{sample}-noPCRdup.bam",
        chrsizes = config["combinedgenome"]["chrsizes"],
    output:
        "alignment/{sample}-exponly.bam"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"]
    log: "logs/build_macs2_input/build_macs2_input-{sample}.log"
    shell: """
        (samtools index {input.bam}) &> {log}
        (samtools view -b {input.bam} $(grep {params.exp_prefix} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') > {output}) &>> {log}
        """

#TODO: put parameters in config file
rule macs2:
    input:
        bam = lambda wildcards: expand("alignment/{sample}-exponly.bam", sample= {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.group and v["pass-qc"] == "pass")}),
        chrsizes = config["genome"]["chrsizes"]
    output:
        xls = "peakcalling/macs/{group}_peaks.xls",
        peaks = "peakcalling/macs/{group}_peaks.narrowPeak",
        summits = "peakcalling/macs/{group}_summits.bed",
        script = "peakcalling/macs/{group}_model.r",
        pdf = "peakcalling/macs/{group}_model.pdf",
        treat_bg = "peakcalling/macs/{group}_treat_pileup.bdg",
        cntrl_bg = "peakcalling/macs/{group}_control_lambda.bdg"
    params:
        bw = 160,
        slocal = 200,
        llocal = 1000,
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        qscore = 0.01
    conda:
        "envs/python2.yaml"
    log: "logs/macs2/macs2-{group}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --keep-dup all --bdg -n peakcalling/macs/{wildcards.group} --SPMR --bw {params.bw} --slocal {params.slocal} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript peakcalling/macs/{wildcards.group}_model.r) &>> {log}
        (sed -i -e 's/{params.exp_prefix}//g; s/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}_peaks.narrowPeak) &>> {log}
        (sed -i -e 's/{params.exp_prefix}//g; s/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}_summits.bed) &>> {log}
        (sed -i -e 's/{params.exp_prefix}//g' peakcalling/macs/{wildcards.group}_treat_pileup.bdg) &>> {log}
        (sed -i -e 's/{params.exp_prefix}//g' peakcalling/macs/{wildcards.group}_control_lambda.bdg) &>> {log}
        """

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    output:
        SIplmin = temp("coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-plmin.bedgraph"),
        SIpl = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-plus.bedgraph",
        SImin = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-minus.bedgraph",
        plmin = temp("coverage/counts/{sample}-{factor}-chipnexus-counts-plmin.bedgraph"),
        plus = "coverage/counts/{sample}-{factor}-chipnexus-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-{factor}-chipnexus-counts-minus.bedgraph"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
        si_prefix = config["combinedgenome"]["spikein_prefix"]
    log: "logs/get_coverage/get_coverage-{sample}.log"
    shell: """
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIplmin}) &> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SIpl}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.si_prefix} | sed 's/{params.si_prefix}//g' | sort -k1,1 -k2,2n > {output.SImin}) &>> {log}
        (genomeCoverageBed -bga -5 -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &>> {log}
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.exp_prefix} | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """

rule normalize:
    input:
        plus = "coverage/counts/{sample}-{factor}-chipnexus-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-{factor}-chipnexus-counts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-{factor}-chipnexus-counts-plmin.bedgraph",
        SIplmin = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-plmin.bedgraph",
        qfrags = "coverage/counts/{sample}-{factor}-chipnexus-qfragcounts.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}-{factor}-chipnexus-spikenorm-plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}-{factor}-chipnexus-spikenorm-minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-minus.bedgraph",
        qfrag_spikenorm = "coverage/spikenorm/{sample}-{factor}-chipnexus-qfrags-spikenorm.bedgraph",
        qfrag_libnorm= "coverage/libsizenorm/{sample}-{factor}-chipnexus-qfrags-libsizenorm.bedgraph",
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (scripts/libsizenorm.awk {input.SIplmin} {input.plus} > {output.spikePlus}) &> {log}
        (scripts/libsizenorm.awk {input.SIplmin} {input.minus} > {output.spikeMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.plus} > {output.libnormPlus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.minus} > {output.libnormMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.SIplmin} {input.qfrags} > {output.qfrag_spikenorm}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.qfrags} > {output.qfrag_libnorm}) &>> {log}
        """

rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
        si = config["genome"]["si-chrsizes"]
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        si= os.path.splitext(config["genome"]["si-chrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.si} > {output.si}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-minus.bedgraph"
    output:
        sense = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus}> {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus}> {output.antisense}) &>> {log}
        """

rule make_stranded_sicounts_bedgraph:
    input:
        plus = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-plus.bedgraph",
        minus = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-minus.bedgraph"
    output:
        sense = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-SENSE.bedgraph"
    log: "logs/make_stranded_sicounts_bedgraph/make_stranded_sicounts_bedgraph-{sample}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus}> {output.sense}) &> {log}
        """

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "../genome/annotations/stranded/{annotation}-STRANDED.bed"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule map_counts_to_windows:
    input:
        bg = "coverage/counts/{sample}-{factor}-chipnexus-counts-STRANDED.bedgraph",
        si_bg = "coverage/counts/spikein/{sample}-{factor}-chipnexus-SI-counts-SENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        si_chrsizes = os.path.splitext(config["genome"]["si-chrsizes"])[0] + "-STRANDED.tsv"
    output:
        exp = temp("coverage/counts/{sample}-{factor}-windowcounts.bedgraph"),
        si= temp("coverage/counts/spikein/{sample}-{factor}-SI-windowcounts.bedgraph")
    params:
        windowsize = config["corr-windowsize"]
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {params.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        bedtools makewindows -g {input.si_chrsizes} -w {params.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.si_bg} -c 4 -o sum > {output.si}
        """

rule join_window_counts:
    input:
        exp = expand("coverage/counts/{sample}-{factor}-windowcounts.bedgraph", sample=SAMPLES, factor=config["factor"]),
        si= expand("coverage/counts/spikein/{sample}-{factor}-SI-windowcounts.bedgraph", sample=SAMPLES, factor=config["factor"])
    output:
        exp = "coverage/counts/union-bedgraph-{factor}-chipnexus-allwindowcounts.tsv",
        si = "coverage/counts/spikein/union-bedgraph-{factor}-chipnexus-SI-allwindowcounts.tsv",
    params:
        names = list(SAMPLES.keys())
    shell: """
        bedtools unionbedg -i {input.exp} -names {params.names} -header > .union-bedgraph-exp.temp
        cut -f1-3 .union-bedgraph-exp.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}' > .positions.txt
        cut -f4- .union-bedgraph-exp.temp > .values.txt
        paste .positions.txt .values.txt > {output.exp}
        rm .union-bedgraph-exp.temp .positions.txt .values.txt

        bedtools unionbedg -i {input.si} -names {params.names} -header > .union-bedgraph-si.temp
        cut -f1-3 .union-bedgraph-si.temp | awk 'BEGIN{{FS="\t"; OFS="-"}}{{print $1, $2, $3}}' > .positions.txt
        cut -f4- .union-bedgraph-si.temp > .values.txt
        paste .positions.txt .values.txt > {output.si}
        rm .union-bedgraph-si.temp .positions.txt .values.txt
        """

rule deseq_initial_qc:
    input:
        exp = expand("coverage/counts/union-bedgraph-{factor}-chipnexus-allwindowcounts.tsv", factor = config["factor"]),
        si = expand("coverage/counts/spikein/union-bedgraph-{factor}-chipnexus-SI-allwindowcounts.tsv", factor = config["factor"])
    params:
        # alpha = config["deseq"]["fdr"],
        alpha = 0.1,
        samples = list(SAMPLES.keys()),
        samplegroups = [SAMPLES[x]["group"] for x in SAMPLES],
        nospikein = list({k:v for (k,v) in SAMPLES.items() if (v["spikein"] == "n")}.keys())
    output:
        exp_size_v_sf = protected("qual_ctrl/all/all-libsize-v-sizefactor-experimental.png"),
        si_size_v_sf = protected("qual_ctrl/all/all-libsize-v-sizefactor-spikein.png"),
        si_pct = protected("qual_ctrl/all/all-spikein-pct.png"),
        corrplot_spikenorm = protected("qual_ctrl/all/all-pairwise-correlation-spikenorm.png"),
        corrplot_libsizenorm = protected("qual_ctrl/all/all-pairwise-correlation-libsizenorm.png"),
        #count_heatmap_spikenorm = protected("qual_ctrl/all/all-de-bases-heatmap-spikenorm.png"),
        #count_heatmap_libsizenorm = protected("qual_ctrl/all/all-de-bases-heatmap-libsizenorm.png"),
        dist_heatmap_spikenorm = protected("qual_ctrl/all/all-sample-dists-spikenorm.png"),
        dist_heatmap_libsizenorm = protected("qual_ctrl/all/all-sample-dists-libsizenorm.png"),
        pca_spikenorm = protected("qual_ctrl/all/all-pca-spikenorm.png"),
        scree_spikenorm = protected("qual_ctrl/all/all-pca-scree-spikenorm.png"),
        pca_libsizenorm = protected("qual_ctrl/all/all-pca-libsizenorm.png"),
        scree_libsizenorm = protected("qual_ctrl/all/all-pca-scree-libsizenorm.png")
    script:
       "scripts/deseq2_qc.R"

#sum reads from both strands for visualization purposes
rule make_combined_bedgraph:
    input:
        bg = expand("coverage/{{norm}}/{{sample}}-{{factor}}-chipnexus-{{norm}}-{strand}.bedgraph", strand=["plus","minus"]),
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-COMBINED.bedgraph"
    shell: """
        bedtools unionbedg -i {input.bg} -g {input.chrsizes} | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4+$5}}' > {output}
        """

rule bedgraph_to_bigwig:
    input:
        # bg = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-COMBINED.bedgraph",
        bg = "coverage/{norm}/{sample}-{factor}-chipnexus-qfrags-{norm}.bedgraph",
        # bg = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-{strand}.bedgraph",
        chrsizes = config["genome"]["chrsizes"]
    output:
        "coverage/{norm}/bw/{sample}-{factor}-chipnexus-qfrags-{norm}.bw"
    shell: """
        bedGraphToBigWig {input.bg} {input.chrsizes} {output}
        """

rule bedgraph_to_bigwig_stranded:
    input:
        # bg = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-COMBINED.bedgraph",
        sense = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-ANTISENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        sense = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-SENSE.bw",
        antisense = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-ANTISENSE.bw"
    shell: """
        bedGraphToBigWig {input.sense} {input.chrsizes} {output.sense}
        bedGraphToBigWig {input.antisense} {input.chrsizes} {output.antisense}
        """

rule deeptools_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
        bw = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-qfrags-{norm}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-QFRAGS.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-QFRAGS.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{factor}-{norm}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        """

rule deeptools_matrix_stranded:
    input:
        annotation = "../genome/annotations/stranded/{annotation}-STRANDED.bed",
        sense = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-SENSE.bw",
        antisense = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-ANTISENSE.bw"
    output:
        dtfile_sense = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-SENSE.mat.gz"),
        matrix_sense = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-SENSE.tsv"),
        dtfile_antisense = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-ANTISENSE.mat.gz"),
        matrix_antisense = temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-ANTISENSE.tsv")
    params:
        refpoint = lambda wildcards: config["annotations"][wildcards.annotation]["refpoint"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/deeptools/computeMatrix-{annotation}-{sample}-{factor}-{norm}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.sense} --referencePoint {params.refpoint} -out {output.dtfile_sense} --outFileNameMatrix {output.matrix_sense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}
        (computeMatrix reference-point -R {input.annotation} -S {input.antisense} --referencePoint {params.refpoint} -out {output.dtfile_antisense} --outFileNameMatrix {output.matrix_antisense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &>> {log}
        """

rule gzip_deeptools_table:
    input:
        tsv = "datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-{strand}.tsv"
    output:
        "datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-{strand}.tsv.gz"
    shell: """
        pigz -f {input.tsv}
        """

rule melt_matrix:
    input:
        matrix = "datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-{strand}.tsv.gz"
    output:
        temp("datavis/{annotation}/{norm}/{annotation}-{sample}-{factor}-chipnexus-{norm}-{strand}-melted.tsv.gz")
    params:
        name = lambda wildcards : wildcards.sample,
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"]
    script:
        "scripts/melt_matrix.R"

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}-{sample}-{{factor}}-chipnexus-{{norm}}-{{strand}}-melted.tsv.gz", sample=SAMPLES)
    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-{strand}.tsv.gz"
    shell: """
        cat {input} > {output}
        """

rule r_heatmaps:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-QFRAGS.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-heatmap-bysample.png",
        heatmap_group= "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-heatmap-bygroup.png",
    params:
        # binsize = lambda wildcards : config["annotations"][wildcards.annotation]["binsize"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        heatmap_height = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_height"],
        figwidth = lambda wildcards : config["annotations"][wildcards.annotation]["figwidth"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        factor = config["factor"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotHeatmaps.R"

rule r_metagenes:
    input:
        plus = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-SENSE.tsv.gz",
        minus = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-ANTISENSE.tsv.gz"
    output:
        meta_sample = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bysample.png",
        meta_group = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bygroup.png"
    params:
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        figwidth = lambda wildcards : config["annotations"][wildcards.annotation]["figwidth"],
        factor = config["factor"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotNexusMetagene.R"

rule get_putative_intragenic:
    input:
        peaks = "peakcalling/macs/{group}_peaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        tsv = "peakcalling/macs/{group}_peaks-intragenic.tsv",
        narrowpeak = "peakcalling/macs/{group}_peaks-intragenic.narrowPeak"
    shell: """
        sort -k1,1 -k2,2n -k3,3n -k9,9nr {input.peaks} | sort -k1,1 -k2,2n -k3,3n -u | bedtools intersect -a stdin -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -f 1 -wo | tee {output.tsv} | cut -f1-10 > {output.narrowpeak}
        """

rule intragenic_overlap:
    input:
        expand("peakcalling/macs/{group}_peaks-intragenic.narrowPeak", group=GROUPS)
    output:
        "peakcalling/macs/intragenic-peaks-overlap.tsv"
    shell: """
        bedtools multiinter -header -names {GROUPS} -cluster -i {input} > {output}
        """


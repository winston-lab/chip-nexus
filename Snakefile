#!/usr/bin/env python
import os
from math import log2, log10

configfile: "config.yaml"

SAMPLES = config["samples"]
PASSING = {k:v for (k,v) in SAMPLES.items() if v["pass-qc"] == "pass"}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = config["comparisons"]["libsizenorm"]["controls"]
conditiongroups = config["comparisons"]["libsizenorm"]["conditions"]
controlgroups_si = config["comparisons"]["spikenorm"]["controls"]
conditiongroups_si = config["comparisons"]["spikenorm"]["conditions"]

CATEGORIES = ["genic", "intragenic", "intergenic"]

localrules: all,
            make_stranded_annotations,
            make_combined_bedgraph,
            make_combined_si_bedgraph,
            get_si_pct, cat_si_pct,
            index_bam,
            bedgraph_to_bigwig,

rule all:
    input:
        #FastQC
        expand("qual_ctrl/fastqc/raw/{sample}", sample=SAMPLES),
        expand("qual_ctrl/fastqc/cleaned/{sample}/{sample}-clean_fastqc.zip", sample=SAMPLES),
        #alignment
        expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","qfrags"]),
        #initial QC
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-{{factor}}-chipnexus-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], factor=config["factor"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}-{{factor}}-chipnexus-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], factor=config["factor"]),
        #qnexus
        expand("peakcalling/qnexus/{sample}-{factor}-Q-treatment.bedgraph", sample=SAMPLES, factor=config["factor"]),
        #macs2
        expand("peakcalling/macs/{group}-{species}_peaks.narrowPeak", group = GROUPS, species=["Scer_","Spom_"]),
        #datavis
        # expand(expand("datavis/{{annotation}}/libsizenorm/{{factor}}-chipnexus-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}-heatmap-bygroup.svg", zip, condition=conditiongroups, control=controlgroups), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        # expand(expand("datavis/{{annotation}}/spikenorm/{{factor}}-chipnexus-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}-heatmap-bygroup.svg", zip, condition=conditiongroups_si, control=controlgroups_si), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        # expand("datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bygroup.svg", annotation = config["annotations"], norm = ["libsizenorm", "spikenorm"], factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","down"], fmt=["tsv", "bed"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","down"], fmt=["tsv", "bed"]),
        # expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","down"], category=CATEGORIES),
        # expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","down"], category=CATEGORIES),
        expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES)

def plotcorrsamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        if wildcards.norm=="libsizenorm": #condition==all,norm==lib
            return list(dd.keys())
        else: #condition==all,norm==spike
            return list({k:v for (k,v) in dd.items() if v["spikein"]=="y"}.keys())
    elif wildcards.norm=="libsizenorm": #condition!=all;norm==lib
        return list({k:v for (k,v) in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition}.keys())
    else: #condition!=all;norm==spike
        return list({k:v for (k,v) in dd.items() if (v["group"]==wildcards.control or v["group"]==wildcards.condition) and v["spikein"]=="y"}.keys())

rule fastqc_raw:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        "qual_ctrl/fastqc/raw/{sample}"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        (mkdir -p {output}) &> {log}
        (fastqc -o {output} --noextract -t {threads} {input}) &>> {log}
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
        (pigz -f fastq/cleaned/{wildcards.sample}-clean.fastq) &>> {log}
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
        (mkdir -p qual_ctrl/fastqc/cleaned/{wildcards.sample}) &> {log}
        (fastqc -o qual_ctrl/fastqc/cleaned/{wildcards.sample} --noextract -t {threads} {input}) &>> {log}
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

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    output:
        SIplmin = "coverage/sicounts/{sample}-{factor}-chipnexus-sicounts-combined.bedgraph",
        SIpl = "coverage/sicounts/{sample}-{factor}-chipnexus-sicounts-plus.bedgraph",
        SImin = "coverage/sicounts/{sample}-{factor}-chipnexus-sicounts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-{factor}-chipnexus-counts-combined.bedgraph",
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
        coverage = "coverage/counts/{sample}-{factor}-chipnexus-counts-qfrags.bedgraph"
    params:
        exp_prefix = config["combinedgenome"]["experimental_prefix"],
    log: "logs/qnexus/qnexus-{sample}-{factor}.log"
    shell: """
        (../programs/Q/bin/Q --nexus-mode -t {input} -o peakcalling/qnexus/{wildcards.sample}-{wildcards.factor} -v -wbt) &> {log}
        # (sed -i 's/peakcalling\///g' peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-qfrag-binding-characteristics.R) &>> {log}
        (Rscript peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-qfrag-binding-characteristics.R) &>> {log}
        (grep {params.exp_prefix} peakcalling/qnexus/{wildcards.sample}-{wildcards.factor}-Q-treatment.bedgraph | sed 's/{params.exp_prefix}//g' | sort -k1,1 -k2,2n > {output.coverage}) &>> {log}
        """

rule normalize:
    input:
        plus = "coverage/counts/{sample}-{factor}-chipnexus-counts-plus.bedgraph",
        minus = "coverage/counts/{sample}-{factor}-chipnexus-counts-minus.bedgraph",
        plmin = "coverage/counts/{sample}-{factor}-chipnexus-counts-combined.bedgraph",
        SIplmin = "coverage/sicounts/{sample}-{factor}-chipnexus-sicounts-combined.bedgraph",
        qfrags = "coverage/counts/{sample}-{factor}-chipnexus-counts-qfrags.bedgraph"
    output:
        spikePlus = "coverage/spikenorm/{sample}-{factor}-chipnexus-spikenorm-plus.bedgraph",
        spikeMinus = "coverage/spikenorm/{sample}-{factor}-chipnexus-spikenorm-minus.bedgraph",
        libnormPlus = "coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-plus.bedgraph",
        libnormMinus = "coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-minus.bedgraph",
        qfrag_spikenorm = "coverage/spikenorm/{sample}-{factor}-chipnexus-spikenorm-qfrags.bedgraph",
        qfrag_libnorm= "coverage/libsizenorm/{sample}-{factor}-chipnexus-libsizenorm-qfrags.bedgraph",
    log: "logs/normalize/normalize-{sample}.log"
    shell: """
        (scripts/libsizenorm.awk {input.SIplmin} {input.plus} > {output.spikePlus}) &> {log}
        (scripts/libsizenorm.awk {input.SIplmin} {input.minus} > {output.spikeMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.plus} > {output.libnormPlus}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.minus} > {output.libnormMinus}) &>> {log}
        (scripts/libsizenorm.awk {input.SIplmin} {input.qfrags} > {output.qfrag_spikenorm}) &>> {log}
        (scripts/libsizenorm.awk {input.plmin} {input.qfrags} > {output.qfrag_libnorm}) &>> {log}
        """

rule get_si_pct:
    input:
        plmin = "coverage/counts/{sample}-{factor}-chipnexus-counts-combined.bedgraph",
        SIplmin = "coverage/sicounts/{sample}-{factor}-chipnexus-sicounts-combined.bedgraph"
    output:
        temp("qual_ctrl/all/{sample}-{factor}-spikeincounts.tsv")
    params:
        group = lambda wildcards: SAMPLES[wildcards.sample]["group"]
    log: "logs/get_si_pct/get_si_pct-{sample}-{factor}.log"
    shell: """
        (echo -e "{wildcards.sample}\t{params.group}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {input.SIplmin} {input.plmin}) > {output}) &> {log}
        """

rule cat_si_pct:
    input:
        expand("qual_ctrl/all/{sample}-{factor}-spikeincounts.tsv", sample=SAMPLES, factor=config["factor"])
    output:
        "qual_ctrl/all/spikein-counts.tsv"
    log: "logs/cat_si_pct.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_si_pct:
    input:
        "qual_ctrl/all/spikein-counts.tsv"
    output:
        plot = "qual_ctrl/{status}/{status}-spikein-plots.svg",
        stats = "qual_ctrl/{status}/{status}-spikein-stats.tsv"
    params:
        samplelist = lambda wildcards : list({k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}.keys()) if wildcards.status=="all" else list({k:v for (k,v) in PASSING.items() if v["spikein"]=="y"}.keys()),
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
    script: "scripts/plotsipct.R"

rule index_bam:
    input:
        bam = "alignment/{sample}-noPCRdup.bam",
    output:
        "alignment/{sample}-noPCRdup.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input.bam}) &> {log}
        """

rule build_macs2_input:
    input:
        bam = "alignment/{sample}-noPCRdup.bam",
        bai = "alignment/{sample}-noPCRdup.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    output:
        "alignment/{sample}-{species}only.bam"
    log: "logs/build_macs2_input/build_macs2_input-{sample}-{species}.log"
    shell: """
        (samtools view -b {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') > {output}) &> {log}
        """

rule macs2:
    input:
        bam = lambda wildcards: expand("alignment/{sample}-{species}only.bam", sample= {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.group and v["pass-qc"] == "pass")}, species=wildcards.species),
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if wildcards.species==config["combinedgenome"]["experimental_prefix"] else config["genome"]["si-chrsizes"],
    output:
        xls = "peakcalling/macs/{group}-{species}_peaks.xls",
        peaks = "peakcalling/macs/{group}-{species}_peaks.narrowPeak",
        summits = "peakcalling/macs/{group}-{species}_summits.bed",
        script = "peakcalling/macs/{group}-{species}_model.r",
        pdf = "peakcalling/macs/{group}-{species}_model.pdf",
        treat_bg = "peakcalling/macs/{group}-{species}_treat_pileup.bdg",
        cntrl_bg = "peakcalling/macs/{group}-{species}_control_lambda.bdg"
    params:
        bw = config["macs2"]["bw"],
        slocal = config["macs2"]["slocal"],
        llocal = config["macs2"]["llocal"],
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.species==config["combinedgenome"]["experimental_prefix"] else config["combinedgenome"]["spikein_prefix"],
        qscore = config["macs2"]["fdr"]
    conda:
        "envs/python2.yaml"
    log: "logs/macs2/macs2-{group}-{species}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --keep-dup all --bdg -n peakcalling/macs/{wildcards.group}-{wildcards.species} --SPMR --bw {params.bw} --slocal {params.slocal} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript peakcalling/macs/{wildcards.group}-{wildcards.species}_model.r) &>> {log}
        (sed -i -e 's/{params.prefix}//g; s/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_peaks.narrowPeak) &>> {log}
        (sed -i -e 's/{params.prefix}//g; s/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_summits.bed) &>> {log}
        (sed -i -e 's/{params.prefix}//g' peakcalling/macs/{wildcards.group}-{wildcards.species}_treat_pileup.bdg) &>> {log}
        (sed -i -e 's/{params.prefix}//g' peakcalling/macs/{wildcards.group}-{wildcards.species}_control_lambda.bdg) &>> {log}
        """

# rule get_putative_intragenic:
#     input:
#         peaks = "peakcalling/macs/{group}_peaks.narrowPeak",
#         orfs = config["genome"]["orf-annotation"],
#         genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
#     output:
#         tsv = "peakcalling/macs/{group}_peaks-intragenic.tsv",
#         narrowpeak = "peakcalling/macs/{group}_peaks-intragenic.narrowPeak"
#     shell: """
#         sort -k1,1 -k2,2n -k3,3n -k9,9nr {input.peaks} | sort -k1,1 -k2,2n -k3,3n -u | bedtools intersect -a stdin -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -f 1 -wo | tee {output.tsv} | cut -f1-10 > {output.narrowpeak}
#         """
rule classify_peaks_genic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "peaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/genic/{group}-exp-peaks-genic.tsv"
    shell: """
        sed 's/{params.prefix}//g' {input.peaks} | bedtools intersect -a stdin -b {input.annotation} -wo | cut --complement -f11-13,15-17 > {output}
        """

rule classify_peaks_intragenic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "peaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/intragenic/{group}-exp-peaks-intragenic.tsv"
    shell: """
        sed 's/{params.prefix}//g' {input.peaks} | bedtools intersect -a stdin -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -wo | awk 'BEGIN{{FS=OFS="\t"}} $16=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $16=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_intergenic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "peaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/intergenic/{group}-exp-peaks-intergenic.tsv"
    shell: """
        sed 's/{params.prefix}//g' {input.peaks} | bedtools intersect -a stdin -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v | bedtools intersect -a stdin -b {input.annotation} -wa > {output}
        """

rule combine_peaks:
    input:
        expand("peakcalling/macs/{group}-{{species}}_peaks.narrowPeak", group=GROUPS)
    output:
        "peakcalling/macs/allpeaks-{species}.bed"
    shell: """
        bedtools multiinter -i {input} | bedtools merge -i stdin | sort -k1,1 -k2,2n > {output}
        """

rule map_counts_to_peaks:
    input:
        bed = "peakcalling/macs/allpeaks-{species}.bed",
        bg = lambda wildcards: "coverage/counts/" + wildcards.sample +"-"+ wildcards.factor+"-chipnexus-counts-combined.bedgraph" if wildcards.species==config["combinedgenome"]["experimental_prefix"] else "coverage/sicounts/" + wildcards.sample + "-" + wildcards.factor + "-chipnexus-sicounts-combined.bedgraph"
    output:
        temp("coverage/counts/{sample}-{factor}-{species}-peak-counts.tsv")
    shell: """
        bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}
        """

rule join_peak_counts:
    input:
        expand("coverage/counts/{sample}-{factor}-{{species}}-peak-counts.tsv", sample=SAMPLES, factor=config["factor"])
    output:
        "coverage/counts/union-bedgraph-allpeakcounts-{species}.tsv.gz"
    params:
        names = list(SAMPLES.keys())
    shell: """
        bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output}
        """

rule call_de_peaks:
    input:
        expcounts = "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["experimental_prefix"] + ".tsv.gz",
        sicounts = lambda wildcards: "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["spikein_prefix"] + ".tsv.gz" if wildcards.norm=="spikenorm" else "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["experimental_prefix"] + ".tsv.gz"
    params:
        samples = lambda wildcards : list({k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}.keys()),
        groups = lambda wildcards : [PASSING[x]["group"] for x in {k:v for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)}],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"]),
    output:
        results = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv",
        normcounts = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-peakcounts-sfnorm-{norm}.tsv",
        rldcounts = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-peakcounts-rlog-{norm}.tsv",
        qcplots = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-qcplots-{norm}.svg"
    script:
        "scripts/call_de_peaks.R"

rule separate_de_peaks:
    input:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv"
    output:
        up = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-up.tsv",
        down = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-down.tsv"
    params:
        fdr = -log10(config["deseq"]["fdr"])
    log: "logs/separate_de_peaks/separate_de_peaks-{condition}-v-{control}-{norm}.log"
    shell: """
        (awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}} NR==1{{print > "{output.up}"; print > "{output.down}" }} NR>1 && $10>afdr && $7>0 {{print > "{output.up}"}} NR>1 && $10>afdr && $7<0 {{print > "{output.down}"}}' {input}) &> {log}
        """

rule de_peaks_to_bed:
    input:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.tsv",
    output:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.bed",
    log: "logs/de_peaks_to_bed/de_peaks_to_bed-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        tail -n +2 {input} | awk 'BEGIN{{FS=OFS="\t"}}{{print $2, $4, $5, $1, $7":"$11, $3}}' > {output}
        """

rule build_genic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        windowsize = config["genic-windowsize"]
    log : "logs/build_genic_annotation.log"
    shell: """
        (python scripts/make_genic_annotation.py -t {input.transcripts} -o {input.orfs} -d {params.windowsize} -g {input.chrsizes} -p {output}) &> {log}
        """

rule get_putative_genic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_binding/{condition}-v-{control}/genic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-genic.tsv"
    log : "logs/get_putative_genic/get_putative_genic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | awk 'BEGIN{{FS="\t|:";OFS="\t"}}{{print $1, $13, $2, $3, $4, $9, $10, $11, $5, $6}}' | sort -k10,10nr | cat <(echo -e "chrom\ttranscript_strand\tpeak_start\tpeak_end\tpeak_name\ttranscript_start\ttranscript_end\ttranscript_name\tpeak_lfc\tpeak_significance") - > {output}) &> {log}
        """

rule build_intergenic_annotation:
    input:
        transcripts = config["genome"]["transcripts"],
        chrsizes = config["genome"]["chrsizes"]
    output:
        os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    params:
        genic_up = config["genic-windowsize"]
    log: "logs/build_intergenic_annotation.log"
    shell: """
        (bedtools slop -s -l {params.genic_up} -r 0 -i {input.transcripts} -g <(sort -k1,1 {input.chrsizes})| sort -k1,1 -k2,2n | bedtools complement -i stdin -g <(sort -k1,1 {input.chrsizes}) > {output}) &> {log}
        """

rule get_putative_intergenic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"
    output:
        "diff_binding/{condition}-v-{control}/intergenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-intergenic.tsv",
    log : "logs/get_putative_intergenic/get_putative_intergenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | awk 'BEGIN{{FS="\t|:";OFS="\t"}}{{print $1, ".", $2, $3, $4, $9, $10, ".", $5, $6}}'| sort -k10,10nr | cat <(echo -e "chrom\tpeak_strand\tpeak_start\tpeak_end\tpeak_name\tregion_start\tregion_end\tregion_name\tregion_name\tpeak_lfc\tpeak_significance") - > {output}) &> {log}
        """

rule get_putative_intragenic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.bed",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    output:
        "diff_binding/{condition}-v-{control}/intragenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-intragenic.tsv",
    log: "logs/get_putative_intragenic/get_putative_intragenic-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -wo | awk 'BEGIN{{FS="\t|:";OFS="\t"}} $13=="+"{{print $1, $13, $2, $3, $4, $9, $10, $11, $5, $6, ((($2+1)+$3)/2)-$9}} $13=="-"{{print $1, $13, $2, $3, $4, $9, $10, $11, $5, $6, $10-((($2+1)+$3)/2)}}' | sort -k10,10nr | cat <(echo -e  "chrom\torf_strand\tpeak_start\tpeak_end\tpeak_name\torf_start\torf_end\torf_name\tpeak_lfc\tpeak_significance\tpeak_dist_from_ATG") - > {output}) &> {log}
        """

rule get_category_bed:
    input:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $1, $3, $4, $5, $10, $2}}' {input} | sort -k1,1 -k2,2n  > {output}) &> {log}
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

rule make_stranded_annotations:
    input:
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "../genome/annotations/stranded/{annotation}-STRANDED.bed"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
        """

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-SENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        exp = temp("coverage/{norm}/{sample}-{factor}-window-coverage-{norm}.bedgraph"),
    params:
        windowsize = config["corr-windowsize"]
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {params.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        """

rule join_window_counts:
    input:
        exp = expand("coverage/{{norm}}/{sample}-{factor}-window-coverage-{{norm}}.bedgraph", sample=SAMPLES, factor=config["factor"]),
    output:
        exp = "coverage/{norm}/union-bedgraph-allwindowcoverage-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz > {output.exp}) &> {log}
        """

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-allwindowcoverage-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}-{factor}-chipnexus-{norm}-correlations.svg"
    params:
        pcount = 0.1,
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule bedgraph_to_bigwig:
    input:
        bg = "coverage/{norm}/{sample}-{factor}-chipnexus-{norm}-{strand}.bedgraph",
        chrsizes = lambda wildcards: os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv" if (wildcards.strand=="SENSE" or wildcards.strand=="ANTISENSE") else config["genome"]["chrsizes"]
    output:
        "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-{strand}.bw"
    shell: """
        bedGraphToBigWig {input.bg} {input.chrsizes} {output}
        """

#TODO: combine the two matrix rules into one
rule deeptools_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"],
        bw = "coverage/{norm}/bw/{sample}-{factor}-chipnexus-{norm}-qfrags.bw"
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
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")

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
    log: "logs/deeptools/computeMatrix-stranded-{annotation}-{sample}-{factor}-{norm}.log"
    run:
        if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.sense} --referencePoint {params.refpoint} -out {output.dtfile_sense} --outFileNameMatrix {output.matrix_sense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.antisense} --referencePoint {params.refpoint} -out {output.dtfile_antisense} --outFileNameMatrix {output.matrix_antisense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --nanAfterEnd --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &>> {log}")
        else:
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.sense} --referencePoint {params.refpoint} -out {output.dtfile_sense} --outFileNameMatrix {output.matrix_sense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}")
            shell("(computeMatrix reference-point -R {input.annotation} -S {input.antisense} --referencePoint {params.refpoint} -out {output.dtfile_antisense} --outFileNameMatrix {output.matrix_antisense} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &>> {log}")

rule gzip_deeptools_matrix:
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
        heatmap_sample = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}-heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}-heatmap-bygroup.svg",
    params:
        samplelist = plotcorrsamples,
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        cluster = lambda wildcards : config["annotations"][wildcards.annotation]["cluster"],
        nclust = lambda wildcards: config["annotations"][wildcards.annotation]["nclusters"],
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
        meta_sample = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bysample.svg",
        meta_group = "datavis/{annotation}/{norm}/{factor}-chipnexus-{annotation}-{norm}-metagene-bygroup.svg"
    params:
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        figwidth = lambda wildcards : config["annotations"][wildcards.annotation]["figwidth"],
        factor = config["factor"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    script:
        "scripts/plotNexusMetagene.R"


rule intragenic_overlap:
    input:
        expand("peakcalling/macs/{group}_peaks-intragenic.narrowPeak", group=GROUPS)
    output:
        "peakcalling/macs/intragenic-peaks-overlap.tsv"
    shell: """
        bedtools multiinter -header -names {GROUPS} -cluster -i {input} > {output}
        """


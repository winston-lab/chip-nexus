#!/usr/bin/env python
import os
from math import log2, log10
import itertools

configfile: "config.yaml"
subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

FACTOR = config["factor"]

SAMPLES = config["samples"]
sisamples = {k:v for k,v in SAMPLES.items() if v["spikein"]=="y"}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"] == "pass"}
sipassing = {k:v for k,v in PASSING.items() if v["spikein"]=="y"}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items()]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items()]
controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items()]
conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items()]

CATEGORIES = ["genic", "intragenic", "intergenic"]

FIGURES = config["figures"]

localrules:
    all,
    fastqc_aggregate,
    make_stranded_annotations,
    classify_peaks_genic, classify_peaks_intragenic, classify_peaks_intergenic,
    separate_de_peaks,
    get_de_genic, get_de_intragenic, get_de_intergenic,
    de_peaks_to_bed, separate_sig_de, get_de_category_bed,
    cat_matrices,
    make_ratio_annotation, cat_ratio_counts

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

rule all:
    input:
        #FastQC
        'qual_ctrl/fastqc/' + FACTOR + '-chipnexus-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam", sample=SAMPLES, factor=FACTOR),
        #macs2
        expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.narrowPeak", group=GROUPS, species=["experimental", "spikein"], factor=FACTOR),
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        #read processing stats
        "qual_ctrl/read_processing/" + FACTOR + "-chipnexus_read-processing-loss.svg",
        expand("qual_ctrl/spikein/{factor}-chipnexus_spikein-plots-{status}.svg", factor=FACTOR, status=["all","passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        ##categorise peaks
        #expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES),
        #expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        # datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, factor=config["factor"], status=["all", "passing"]),
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, factor=config["factor"], status=["all", "passing"]),
        ##differential binding of peaks
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","unchanged","down"], fmt=["tsv", "bed"]),
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","unchanged", "down"], fmt=["tsv", "bed"]),
        ##categorize DB peaks
        #expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","unchanged","down"], category=CATEGORIES),
        #expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","unchanged","down"], category=CATEGORIES),
        ##DB summary
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-libsizenorm-diffbind-summary.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-spikenorm-diffbind-summary.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        #expand(expand("ratios/{{ratio}}/{condition}-v-{control}/{{factor}}-chipnexus-{{ratio}}_{{status}}_{condition}-v-{control}_violin.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"], ratio=config["ratios"], status=["all", "passing"])

def get_condition_control_samples(wc):
    if wc.condition=="all":
        if wc.norm=="libsizenorm": #condition==all,norm==lib
            return list(SAMPLES.keys())
        else: #condition==all,norm==spike
            return list(sisamples.keys())
    elif wc.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in PASSING.items() if v["group"] in [wc.control, wc.condition]]
    else: #condition!=all;norm==spike
        return [k for k,v in sipassing.items() if v["group"] in [wc.control, wc.condition]]

def cluster_samples(status, norm, cluster_groups, cluster_strands):
    ll = []
    dd = SAMPLES if status=="all" else PASSING
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in dd.items() if v["group"]==group] if norm=="libsizenorm" else [k for k,v in dd.items() if v["group"]==group and v["spikein"]=="y"]
        if strand=="protection":
            ll.append([sample + "-" + "protection" for sample in sublist])
        if strand in ["sense", "both"]:
            ll.append([sample + "-" + "sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([sample + "-" + "antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))


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

rule map_to_windows:
    input:
        bg = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-SENSE.bedgraph",
        chrsizes = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
    output:
        exp = temp("coverage/{norm}/{sample}_{factor}-window-{windowsize}-coverage-{norm}.bedgraph"),
    shell: """
        bedtools makewindows -g {input.chrsizes} -w {wildcards.windowsize} | LC_COLLATE=C sort -k1,1 -k2,2n | bedtools map -a stdin -b {input.bg} -c 4 -o sum > {output.exp}
        """

rule join_window_counts:
    input:
        exp = expand("coverage/{{norm}}/{sample}_{factor}-window-{{windowsize}}-coverage-{{norm}}.bedgraph", sample=SAMPLES, factor=config["factor"]),
    output:
        exp = "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz",
    params:
        names = list(SAMPLES.keys())
    log: "logs/join_window_counts/join_window_counts-{norm}.log"
    shell: """
        (bedtools unionbedg -i {input.exp} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output.exp}) &> {log}
        """

rule plotcorrelations:
    input:
        "coverage/{norm}/union-bedgraph-window-{windowsize}-{norm}.tsv.gz"
    output:
        "qual_ctrl/{status}/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{status}-window-{windowsize}-{norm}-correlations.svg"
    params:
        pcount = lambda wc: 0.01*int(wc.windowsize),
        samplelist = get_condition_control_samples
    script:
        "scripts/plotcorr.R"

rule make_ratio_annotation:
    input:
        lambda wc: config["ratios"][wc.ratio]["path"]
    params:
        totalsize = lambda wc: config["ratios"][wc.ratio]["numerator"]["upstream"] + config["ratios"][wc.ratio]["numerator"]["dnstream"] + config["ratios"][wc.ratio]["denominator"]["upstream"] + config["ratios"][wc.ratio]["denominator"]["dnstream"],
    output:
        "ratios/{ratio}/{ratio}.bed"
    log: "logs/make_ratio_annotation/make_ratio_annotation-{ratio}.log"
    shell:  """
        (awk 'BEGIN{{FS=OFS="\t"}} ($3-$2)>={params.totalsize}' {input} > {output}) &> {log}
        """

rule ratio_counts:
    input:
        annotation = "ratios/{ratio}/{ratio}.bed",
        bw = "coverage/libsizenorm/{sample}_" + config["factor"] + "-chipnexus-libsizenorm-midpoints.bw"
    output:
        dtfile = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}.mat.gz"),
        matrix = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}.tsv"),
        melted = temp("ratios/{ratio}/{ratio}_{fractype}_{sample}-melted.tsv.gz"),
    params:
        group = lambda wc : SAMPLES[wc.sample]["group"],
        upstream = lambda wc: config["ratios"][wc.ratio][wc.fractype]["upstream"],
        dnstream = lambda wc: config["ratios"][wc.ratio][wc.fractype]["dnstream"],
        refpoint = lambda wc: config["ratios"][wc.ratio][wc.fractype]["refpoint"]
    threads: config["threads"]
    log: "logs/ratio_counts/ratio_counts-{ratio}-{fractype}-{sample}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize $(echo {params.upstream} + {params.dnstream} | bc) --averageTypeBins sum -p {threads}) &> {log}
        (Rscript scripts/melt_matrix.R -i {output.matrix} -r TSS --group {params.group} -s {wildcards.sample} -a none -b $(echo {params.upstream} + {params.dnstream} | bc) -u {params.upstream} -o {output.melted}) &>> {log}
        """

rule cat_ratio_counts:
    input:
        expand("ratios/{{ratio}}/{{ratio}}_{{fractype}}_{sample}-melted.tsv.gz", sample=SAMPLES)
    output:
        "ratios/{ratio}/allsamples_{ratio}_{fractype}.tsv.gz"
    log: "logs/cat_ratio_counts/cat_ratio_counts-{ratio}-{fractype}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

def ratiosamples(wc):
    dd = SAMPLES if wc.status=="all" else PASSING
    if wc.condition=="all":
        return list(dd.keys())
    else:
        return [k for k,v in dd.items() if v["group"]==wc.control or v["group"]==wc.condition]

rule plot_ratios:
    input:
        numerator = "ratios/{ratio}/allsamples_{ratio}_numerator.tsv.gz",
        denominator = "ratios/{ratio}/allsamples_{ratio}_denominator.tsv.gz",
    output:
        violin = "ratios/{ratio}/{condition}-v-{control}/{factor}-chipnexus-{ratio}_{status}_{condition}-v-{control}_violin.svg",
        ecdf = "ratios/{ratio}/{condition}-v-{control}/{factor}-chipnexus-{ratio}_{status}_{condition}-v-{control}_ecdf.svg"
    params:
        num_size = lambda wc: config["ratios"][wc.ratio]["numerator"]["upstream"] + config["ratios"][wc.ratio]["numerator"]["dnstream"],
        den_size = lambda wc: config["ratios"][wc.ratio]["denominator"]["upstream"] + config["ratios"][wc.ratio]["denominator"]["dnstream"],
        pcount = 1e-3,
        samplelist = ratiosamples,
        ratio_label = lambda wc: config["ratios"][wc.ratio]["ratio_name"],
        num_label = lambda wc: config["ratios"][wc.ratio]["numerator"]["region_label"],
        den_label = lambda wc: config["ratios"][wc.ratio]["denominator"]["region_label"],
        annotation_label = lambda wc: config["ratios"][wc.ratio]["label"]
    script:
        "scripts/ratio.R"

include: "rules/chip-nexus_clean_reads.smk"
include: "rules/chip-nexus_alignment.smk"
include: "rules/chip-nexus_fastqc.smk"
include: "rules/chip-nexus_library-processing-summary.smk"
include: "rules/chip-nexus_peakcalling.smk"
include: "rules/chip-nexus_genome_coverage.smk"
include: "rules/chip-nexus_datavis.smk"


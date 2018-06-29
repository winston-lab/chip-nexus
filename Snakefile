#!/usr/bin/env python

import os
from math import log2, log10
import itertools

configfile: "config.yaml"
subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

FACTOR = config["factor"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = [v for k,v in config["comparisons"]["libsizenorm"].items()]
conditiongroups = [k for k,v in config["comparisons"]["libsizenorm"].items()]
controlgroups_si = [v for k,v in config["comparisons"]["spikenorm"].items()]
conditiongroups_si = [k for k,v in config["comparisons"]["spikenorm"].items()]

CATEGORIES = ["genic", "intragenic", "intergenic"]

FIGURES = config["figures"]

status_norm_sample_dict = {
    "all":
        {   "libsizenorm" : SAMPLES,
            "spikenorm" : SISAMPLES
        },
    "passing":
        {   "libsizenorm" : PASSING,
            "spikenorm" : SIPASSING
        }
    }

def get_samples(status, norm, groups):
    if "all" in groups:
        return(list(status_norm_sample_dict[status][norm].keys()))
    else:
        return([k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in groups])

def cluster_samples(status, norm, cluster_groups, cluster_strands):
    ll = []
    for group, strand in zip(cluster_groups, cluster_strands):
        sublist = [k for k,v in status_norm_sample_dict[status][norm].items() if v["group"] in cluster_groups]
        if strand=="protection":
            ll.append([f"{sample}-protection" for sample in sublist])
        if strand in ["sense", "both"]:
            ll.append([f"{sample}-sense" for sample in sublist])
        if strand in ["antisense", "both"]:
            ll.append([f"{sample}-antisense" for sample in sublist])
    return(list(itertools.chain(*ll)))

include: "rules/chip-nexus_clean_reads.smk"
include: "rules/chip-nexus_alignment.smk"
include: "rules/chip-nexus_fastqc.smk"
include: "rules/chip-nexus_library-processing-summary.smk"
include: "rules/chip-nexus_sample_similarity.smk"
include: "rules/chip-nexus_peakcalling.smk"
include: "rules/chip-nexus_classify_peaks.smk"
include: "rules/chip-nexus_differential_binding.smk"
include: "rules/chip-nexus_genome_coverage.smk"
include: "rules/chip-nexus_datavis.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules:
    all,
    make_stranded_annotations,
    make_ratio_annotation, cat_ratio_counts

rule all:
    input:
        #FastQC
        f'qual_ctrl/fastqc/{FACTOR}-chipnexus-per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam", sample=SAMPLES, factor=FACTOR),
        #macs2
        expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.narrowPeak", group=GROUPS, species=["experimental", "spikein"], factor=FACTOR),
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=FACTOR, norm=["counts","libsizenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SISAMPLES, factor=FACTOR, norm=["sicounts", "spikenorm"], strand=["plus","minus","protection","midpoints"]),
        #read processing stats
        f"qual_ctrl/read_processing/{FACTOR}-chipnexus_read-processing-loss.svg",
        #scatterplots
        expand("qual_ctrl/spikein/{factor}-chipnexus_spikein-plots-{status}.svg", factor=FACTOR, status=["all","passing"]) if SISAMPLES else [],
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all","passing"], factor=FACTOR, windowsize=config["corr-windowsizes"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all","passing"], factor=FACTOR, windowsize=config["corr-windowsizes"]) if SISAMPLES else [],
        ##categorise peaks
        expand("peakcalling/macs/{group}/{group}_experimental-{factor}-chipnexus_peaks-{category}.narrowpeak", group=GROUPS, factor=FACTOR, category=CATEGORIES),
        #expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        # datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, factor=FACTOR, status=["all", "passing"]) if config["plot_figures"] else [],
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, factor=FACTOR, status=["all", "passing"]) if config["plot_figures"] and SISAMPLES else [],
        #differential binding of peaks
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-qc-plots.svg", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-qc-plots.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR),
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-results-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR, direction=["all","up","unchanged","down"]),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-results-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR, direction=["all","up","unchanged","down"]),
        #categorize DB peaks
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", condition=conditiongroups, control=controlgroups), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{{category}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-results-{{category}}-{{direction}}.narrowpeak", condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR, direction=["all","up","unchanged","down"], category=CATEGORIES),
        ##DB summary
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-libsizenorm-diffbind-summary.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        #expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-spikenorm-diffbind-summary.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        #expand(expand("ratios/{{ratio}}/{condition}-v-{control}/{{factor}}-chipnexus-{{ratio}}_{{status}}_{condition}-v-{control}_violin.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"], ratio=config["ratios"], status=["all", "passing"])

rule make_stranded_genome:
    input:
        exp = config["genome"]["chrsizes"],
        si = config["genome"]["sichrsizes"]
    output:
        exp = os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv",
        si= os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv",
    log: "logs/make_stranded_genome.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.exp} > {output.exp}) &> {log}
        (awk 'BEGIN{{FS=OFS="\t"}}{{print $1"-plus", $2}}{{print $1"-minus", $2}}' {input.si} > {output.si}) &> {log}
        """

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


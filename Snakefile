#!/usr/bin/env python
import os
import re
import itertools
from math import log2, log10

configfile: "config.yaml"

subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

configfile: build_annotations("config.yaml")

FACTOR = config["factor"]

SAMPLES = config["samples"]
SISAMPLES = {k:v for k,v in SAMPLES.items() if v["spikein"]}
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"]}
SIPASSING = {k:v for k,v in PASSING.items() if v["spikein"]}
GROUPS = set(v["group"] for (k,v) in SAMPLES.items())

controlgroups = list(itertools.chain(*[d.values() for d in config["comparisons"]["libsizenorm"]]))
conditiongroups = list(itertools.chain(*[d.keys() for d in config["comparisons"]["libsizenorm"]]))

comparisons_si = config["comparisons"]["spikenorm"]
if comparisons_si:
    controlgroups_si = list(itertools.chain(*[d.values() for d in config["comparisons"]["spikenorm"]]))
    conditiongroups_si = list(itertools.chain(*[d.keys() for d in config["comparisons"]["spikenorm"]]))

FIGURES = config["figures"]

wildcard_constraints:
    sample = "|".join(re.escape(x) for x in list(SAMPLES.keys())),
    group = "|".join(re.escape(x) for x in GROUPS),
    control = "|".join(set(re.escape(x) for x in controlgroups + (controlgroups_si if comparisons_si else []) + ["all"])),
    condition = "|".join(set(re.escape(x) for x in conditiongroups + (conditiongroups_si if comparisons_si else []) + ["all"])),
    species = "experimental|spikein",
    read_status = "raw|cleaned|aligned_noPCRdup|unaligned",
    figure = "|".join(re.escape(x) for x in list(FIGURES.keys())),
    annotation = "|".join(re.escape(x) for x in set(itertools.chain(*[FIGURES[figure]["annotations"].keys() for figure in FIGURES]))),
    status = "all|passing",
    counttype= "counts|sicounts",
    norm = "counts|sicounts|libsizenorm|spikenorm",
    strand = "SENSE|ANTISENSE|plus|minus|midpoints|protection",
    windowsize = "\d+",
    direction = "all|up|unchanged|down",
    factor=FACTOR

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
include: "rules/chip-nexus_differential_binding.smk"
include: "rules/chip-nexus_genome_coverage.smk"
include: "rules/chip-nexus_datavis.smk"

onsuccess:
    shell("(./mogrify.sh) > mogrify.log")

localrules: all

def statuscheck(dict1, dict2):
    return(["passing"] if dict1 == dict2 else ["all", "passing"])

def conditioncheck(conditionlist):
    return(conditionlist if len(conditionlist)==1 else conditionlist + ["all"])

rule all:
    input:
        #require config file so that it gets archived
        "config.yaml",
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
        expand("qual_ctrl/spikein/{factor}-chipnexus_spikein-plots-{status}.svg", factor=FACTOR, status=statuscheck(SISAMPLES, SIPASSING)) if SISAMPLES else [],
        #scatterplots
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), status=statuscheck(SAMPLES, PASSING), factor=FACTOR, windowsize=config["scatterplot_binsizes"]),
        expand(expand("qual_ctrl/scatter_plots/{condition}-v-{control}/{{status}}/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-scatterplots-{{status}}-window-{{windowsize}}.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), status=statuscheck(SISAMPLES, SIPASSING), factor=FACTOR, windowsize=config["scatterplot_binsizes"]) if SISAMPLES and comparisons_si else [],
        # datavis
        expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditioncheck(conditiongroups), control=conditioncheck(controlgroups)), figure=FIGURES, factor=FACTOR, status=statuscheck(SAMPLES, PASSING)) if config["plot_figures"] else [],
        expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus_{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditioncheck(conditiongroups_si), control=conditioncheck(controlgroups_si)), figure=FIGURES, factor=FACTOR, status=statuscheck(SISAMPLES, SIPASSING)) if config["plot_figures"] and SISAMPLES and comparisons_si else [],
        #differential binding of peaks
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-qc-plots.svg", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-qc-plots.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR) if comparisons_si else [],
        expand(expand("diff_binding/{condition}-v-{control}/libsizenorm/{condition}-v-{control}_{{factor}}-chipnexus-libsizenorm-diffbind-results-{{direction}}.narrowpeak", zip, condition=conditiongroups, control=controlgroups), factor=FACTOR, direction=["all","up","unchanged","down"]),
        expand(expand("diff_binding/{condition}-v-{control}/spikenorm/{condition}-v-{control}_{{factor}}-chipnexus-spikenorm-diffbind-results-{{direction}}.narrowpeak", zip, condition=conditiongroups_si, control=controlgroups_si), factor=FACTOR, direction=["all","up","unchanged","down"]) if comparisons_si else [],


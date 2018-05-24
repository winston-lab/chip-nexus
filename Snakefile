#!/usr/bin/env python
import os
from math import log2, log10
import itertools

configfile: "config.yaml"
subworkflow build_annotations:
    workdir: config["genome"]["annotation_workflow"]

FACTOR = config["factor"]

SAMPLES = config["samples"]
PASSING = {k:v for k,v in SAMPLES.items() if v["pass-qc"] == "pass"}
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
    get_si_pct, cat_si_pct,
    index_bam,
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
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        ##initial QC
        #"qual_ctrl/read_processing-loss.svg",
        #expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        #macs2
        expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.narrowPeak", group=GROUPS, species=["experimental", "spikein"], factor=FACTOR),
        ##categorise peaks
        #expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES),
        #expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        ## datavis
        #expand(expand("datavis/{{figure}}/libsizenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus-{{figure}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), figure=FIGURES, factor=config["factor"], status=["all", "passing"]),
        #expand(expand("datavis/{{figure}}/spikenorm/{condition}-v-{control}/{{status}}/{{factor}}-chipnexus-{{figure}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), figure=FIGURES, factor=config["factor"], status=["all", "passing"]),
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

def plotcorrsamples(wc):
    if wc.condition=="all":
        if wc.norm=="libsizenorm": #condition==all,norm==lib
            return list(SAMPLES.keys())
        else: #condition==all,norm==spike
            return [k for k,v in SAMPLES.items() if v["spikein"]=="y"]
    elif wc.norm=="libsizenorm": #condition!=all;norm==lib
        return [k for k,v in PASSING.items() if v["group"] in (wc.control, wc.condition)]
    else: #condition!=all;norm==spike
        return [k for k,v in PASSING.items() if v["group"] in (wc.control, wc.condition) and v["spikein"]=="y"]

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

rule make_stranded_annotations:
    input:
        lambda wc : FIGURES[wc.figure]["annotations"][wc.annotation]["path"]
    output:
        "{annopath}/stranded/{figure}_{annotation}-STRANDED.{ext}"
    log : "logs/make_stranded_annotations/make_stranded_annotations-{annotation}.log"
    shell: """
        (bash scripts/makeStrandedBed.sh {input} > {output}) &> {log}
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
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule compute_matrix:
    input:
        annotation = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["path"] if wc.strand=="protection" else os.path.dirname(FIGURES[wc.figure]["annotations"][wc.annotation]["path"]) + "/stranded/" + wc.figure + "_" + wc.annotation + "-STRANDED" + os.path.splitext(FIGURES[wc.figure]["annotations"][wc.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{figure}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{figure}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}.tsv"),
        melted = "datavis/{figure}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}-melted.tsv.gz",
    params:
        group = lambda wc : SAMPLES[wc.sample]["group"],
        refpoint = lambda wc: "TSS" if FIGURES[wc.figure]["parameters"]["type"]=="scaled" else FIGURES[wc.figure]["parameters"]["refpoint"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"] + FIGURES[wc.figure]["parameters"]["binsize"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        binsize = lambda wc: FIGURES[wc.figure]["parameters"]["binsize"],
        binstat = lambda wc: FIGURES[wc.figure]["parameters"]["binstat"],
        nan_afterend = lambda wc: [] if FIGURES[wc.figure]["parameters"]["type"]=="scaled" or not FIGURES[wc.figure]["parameters"]["nan_afterend"] else "--nanAfterEnd",
        anno_label = lambda wc: FIGURES[wc.figure]["annotations"][wc.annotation]["label"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{figure}_{annotation}_{sample}_{factor}-{norm}-{strand}.log"
    run:
        if FIGURES[wildcards.figure]["parameters"]["type"]=="absolute":
            shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} {params.nan_afterend} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {params.scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {params.refpoint} -g {params.group} -s {wildcards.sample} -a {params.anno_label} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        lambda wc: expand("datavis/{figure}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}-melted.tsv.gz", annotation=list(FIGURES[wc.figure]["annotations"].keys()), sample=SAMPLES, factor=wc.factor, figure=wc.figure, norm=wc.norm, strand=wc.strand)
    output:
        "datavis/{figure}/{norm}/{figure}-allsamples-allannotations-{factor}-chipnexus-{norm}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{figure}_{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_figures:
    input:
        matrices = expand("datavis/{{figure}}/{{norm}}/{{figure}}-allsamples-allannotations-{{factor}}-chipnexus-{{norm}}-{strand}.tsv.gz", strand=["protection", "SENSE", "ANTISENSE"]),
        annotations = lambda wc: [v["path"] for k,v in FIGURES[wc.figure]["annotations"].items()]
    output:
        heatmap_sample = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample.svg",
        heatmap_group = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup.svg",
        metagene_sample_protection = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-protection.svg",
        metagene_sample_stranded = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bysample-coverage.svg",
        metagene_sample_overlay_protection = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-protection.svg",
        metagene_sample_overlay_stranded = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-sampleoverlay-coverage.svg",
        metagene_group_protection = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-protection.svg",
        metagene_group_stranded = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-bygroup-coverage.svg",
        metagene_sampleclust_protection = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-protection.svg",
        metagene_sampleclust_stranded = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustersample-coverage.svg",
        metagene_groupclust_protection = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-protection.svg",
        metagene_groupclust_stranded = "datavis/{figure}/{norm}/{condition}-v-{control}/{status}/{factor}-chipnexus-{figure}-{norm}-{status}_{condition}-v-{control}_metagene-byclustergroup-coverage.svg",
    params:
        # abusing snakemake a bit here...using params as output paths in order to use lambda functions
        annotations_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + annotation + "_cluster-" + str(cluster) + ".bed" for annotation in FIGURES[wc.figure]["annotations"] for cluster in range(1, FIGURES[wc.figure]["annotations"][annotation]["n_clusters"]+1)],
        clusters_out = lambda wc: ["datavis/" + wc.figure + "/" + wc.norm + "/" + wc.condition + "-v-" + wc.control + "/" + wc.status + "/" + annotation + ".pdf" for annotation in FIGURES[wc.figure]["annotations"]],
        samplelist = plotcorrsamples,
        plottype = lambda wc: FIGURES[wc.figure]["parameters"]["type"],
        upstream = lambda wc: FIGURES[wc.figure]["parameters"]["upstream"],
        dnstream = lambda wc: FIGURES[wc.figure]["parameters"]["dnstream"],
        scaled_length = lambda wc: 0 if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["scaled_length"],
        pct_cutoff = lambda wc: FIGURES[wc.figure]["parameters"]["pct_cutoff"],
        log_transform = lambda wc: str(FIGURES[wc.figure]["parameters"]["log_transform"]).upper(),
        pcount = lambda wc: 0 if not FIGURES[wc.figure]["parameters"]["log_transform"] else FIGURES[wc.figure]["parameters"]["pseudocount"],
        trim_pct = lambda wc: FIGURES[wc.figure]["parameters"]["trim_pct"],
        refpointlabel = lambda wc: FIGURES[wc.figure]["parameters"]["refpointlabel"],
        endlabel = lambda wc:  "HAIL SATAN" if FIGURES[wc.figure]["parameters"]["type"]=="absolute" else FIGURES[wc.figure]["parameters"]["three_prime_label"],
        cmap = lambda wc: FIGURES[wc.figure]["parameters"]["heatmap_colormap"],
        sortmethod = lambda wc: FIGURES[wc.figure]["parameters"]["arrange"],
        cluster_scale = lambda wc: "FALSE" if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else str(FIGURES[wc.figure]["parameters"]["cluster_scale"]).upper(),
        cluster_samples = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else cluster_samples(wc.status, wc.norm, FIGURES[wc.figure]["parameters"]["cluster_conditions"], FIGURES[wc.figure]["parameters"]["cluster_strands"]),
        cluster_five = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_five"],
        cluster_three = lambda wc: [] if FIGURES[wc.figure]["parameters"]["arrange"] != "cluster" else FIGURES[wc.figure]["parameters"]["cluster_three"],
        k = lambda wc: [v["n_clusters"] for k,v in FIGURES[wc.figure]["annotations"].items()],
    script:
        "scripts/plot_nexus_figures.R"

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
include: "rules/chip-nexus_align.smk"
include: "rules/chip-nexus_fastqc.smk"
include: "rules/chip-nexus_peakcalling.smk"
include: "rules/chip-nexus_genome_coverage.smk"


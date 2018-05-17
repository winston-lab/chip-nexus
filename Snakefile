#!/usr/bin/env python
import os
from math import log2, log10
import itertools

configfile: "config.yaml"

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
        ##FastQC
        #'qual_ctrl/fastqc/per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}_{factor}-chipnexus-noPCRduplicates.bam", sample=SAMPLES, factor=FACTOR),
        ##coverage
        #expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        #expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        ##initial QC
        #"qual_ctrl/read_processing-loss.svg",
        #expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        #expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        ##macs2
        #expand("peakcalling/macs/{group}-{species}_peaks.narrowPeak", group = GROUPS, species=[config["combinedgenome"]["experimental_prefix"],config["combinedgenome"]["spikein_prefix"]]),
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

rule fastqc_raw:
    input:
        lambda wc: SAMPLES[wc.sample]["fastq"]
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt"
    threads: config["threads"]
    log: "logs/fastqc/raw/fastqc-raw-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/raw/{wildcards.sample}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/raw/{wildcards.sample} {input}) &>> {log}
        """

rule fastqc_cleaned:
    input:
        "fastq/{fqtype}/{sample}-{fqtype}.fastq.gz"
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/{fqtype}/{sample}-{fqtype}_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/{fqtype}/fastqc-{fqtype}-{sample}.log"
    wildcard_constraints:
        fqtype = "cleaned|unaligned"
    shell: """
        (mkdir -p qual_ctrl/fastqc/{wildcards.fqtype}) &> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/{wildcards.fqtype} {input}) &>> {log}
        """

rule fastqc_aligned:
    input:
        "alignment/{sample}-noPCRdup.bam"
    params:
        adapter = config["cutadapt"]["adapter"]
    output:
        "qual_ctrl/fastqc/aligned_noPCRdup/{sample}-aligned_noPCRdup_fastqc/fastqc_data.txt",
    threads : config["threads"]
    log: "logs/fastqc/aligned_noPCRdup/fastqc-aligned_noPCRdup-{sample}.log"
    shell: """
        (mkdir -p qual_ctrl/fastqc/aligned_noPCRdup) &> {log}
        (bedtools bamtofastq -fq qual_ctrl/fastqc/aligned_noPCRdup/{wildcards.sample}-aligned_noPCRdup.fastq -i {input}) &>> {log}
        (fastqc -a <(echo -e "adapter\t{params.adapter}") --nogroup --extract -t {threads} -o qual_ctrl/fastqc/aligned_noPCRdup qual_ctrl/fastqc/aligned_noPCRdup/{wildcards.sample}-aligned_noPCRdup.fastq) &>> {log}
        (rm qual_ctrl/fastqc/aligned_noPCRdup/{wildcards.sample}-aligned_noPCRdup.fastq) &>> {log}
        """

rule fastqc_aggregate:
    input:
        raw = expand("qual_ctrl/fastqc/raw/{sample}/{fname}/fastqc_data.txt", zip, sample=SAMPLES, fname=[os.path.split(v["fastq"])[1].split(".fastq")[0] + "_fastqc" for k,v in SAMPLES.items()]),
        cleaned = expand("qual_ctrl/fastqc/cleaned/{sample}-cleaned_fastqc/fastqc_data.txt", sample=SAMPLES),
        aligned_noPCRdup = expand("qual_ctrl/fastqc/aligned_noPCRdup/{sample}-aligned_noPCRdup_fastqc/fastqc_data.txt", sample=SAMPLES),
        unaligned = expand("qual_ctrl/fastqc/unaligned/{sample}-unaligned_fastqc/fastqc_data.txt", sample=SAMPLES),
    output:
        'qual_ctrl/fastqc/per_base_quality.tsv',
        'qual_ctrl/fastqc/per_tile_quality.tsv',
        'qual_ctrl/fastqc/per_sequence_quality.tsv',
        'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        'qual_ctrl/fastqc/per_sequence_gc.tsv',
        'qual_ctrl/fastqc/per_base_n.tsv',
        'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        'qual_ctrl/fastqc/adapter_content.tsv',
        'qual_ctrl/fastqc/kmer_content.tsv'
    run:
        shell("rm -f {output}")
        #for each statistic
        for outpath, stat, header in zip(output, ["Per base sequence quality", "Per tile sequence quality", "Per sequence quality scores", "Per base sequence content", "Per sequence GC content", "Per base N content", "Sequence Length Distribution", "Total Deduplicated Percentage", "Adapter Content", "Kmer Content"], ["base\tmean\tmedian\tlower_quartile\tupper_quartile\tten_pct\tninety_pct\tsample\tstatus", "tile\tbase\tmean\tsample\tstatus",
        "quality\tcount\tsample\tstatus", "base\tg\ta\tt\tc\tsample\tstatus", "gc_content\tcount\tsample\tstatus", "base\tn_count\tsample\tstatus", "length\tcount\tsample\tstatus", "duplication_level\tpct_of_deduplicated\tpct_of_total\tsample\tstatus", "position\tpct\tsample\tstatus",
        "sequence\tcount\tpval\tobs_over_exp_max\tmax_position\tsample\tstatus" ]):
            for input_type in ["raw", "cleaned", "aligned_noPCRdup", "unaligned"]:
                for sample_id, fqc in zip(SAMPLES.keys(), input[input_type]):
                    shell("""awk -v sample_id={sample_id} -v input_type={input_type} 'BEGIN{{FS=OFS="\t"}} /{stat}/{{flag=1;next}}/>>END_MODULE/{{flag=0}} flag {{print $0, sample_id, input_type}}' {fqc} | tail -n +2 >> {outpath}""")
            shell("""sed -i "1i {header}" {outpath}""")

rule plot_fastqc_summary:
    input:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.tsv',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.tsv',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.tsv',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.tsv',
        per_base_n = 'qual_ctrl/fastqc/per_base_n.tsv',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.tsv',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.tsv',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.tsv',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.tsv',
        # kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        # kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

rule read_processing_numbers:
    input:
        adapter = expand("logs/remove_adapter/remove_adapter-{sample}.log", sample=SAMPLES),
        align = expand("logs/align/align-{sample}.log", sample=SAMPLES),
        nodups = expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES)
    output:
        "qual_ctrl/read_processing_summary.tsv"
    log: "logs/read_processing_summary.log"
    run:
        shell("""(echo -e "sample\traw\tcleaned\tmapped\tunique_map\tnoPCRdup" > {output}) &> {log}""")
        for sample, adapter, align, nodups in zip(SAMPLES.keys(), input.adapter, input.align, input.nodups):
            shell("""(grep -e "Total reads processed:" -e "Reads written" {adapter} | cut -d: -f2 | sed 's/,//g' | awk 'BEGIN{{ORS="\t"; print "{sample}"}}{{print $1}}' >> {output}) &> {log}""")
            shell("""(grep -e "aligned 0 times" -e "aligned exactly 1 time" {align} | awk 'BEGIN{{ORS="\t"}} {{print $1}}' >> {output}) &> {log}""")
            shell("""(samtools view -c {nodups} | awk '{{print $1}}' >> {output}) &> {log}""")
        shell("""(awk 'BEGIN{{FS=OFS="\t"}} NR==1; NR>1{{$4=$3-$4; print $0}}' {output} > qual_ctrl/.readnumbers.temp; mv qual_ctrl/.readnumbers.temp {output}) &> {log}""")

rule plot_read_processing:
    input:
        "qual_ctrl/read_processing_summary.tsv"
    output:
        surv_abs_out = "qual_ctrl/read_processing-survival-absolute.svg",
        surv_rel_out = "qual_ctrl/read_processing-survival-relative.svg",
        loss_out  = "qual_ctrl/read_processing-loss.svg",
    script: "scripts/processing_summary.R"


rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    params:
        prefix = lambda wc: config["combinedgenome"]["experimental_prefix"] if wc.counttype=="counts" else config["combinedgenome"]["spikein_prefix"],
    output:
        plus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-minus.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/get_coverage/get_coverage-{sample}-{counttype}.log"
    shell: """
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """

rule get_protection:
    input:
        tsv = lambda wc: expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["experimental_prefix"]) if wc.counttype=="counts" else expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["spikein_prefix"]),
        bam = lambda wc: "alignment/" + wc.sample + "-" + config["combinedgenome"]["experimental_prefix"] +"only.bam" if wc.counttype=="counts" else "alignment/" + wc.sample + "-" + config["combinedgenome"]["spikein_prefix"] +"only.bam"
    output:
        coverage = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-protection.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    shell: """
        median_fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4 | sort -k1,1n | awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 2.0;}} }}' | xargs printf "%.*f\n" 0)
        genomeCoverageBed -bga -fs $median_fragsize -scale $(echo 1/$median_fragsize | bc -l) -ibam {input.bam} | LC_COLLATE=C sort -k1,1 -k2,2n > {output.coverage}
        """

rule get_midpoint_coverage:
    input:
        tsv = lambda wc: expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["experimental_prefix"]) if wc.counttype=="counts" else expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["spikein_prefix"]),
        chrsizes = lambda wc: config["genome"]["chrsizes"] if wc.counttype=="counts" else config["genome"]["si-chrsizes"],
        plus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-minus.bedgraph"
    output:
        "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-midpoints.bedgraph"
    wildcard_constraints:
        counttype="counts|sicounts"
    shell: """
        half_median_fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4 | sort -k1,1n | awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]/2.0}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 4.0;}} }}' | xargs printf "%.*f\n" 0)
        bedtools unionbedg -i <(bedtools shift -i {input.plus} -g {input.chrsizes} -s $half_median_fragsize) <(bedtools shift -i {input.minus} -g {input.chrsizes} -s -$half_median_fragsize) -g {input.chrsizes} -empty | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4+$5}}' > {output}
        """

rule normalize:
    input:
        counts = "coverage/counts/{sample}_{factor}-chipnexus-counts-{strand}.bedgraph",
        midpoints = lambda wc: "coverage/counts/" + wc.sample + "_" + wc.factor + "-chipnexus-counts-midpoints.bedgraph" if wc.norm=="libsizenorm" else "coverage/sicounts/" + wc.sample + "_" + wc.factor + "-chipnexus-sicounts-midpoints.bedgraph"
    params:
        scalefactor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
    output:
        normalized = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph",
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus|protection|midpoints"
    log: "logs/normalize/normalize-{sample}-{norm}-{strand}.log"
    shell: """
        (bash scripts/libsizenorm.sh {input.midpoints} {input.counts} {params.scalefactor} > {output.normalized}) &> {log}
        """

rule macs2:
    input:
        bam = lambda wc: expand("alignment/{sample}-" + wc.species + "only.bam", sample=[k for k,v in PASSING.items() if v["group"]==wc.group]),
        chrsizes = lambda wc: config["genome"]["chrsizes"] if wc.species==config["combinedgenome"]["experimental_prefix"] else config["genome"]["si-chrsizes"],
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
        llocal = config["macs2"]["llocal"],
        qscore = config["macs2"]["fdr"]
    conda:
        "envs/macs2.yaml"
    log: "logs/macs2/macs2-{group}-{species}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --keep-dup all --bdg -n peakcalling/macs/{wildcards.group}-{wildcards.species} --SPMR --bw {params.bw} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript peakcalling/macs/{wildcards.group}-{wildcards.species}_model.r) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_peaks.narrowPeak) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_summits.bed) &>> {log}
        """

def selectchrom(wc):
    if wc.strand not in ["SENSE","ANTISENSE"]:
        if wc.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wc.norm=="sicounts":
        return os.path.splitext(config["genome"]["sichrsizes"])[0] + "-STRANDED.tsv"
    return os.path.splitext(config["genome"]["chrsizes"])[0] + "-STRANDED.tsv"

rule bedgraph_to_bigwig:
    input:
        bg = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph",
        chrsizes = selectchrom
    output:
        "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw"
    wildcard_constraints:
        strand="plus|minus|midpoints|protection|SENSE|ANTISENSE"
    shell: """
        bedGraphToBigWig {input.bg} {input.chrsizes} {output}
        """

rule get_si_pct:
    input:
        plmin = "coverage/counts/{sample}_{factor}-chipnexus-counts-midpoints.bedgraph",
        SIplmin = "coverage/sicounts/{sample}_{factor}-chipnexus-sicounts-midpoints.bedgraph"
    output:
        temp("qual_ctrl/all/{sample}_{factor}-spikeincounts.tsv")
    params:
        group = lambda wc: SAMPLES[wc.sample]["group"]
    log: "logs/get_si_pct/get_si_pct-{sample}_{factor}.log"
    shell: """
        (echo -e "{wildcards.sample}\t{params.group}\t" $(awk 'BEGIN{{FS=OFS="\t"; ex=0; si=0}}{{if(NR==FNR){{si+=$4}} else{{ex+=$4}}}} END{{print ex+si, ex, si}}' {input.SIplmin} {input.plmin}) > {output}) &> {log}
        """

rule cat_si_pct:
    input:
        expand("qual_ctrl/all/{sample}_{factor}-spikeincounts.tsv", sample=SAMPLES, factor=config["factor"])
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
        samplelist = lambda wc : [k for k,v in SAMPLES.items() if v["spikein"]=="y"] if wc.status=="all" else [k for k,v in PASSING.items() if v["spikein"]=="y"],
        conditions = conditiongroups_si,
        controls = controlgroups_si
    script: "scripts/plotsipct.R"

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

rule classify_peaks_genic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "_peaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/genic/{group}-exp-peaks-genic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.annotation} -wo | cut --complement -f11-13,15-17 > {output}
        """

rule classify_peaks_intragenic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "_peaks.narrowPeak",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/intragenic/{group}-exp-peaks-intragenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -wo | awk 'BEGIN{{FS=OFS="\t"}} $16=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $16=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
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

rule classify_peaks_intergenic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "_peaks.narrowPeak",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed"
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/intergenic/{group}-exp-peaks-intergenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} -wa -v | bedtools intersect -a stdin -b {input.genic_anno} -wa -v | bedtools intersect -a stdin -b {input.annotation} -wa > {output}
        """

rule peakstats:
    input:
        expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES),
    output:
        table = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peaknumbers.tsv",
        size = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peaksizes.svg",
        dist = "peakcalling/macs/{condition}-v-{control}-{factor}-chipnexus-peakdistances.svg"
    params:
        groups = lambda wc: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wc.condition=="all" else [wc.control, wc.condition],
        prefix = config["combinedgenome"]["experimental_prefix"]
    script:
        "scripts/peakstats.R"

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
        bg = lambda wc: "coverage/counts/" + wc.sample +"_"+ wc.factor+"-chipnexus-counts-midpoints.bedgraph" if wc.species==config["combinedgenome"]["experimental_prefix"] else "coverage/sicounts/" + wc.sample + "_" + wc.factor + "-chipnexus-sicounts-midpoints.bedgraph"
    output:
        temp("coverage/counts/{sample}_{factor}-{species}-peak-counts.tsv")
    shell: """
        bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}
        """

rule join_peak_counts:
    input:
        expand("coverage/counts/{sample}_{factor}-{{species}}-peak-counts.tsv", sample=SAMPLES, factor=config["factor"])
    output:
        "coverage/counts/union-bedgraph-allpeakcounts-{species}.tsv.gz"
    params:
        names = list(SAMPLES.keys())
    shell: """
        bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}
        """

rule call_de_peaks:
    input:
        expcounts = "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["experimental_prefix"] + ".tsv.gz",
        sicounts = lambda wc: "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["spikein_prefix"] + ".tsv.gz" if wc.norm=="spikenorm" else "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["experimental_prefix"] + ".tsv.gz"
    params:
        samples = lambda wc : [k for k,v in PASSING.items() if v["group"] in [wc.control, wc.condition]],
        groups = lambda wc: [v["group"] for k,v in PASSING.items() if v["group"] in [wc.control, wc.condition]],
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
        unchanged = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-unchanged.tsv",
        down = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-down.tsv"
    params:
        fdr = -log10(config["deseq"]["fdr"])
    log: "logs/separate_de_peaks/separate_de_peaks-{condition}-v-{control}-{norm}.log"
    shell: """
        (awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}} NR==1{{print > "{output.up}";print > "{output.unchanged}"; print > "{output.down}" }} NR>1 && $11>afdr && 11 != "NA" && $7>0 {{print > "{output.up}"}} NR>1 && $11>afdr && 11 != "NA" && $7<0 {{print > "{output.down}"}} NR>1 && ($11<=afdr || $11=="NA"){{print > "{output.unchanged}"}}' {input}) &> {log}
        """

rule de_peaks_to_bed:
    input:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.tsv",
    output:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}.bed",
    log: "logs/de_peaks_to_bed/de_peaks_to_bed-{condition}-v-{control}-{norm}-{direction}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $2, $4, $5, $1, $7":"$11, $3}}' {input} > {output}) &> {log}
        """

rule get_de_genic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv"
    output:
        "diff_binding/{condition}-v-{control}/genic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-genic.tsv"
    log : "logs/get_de_genic/get_de_genic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.annotation} -wo | awk 'BEGIN{{FS=OFS="\t"}} {{print $4, $8, $9, $10}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\ttranscript_start\ttranscript_end\ttranscript_name") - > {output}) &> {log}
        """

rule get_de_intragenic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.bed",
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv",
    output:
        "diff_binding/{condition}-v-{control}/intragenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intragenic.tsv",
    log: "logs/get_de_intragenic/get_de_intragenic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -wo | awk 'BEGIN{{FS=OFS="\t"}} $12=="+"{{print $4, $8, $9, $10, ((($2+1)+$3)/2)-$8}} $12=="-"{{print $4, $8, $9, $10, $9-((($2+1)+$3)/2)}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tORF_start\tORF_end\tORF_name\tdist_ATG_to_peak") - > {output}) &> {log}
        """

rule get_de_intergenic:
    input:
        peaks = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.bed",
        annotation = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed",
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed",
        totalresults = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv",
    output:
        "diff_binding/{condition}-v-{control}/intergenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intergenic.tsv",
    log : "logs/get_de_intergenic/get_de_intergenic-{condition}-v-{control}-{norm}.log"
    shell: """
        (bedtools intersect -a {input.peaks} -b {input.transcripts} {input.orfs} {input.genic_anno} -wa -v | bedtools intersect -a stdin -b {input.annotation} -wo | awk 'BEGIN{{FS=OFS="\t"}} {{print $4, $8, $9}}' | sort -k1,1 | join -t $'\t' <(tail -n +2 {input.totalresults} | cut --complement -f8-10 | sort -k1,1) - | sort -k8,8nr | cat <(echo -e "peak_name\tchrom\tstrand\tpeak_start\tpeak_end\tmeanExpr\tlog2FoldChange\tlogpadj\t{wildcards.condition}\t{wildcards.control}\tregion_start\tregion_end") - > {output}) &> {log}
        """

rule separate_sig_de:
    input:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-{category}.tsv"
    output:
        up = "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-up-{category}.tsv",
        unchanged = "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-unchanged-{category}.tsv",
        down = "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-down-{category}.tsv"
    params:
        fdr = -log10(config["deseq"]["fdr"])
    log: "logs/separate_sig_de/separate_sig_de-{condition}-v-{control}-{norm}-{category}.log"
    shell: """
        awk -v afdr={params.fdr} 'BEGIN{{FS=OFS="\t"}}NR==1{{print > "{output.up}";print > "{output.unchanged}"; print > "{output.down}"}} NR>1 && $7>0 && $8>afdr && $8 != "NA"{{print > "{output.up}"}} NR>1 && $7<0 && $8>afdr && $8 !="NA" {{print > "{output.down}"}} NR>1 && ($8<=afdr || $8 == "NA"){{print > "{output.unchanged}"}}' {input}
        """

rule get_de_category_bed:
    input:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $2, $4, $5, $1, $7":"$8, $3}}' {input} | sort -k1,1 -k2,2n  > {output}) &> {log}
        """

rule summarise_db_results:
    input:
        total = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all.tsv",
        genic = "diff_binding/{condition}-v-{control}/genic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-genic.tsv",
        intragenic = "diff_binding/{condition}-v-{control}/intragenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intragenic.tsv",
        intergenic = "diff_binding/{condition}-v-{control}/intergenic/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-all-intergenic.tsv",
    params:
        lfc = config["deseq"]["fold-change-threshold"],
        alpha = config["deseq"]["fdr"]
    output:
        summary = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-summary.svg",
        maplot = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-maplot.svg",
        volcano = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-volcano.svg",
        volcano_free = "diff_binding/{condition}-v-{control}/{condition}-v-{control}-{factor}-chipnexus-{norm}-diffbind-volcano-freescale.svg",
    script: "scripts/de_summary.R"

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
        plus = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-minus.bedgraph"
    output:
        sense = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus}> {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus}> {output.antisense}) &>> {log}
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


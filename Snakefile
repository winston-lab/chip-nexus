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

rule all:
    input:
        #FastQC
        'qual_ctrl/fastqc/per_base_sequence_content.svg',
        #alignment
        expand("alignment/{sample}-noPCRdup.bam", sample=SAMPLES),
        #coverage
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw", sample=SAMPLES, factor=config["factor"], norm=["counts","libsizenorm","spikenorm"], strand=["plus","minus","protection","midpoints"]),
        #initial QC
        "qual_ctrl/read_processing-loss.svg",
        expand("qual_ctrl/{status}/{status}-spikein-plots.svg", status=["all","passing"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-libsizenorm-correlations.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        expand(expand("qual_ctrl/{{status}}/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-{{status}}-window-{{windowsize}}-spikenorm-correlations.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), status=["all", "passing"], factor=config["factor"], windowsize=config["corr-windowsizes"]),
        #macs2
        expand("peakcalling/macs/{group}-{species}_peaks.narrowPeak", group = GROUPS, species=[config["combinedgenome"]["experimental_prefix"],config["combinedgenome"]["spikein_prefix"]]),
        #categorise peaks
        expand("peakcalling/macs/{category}/{group}-exp-peaks-{category}.tsv", group=GROUPS, category=CATEGORIES),
        expand(expand("peakcalling/macs/{condition}-v-{control}-{{factor}}-chipnexus-peaknumbers.tsv", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"]),
        # datavis
        expand(expand("datavis/{{annotation}}/libsizenorm/{condition}-v-{control}/{{factor}}-chipnexus-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        expand(expand("datavis/{{annotation}}/spikenorm/{condition}-v-{control}/{{factor}}-chipnexus-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}_heatmap-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        expand(expand("datavis/{{annotation}}/libsizenorm/{condition}-v-{control}/{{factor}}-chipnexus-{{annotation}}-libsizenorm-{{status}}_{condition}-v-{control}_stranded-metagene-bygroup.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        expand(expand("datavis/{{annotation}}/spikenorm/{condition}-v-{control}/{{factor}}-chipnexus-{{annotation}}-spikenorm-{{status}}_{condition}-v-{control}_stranded-metagene-bygroup.svg", zip, condition=conditiongroups_si+["all"], control=controlgroups_si+["all"]), annotation=config["annotations"], factor=config["factor"], status=["all","passing"]),
        #differential binding of peaks
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-libsizenorm.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-qcplots-spikenorm.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","unchanged","down"], fmt=["tsv", "bed"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}.{{fmt}}", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","unchanged", "down"], fmt=["tsv", "bed"]),
        #categorize DB peaks
        expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-libsizenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"], direction=["up","unchanged","down"], category=CATEGORIES),
        expand(expand("diff_binding/{condition}-v-{control}/{{category}}/{condition}-v-{control}-{{factor}}-chipnexus-results-spikenorm-{{direction}}-{{category}}.bed", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"], direction=["up","unchanged","down"], category=CATEGORIES),
        #DB summary
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-libsizenorm-diffbind-summary.svg", zip, condition=conditiongroups, control=controlgroups), factor=config["factor"]),
        expand(expand("diff_binding/{condition}-v-{control}/{condition}-v-{control}-{{factor}}-chipnexus-spikenorm-diffbind-summary.svg", zip, condition=conditiongroups_si, control=controlgroups_si), factor=config["factor"]),
        expand(expand("ratios/{{ratio}}/{condition}-v-{control}/{{factor}}-chipnexus-{{ratio}}_{{status}}_{condition}-v-{control}_violin.svg", zip, condition=conditiongroups+["all"], control=controlgroups+["all"]), factor=config["factor"], ratio=config["ratios"], status=["all", "passing"])

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

rule remove_adapter:
    input:
        lambda wildcards: SAMPLES[wildcards.sample]["fastq"]
    output:
        fq = temp("fastq/cleaned/{sample}-noadapter.fastq"),
        log =  "logs/remove_adapter/remove_adapter-{sample}.log"
    params:
        adapter = config["cutadapt"]["adapter"],
        trim_qual = config["cutadapt"]["trim_qual"]
    shell: """
       (less {input} | grep -B 1 -A 2 --no-group-separator '^.....CTGA' | cutadapt -e 0.1 -a {params.adapter} --nextseq-trim={params.trim_qual} -m 11 -o {output.fq} -) &> {output.log}
       """

rule remove_molec_barcode:
    input:
        "fastq/cleaned/{sample}-noadapter.fastq"
    output:
        fq = "fastq/cleaned/{sample}-cleaned.fastq.gz",
        barcodes = "qual_ctrl/molec_barcode/barcodes-{sample}.tsv",
        ligation = "qual_ctrl/molec_barcode/ligation-{sample}.tsv"
    threads: config["threads"]
    log: "logs/remove_molec_barcode/removeMBC-{sample}.log"
    shell: """
        (python scripts/extractNexusMolecularBarcode.py {input} fastq/cleaned/{wildcards.sample}-cleaned.fastq {output.barcodes} {output.ligation}) &> {log}
        (pigz -f fastq/cleaned/{wildcards.sample}-cleaned.fastq) &>> {log}
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
        kmer = 'qual_ctrl/fastqc/kmer_content.tsv'
    output:
        seq_len_dist = 'qual_ctrl/fastqc/sequence_length_distribution.svg',
        per_tile = 'qual_ctrl/fastqc/per_tile_quality.svg',
        per_base_qual = 'qual_ctrl/fastqc/per_base_quality.svg',
        per_base_seq = 'qual_ctrl/fastqc/per_base_sequence_content.svg',
        per_seq_gc = 'qual_ctrl/fastqc/per_sequence_gc.svg',
        per_seq_qual = 'qual_ctrl/fastqc/per_sequence_quality.svg',
        adapter_content = 'qual_ctrl/fastqc/adapter_content.svg',
        seq_dup = 'qual_ctrl/fastqc/sequence_duplication_levels.svg',
        kmer = 'qual_ctrl/fastqc/kmer_content.svg',
    script: "scripts/fastqc_summary.R"

rule bowtie2_build:
    input:
        fasta = config["combinedgenome"]["fasta"]
    output:
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/{{basename}}.rev.{num}.bt2", num=[1,2])
    params:
        idx_path = config["bowtie2"]["index-path"]
    log: "logs/bowtie2_build.log"
    shell: """
        (bowtie2-build {input.fasta} {params.idx_path}/{wildcards.basename}) &> {log}
        """

rule align:
    input:
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".{num}.bt2", num=[1,2,3,4]),
        expand(config["bowtie2"]["index-path"] + "/" + config["combinedgenome"]["name"] + ".rev.{num}.bt2", num=[1,2]),
        fastq = "fastq/cleaned/{sample}-cleaned.fastq.gz"
    output:
        bam = temp("alignment/{sample}-unique.bam"),
        aligned_fq = "fastq/aligned/{sample}-aligned.fastq.gz",
        unaligned_fq = "fastq/unaligned/{sample}-unaligned.fastq.gz",
        log = "logs/align/align-{sample}.log"
    params:
        outbase =  config["bowtie2"]["index-path"] + "/" +  config["combinedgenome"]["name"],
        minmapq = config["bowtie2"]["minmapq"]
    threads : config["threads"]
    shell: """
        (bowtie2 -x {params.outbase} -U {input.fastq} --al-gz {output.aligned_fq} --un-gz {output.unaligned_fq} -p {threads} | samtools view -buh -q {params.minmapq} - | samtools sort -T .{wildcards.sample} -@ {threads} -o {output.bam} -) &> {output.log}
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

#indexing is required for separating species by samtools view
rule index_bam:
    input:
        bam = "alignment/{sample}-noPCRdup.bam",
    output:
        "alignment/{sample}-noPCRdup.bam.bai"
    log : "logs/index_bam/index_bam-{sample}.log"
    shell: """
        (samtools index {input.bam}) &> {log}
        """

#peakcalling programs (MACS2 and Qnexus) require BAM input
rule bam_separate_species:
    input:
        bam = "alignment/{sample}-noPCRdup.bam",
        bai = "alignment/{sample}-noPCRdup.bam.bai",
        chrsizes = config["combinedgenome"]["chrsizes"]
    params:
        filterprefix = lambda wildcards: config["combinedgenome"]["spikein_prefix"] if wildcards.species==config["combinedgenome"]["experimental_prefix"] else config["combinedgenome"]["experimental_prefix"],
    output:
        "alignment/{sample}-{species}only.bam"
    threads: config["threads"]
    log: "logs/bam_separate_species/bam_separate_species-{sample}-{species}.log"
    shell: """
        (samtools view -h {input.bam} $(grep {wildcards.species} {input.chrsizes} | awk 'BEGIN{{FS="\t"; ORS=" "}}{{print $1}}') | grep -v -e 'SN:{params.filterprefix}' | sed 's/{wildcards.species}//g' | samtools view -bh -@ {threads} -o {output} -) &> {log}
        """

rule get_coverage:
    input:
        "alignment/{sample}-noPCRdup.bam"
    params:
        prefix = lambda wildcards: config["combinedgenome"]["experimental_prefix"] if wildcards.counttype=="counts" else config["combinedgenome"]["spikein_prefix"],
    output:
        # plmin = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-combined.bedgraph",
        plus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-minus.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/get_coverage/get_coverage-{sample}-{counttype}.log"
        # (genomeCoverageBed -bga -5 -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plmin}) &> {log}
    shell: """
        (genomeCoverageBed -bga -5 -strand + -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.plus}) &>> {log}
        (genomeCoverageBed -bga -5 -strand - -ibam {input} | grep {params.prefix} | sed 's/{params.prefix}//g' | sort -k1,1 -k2,2n > {output.minus}) &>> {log}
        """

# #NOTE: requires Q to be in $PATH
# rule qnexus:
#     input:
#         lambda wildcards: "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["experimental_prefix"] + "only.bam" if wildcards.counttype=="counts" else "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["spikein_prefix"] + "only.bam"
#     output:
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-narrowPeak.bed",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-qfrag-binding-characteristics.R",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-qfrag-binding-characteristics.pdf",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-quality-statistics.tab",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-runinfo.txt",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-summit-info.tab",
#         "peakcalling/qnexus/{sample}_{factor}-{counttype}-Q-treatment.bedgraph",
#         coverage = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-protection.bedgraph"
#     wildcard_constraints:
#         counttype="counts|sicounts"
#     log: "logs/qnexus/qnexus-{sample}-{counttype}.log"
#     shell: """
#         (Q --nexus-mode -t {input} -o peakcalling/qnexus/{wildcards.sample}_{wildcards.factor}-{wildcards.counttype} -v -wbt) &> {log}
#         (Rscript peakcalling/qnexus/{wildcards.sample}_{wildcards.factor}-{wildcards.counttype}-Q-qfrag-binding-characteristics.R) &>> {log}
#         (LC_COLLATE=C sort -k1,1 -k2,2n peakcalling/qnexus/{wildcards.sample}_{wildcards.factor}-{wildcards.counttype}-Q-treatment.bedgraph > {output.coverage}) &>> {log}
#         """

rule get_protection:
    input:
        tsv = lambda wildcards: expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["experimental_prefix"]) if wildcards.counttype=="counts" else expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["spikein_prefix"]),
        bam = lambda wildcards: "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["experimental_prefix"] +"only.bam" if wildcards.counttype=="counts" else "alignment/" + wildcards.sample + "-" + config["combinedgenome"]["spikein_prefix"] +"only.bam"
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
        tsv = lambda wildcards: expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["experimental_prefix"]) if wildcards.counttype=="counts" else expand("peakcalling/macs/{group}-{species}_peaks.xls", group=GROUPS, species=config["combinedgenome"]["spikein_prefix"]),
        chrsizes = lambda wildcards: config["genome"]["chrsizes"] if wildcards.counttype=="counts" else config["genome"]["si-chrsizes"],
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
        midpoints = lambda wildcards: "coverage/counts/" + wildcards.sample + "_" + wildcards.factor + "-chipnexus-counts-midpoints.bedgraph" if wildcards.norm=="libsizenorm" else "coverage/sicounts/" + wildcards.sample + "_" + wildcards.factor + "-chipnexus-sicounts-midpoints.bedgraph"
    params:
        scalefactor = lambda wildcards: config["spikein-pct"] if wildcards.norm=="spikenorm" else 1
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
        bam = lambda wildcards: expand("alignment/{sample}-" + wildcards.species + "only.bam", sample={k for (k,v) in PASSING.items() if (v["group"]==wildcards.group and v["pass-qc"]=="pass")}),
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
        qscore = config["macs2"]["fdr"]
    conda:
        "envs/macs2.yaml"
    log: "logs/macs2/macs2-{group}-{species}.log"
    shell: """
        (macs2 callpeak -t {input.bam} -f BAM -g $(awk '{{sum += $2}} END {{print sum}}' {input.chrsizes}) --keep-dup all --bdg -n peakcalling/macs/{wildcards.group}-{wildcards.species} --SPMR --bw {params.bw} --slocal {params.slocal} --llocal {params.llocal} --call-summits -q {params.qscore}) &> {log}
        (Rscript peakcalling/macs/{wildcards.group}-{wildcards.species}_model.r) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_peaks.narrowPeak) &>> {log}
        (sed -i -e 's/peakcalling\/macs\///g' peakcalling/macs/{wildcards.group}-{wildcards.species}_summits.bed) &>> {log}
        """

def selectchrom(wildcards):
    if wildcards.strand not in ["SENSE","ANTISENSE"]:
        if wildcards.norm=="sicounts":
            return config["genome"]["sichrsizes"]
        return config["genome"]["chrsizes"]
    if wildcards.norm=="sicounts":
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
        group = lambda wildcards: SAMPLES[wildcards.sample]["group"]
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
        samplelist = lambda wildcards : list({k:v for (k,v) in SAMPLES.items() if v["spikein"]=="y"}.keys()) if wildcards.status=="all" else list({k:v for (k,v) in PASSING.items() if v["spikein"]=="y"}.keys()),
        conditions = config["comparisons"]["spikenorm"]["conditions"],
        controls = config["comparisons"]["spikenorm"]["controls"],
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
        groups = lambda wildcards: [g for sublist in zip(controlgroups, conditiongroups) for g in sublist] if wildcards.condition=="all" else [wildcards.control, wildcards.condition],
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
        bg = lambda wildcards: "coverage/counts/" + wildcards.sample +"_"+ wildcards.factor+"-chipnexus-counts-midpoints.bedgraph" if wildcards.species==config["combinedgenome"]["experimental_prefix"] else "coverage/sicounts/" + wildcards.sample + "_" + wildcards.factor + "-chipnexus-sicounts-midpoints.bedgraph"
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
        sicounts = lambda wildcards: "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["spikein_prefix"] + ".tsv.gz" if wildcards.norm=="spikenorm" else "coverage/counts/union-bedgraph-allpeakcounts-" + config["combinedgenome"]["experimental_prefix"] + ".tsv.gz"
    params:
        samples = lambda wildcards : [k for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)],
        groups = lambda wildcards : [PASSING[x]["group"] for x in [k for (k,v) in PASSING.items() if (v["group"]== wildcards.control or v["group"]==wildcards.condition)]],
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
        lambda wildcards : config["annotations"][wildcards.annotation]["path"]
    output:
        "{annopath}/stranded/{annotation}-STRANDED.{ext}"
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
        pcount = lambda wildcards: 0.01*int(wildcards.windowsize),
        samplelist = plotcorrsamples
    script:
        "scripts/plotcorr.R"

rule compute_matrix:
    input:
        annotation = lambda wildcards: config["annotations"][wildcards.annotation]["path"] if wildcards.strand=="protection" else os.path.dirname(config["annotations"][wildcards.annotation]["path"]) + "/stranded/" + wildcards.annotation + "-STRANDED" + os.path.splitext(config["annotations"][wildcards.annotation]["path"])[1],
        bw = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw"
    output:
        dtfile = temp("datavis/{annotation}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}.mat.gz"),
        matrix = temp("datavis/{annotation}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}.tsv"),
        melted = "datavis/{annotation}/{norm}/{annotation}_{sample}_{factor}-chipnexus-{norm}-{strand}-melted.tsv.gz",
    params:
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        upstream = lambda wildcards: config["annotations"][wildcards.annotation]["upstream"] + config["annotations"][wildcards.annotation]["binsize"],
        dnstream = lambda wildcards: config["annotations"][wildcards.annotation]["dnstream"] + config["annotations"][wildcards.annotation]["binsize"],
        binsize = lambda wildcards: config["annotations"][wildcards.annotation]["binsize"],
        sort = lambda wildcards: config["annotations"][wildcards.annotation]["sort"],
        sortusing = lambda wildcards: config["annotations"][wildcards.annotation]["sortby"],
        binstat = lambda wildcards: config["annotations"][wildcards.annotation]["binstat"]
    threads : config["threads"]
    log: "logs/compute_matrix/compute_matrix-{annotation}_{sample}_{factor}-{norm}-{strand}.log"
    run:
        if config["annotations"][wildcards.annotation]["type"]=="absolute":
            refpoint = config["annotations"][wildcards.annotation]["refpoint"]
            if config["annotations"][wildcards.annotation]["nan_afterend"]=="y":
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --nanAfterEnd --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
            else:
                shell("""(computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        else:
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            refpoint = "TSS"
            shell("""(computeMatrix scale-regions -R {input.annotation} -S {input.bw} -out {output.dtfile} --outFileNameMatrix {output.matrix} -m {scaled_length} -b {params.upstream} -a {params.dnstream} --binSize {params.binsize} --sortRegions {params.sort} --sortUsing {params.sortusing} --averageTypeBins {params.binstat} -p {threads}) &> {log}""")
        melt_upstream = params.upstream-params.binsize
        shell("""(Rscript scripts/melt_matrix.R -i {output.matrix} -r {refpoint} --group {params.group} -s {wildcards.sample} -b {params.binsize} -u {melt_upstream} -o {output.melted}) &>> {log}""")

rule cat_matrices:
    input:
        expand("datavis/{{annotation}}/{{norm}}/{{annotation}}_{sample}_{{factor}}-chipnexus-{{norm}}-{{strand}}-melted.tsv.gz", sample=SAMPLES)
    output:
        "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-{strand}.tsv.gz"
    log: "logs/cat_matrices/cat_matrices-{annotation}_{norm}-{strand}.log"
    shell: """
        (cat {input} > {output}) &> {log}
        """

rule plot_heatmaps:
    input:
        matrix = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-protection.tsv.gz"
    output:
        heatmap_sample = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_heatmap-bysample.svg",
        heatmap_group = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_heatmap-bygroup.svg",
    params:
        samplelist = plotcorrsamples,
        mtype = lambda wildcards : config["annotations"][wildcards.annotation]["type"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        pct_cutoff = lambda wildcards : config["annotations"][wildcards.annotation]["pct_cutoff"],
        cluster = lambda wildcards : config["annotations"][wildcards.annotation]["cluster"],
        nclust = lambda wildcards: config["annotations"][wildcards.annotation]["nclusters"],
        heatmap_cmap = lambda wildcards : config["annotations"][wildcards.annotation]["heatmap_colormap"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_nexus_heatmaps.R -i {input.matrix} -s {params.samplelist} -t {params.mtype} -u {params.upstream} -d {params.dnstream} -c {params.pct_cutoff} -z {params.cluster} -k {params.nclust} -r {params.refpointlabel} -f {wildcards.factor} -l {scaled_length} -e {endlabel} -y {params.ylabel} -m {params.heatmap_cmap} -o {output.heatmap_sample} -p {output.heatmap_group}""")

rule plot_metagenes:
    input:
        plus = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-SENSE.tsv.gz",
        minus = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-ANTISENSE.tsv.gz",
        protection = "datavis/{annotation}/{norm}/allsamples-{annotation}-{factor}-chipnexus-{norm}-protection.tsv.gz"
    output:
        smeta_group = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_stranded-metagene-bygroup.svg",
        smeta_sample = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_stranded-metagene-bysample.svg",
        pmeta_group = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metagene-bygroup.svg",
        pmeta_sample = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metagene-bysample.svg",
        pmeta_goverlay = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metagene-overlay-group.svg",
        pmeta_soverlay = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metagene-overlay-sample.svg",
        pmeta_soverlay_bygroup = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metagene-overlay-sample-bygroup.svg",
        meta_heatmap_group = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metaheatmap-bygroup.svg",
        meta_heatmap_sample = "datavis/{annotation}/{norm}/{condition}-v-{control}/{factor}-chipnexus-{annotation}-{norm}-{status}_{condition}-v-{control}_protection-metaheatmap-bysample.svg",
    params:
        samplelist = plotcorrsamples,
        mtype = lambda wildcards : config["annotations"][wildcards.annotation]["type"],
        upstream = lambda wildcards : config["annotations"][wildcards.annotation]["upstream"],
        dnstream = lambda wildcards : config["annotations"][wildcards.annotation]["dnstream"],
        trim_pct = lambda wildcards : config["annotations"][wildcards.annotation]["trim_pct"],
        refpointlabel = lambda wildcards : config["annotations"][wildcards.annotation]["refpointlabel"],
        factor = config["factor"],
        ylabel = lambda wildcards : config["annotations"][wildcards.annotation]["ylabel"]
    run:
        if config["annotations"][wildcards.annotation]["type"]=="scaled":
            scaled_length = config["annotations"][wildcards.annotation]["scaled_length"]
            endlabel = config["annotations"][wildcards.annotation]["three_prime_label"]
        else:
            scaled_length=0
            endlabel = "HAIL SATAN!"
        shell("""Rscript scripts/plot_nexus_metagenes.R --inplus {input.plus} --inminus {input.minus} --inprotection {input.protection} -s {params.samplelist} -t {params.mtype} -u {params.upstream} -d {params.dnstream} -p {params.trim_pct} -r {params.refpointlabel} -f {wildcards.factor} -l {scaled_length} -e {endlabel} -y {params.ylabel} --out1 {output.smeta_group} --out2 {output.smeta_sample} --out3 {output.pmeta_group} --out4 {output.pmeta_sample} --out5 {output.pmeta_goverlay} --out6 {output.pmeta_soverlay} --out7 {output.pmeta_soverlay_bygroup} --out8 {output.meta_heatmap_group} --out9 {output.meta_heatmap_sample}""")

rule make_ratio_annotation:
    input:
        lambda wildcards: config["ratios"][wildcards.ratio]["path"]
    params:
        totalsize = lambda wildcards: config["ratios"][wildcards.ratio]["numerator"]["upstream"] + config["ratios"][wildcards.ratio]["numerator"]["dnstream"] + config["ratios"][wildcards.ratio]["denominator"]["upstream"] + config["ratios"][wildcards.ratio]["denominator"]["dnstream"],
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
        group = lambda wildcards : SAMPLES[wildcards.sample]["group"],
        upstream = lambda wildcards: config["ratios"][wildcards.ratio][wildcards.fractype]["upstream"],
        dnstream = lambda wildcards: config["ratios"][wildcards.ratio][wildcards.fractype]["dnstream"],
        refpoint = lambda wildcards: config["ratios"][wildcards.ratio][wildcards.fractype]["refpoint"]
    threads: config["threads"]
    log: "logs/ratio_counts/ratio_counts-{ratio}-{fractype}-{sample}.log"
    shell: """
        (computeMatrix reference-point -R {input.annotation} -S {input.bw} --referencePoint {params.refpoint} -out {output.dtfile} --outFileNameMatrix {output.matrix} -b {params.upstream} -a {params.dnstream} --binSize $(echo {params.upstream} + {params.dnstream} | bc) --averageTypeBins sum -p {threads}) &> {log}
        (Rscript scripts/melt_matrix.R -i {output.matrix} -r TSS --group {params.group} -s {wildcards.sample} -b $(echo {params.upstream} + {params.dnstream} | bc) -u {params.upstream} -o {output.melted}) &>> {log}
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

def ratiosamples(wildcards):
    dd = SAMPLES if wildcards.status=="all" else PASSING
    if wildcards.condition=="all":
        return list(dd.keys())
    else:
        return [k for k,v in dd.items() if v["group"]==wildcards.control or v["group"]==wildcards.condition]

rule plot_ratios:
    input:
        numerator = "ratios/{ratio}/allsamples_{ratio}_numerator.tsv.gz",
        denominator = "ratios/{ratio}/allsamples_{ratio}_denominator.tsv.gz",
    output:
        violin = "ratios/{ratio}/{condition}-v-{control}/{factor}-chipnexus-{ratio}_{status}_{condition}-v-{control}_violin.svg",
        ecdf = "ratios/{ratio}/{condition}-v-{control}/{factor}-chipnexus-{ratio}_{status}_{condition}-v-{control}_ecdf.svg"
    params:
        num_size = lambda wildcards: config["ratios"][wildcards.ratio]["numerator"]["upstream"] + config["ratios"][wildcards.ratio]["numerator"]["dnstream"],
        den_size = lambda wildcards: config["ratios"][wildcards.ratio]["denominator"]["upstream"] + config["ratios"][wildcards.ratio]["denominator"]["dnstream"],
        pcount = 1e-3,
        samplelist = ratiosamples,
        ratio_label = lambda wildcards: config["ratios"][wildcards.ratio]["ratio_name"],
        num_label = lambda wildcards: config["ratios"][wildcards.ratio]["numerator"]["region_label"],
        den_label = lambda wildcards: config["ratios"][wildcards.ratio]["denominator"]["region_label"],
        annotation_label = lambda wildcards: config["ratios"][wildcards.ratio]["label"]
    script:
        "scripts/ratio.R"

#!/usr/bin/env python

rule classify_peaks_genic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "_peaks.narrowPeak",
        annotation = build_annotations(os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed")
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
        genic_anno = build_annotations(os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed")
    params:
        prefix = config["combinedgenome"]["experimental_prefix"]
    output:
        "peakcalling/macs/intragenic/{group}-exp-peaks-intragenic.tsv"
    shell: """
        bedtools intersect -a {input.peaks} -b {input.genic_anno} -v | bedtools intersect -a stdin -b {input.orfs} -wo | awk 'BEGIN{{FS=OFS="\t"}} $16=="+"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $2+$10-$12}} $16=="-"{{print $1, $2, $3, $4, $5, $6, $7, $8, $9, $10, $14, $13-$2+$10}}' > {output}
        """

rule classify_peaks_intergenic:
    input:
        peaks = "peakcalling/macs/{group}-" + config["combinedgenome"]["experimental_prefix"] + "_peaks.narrowPeak",
        annotation = build_annotations(os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "intergenic-regions.bed"),
        transcripts = config["genome"]["transcripts"],
        orfs = config["genome"]["orf-annotation"],
        genic_anno = build_annotations(os.path.dirname(config["genome"]["transcripts"]) + "/" + config["combinedgenome"]["experimental_prefix"] + "genic-regions.bed")
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

rule get_de_category_bed:
    input:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.tsv"
    output:
        "diff_binding/{condition}-v-{control}/{category}/{condition}-v-{control}-{factor}-chipnexus-results-{norm}-{direction}-{category}.bed"
    log: "logs/get_category_bed/get_category_bed-{condition}-v-{control}-{norm}-{direction}-{category}.log"
    shell: """
        (awk 'BEGIN{{FS=OFS="\t"}} NR>1 {{print $2, $4, $5, $1, $7":"$8, $3}}' {input} | sort -k1,1 -k2,2n  > {output}) &> {log}
        """


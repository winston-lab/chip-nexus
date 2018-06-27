#!/usr/bin/env python

localrules:
    map_counts_to_peaks,
    combine_peak_counts

rule map_counts_to_peaks:
    input:
        bed = "diff_binding/{condition}-v-{control}/{condition}-v-{control}_{type}-{factor}-peaks.bed",
        bg = lambda wc: "coverage/counts/{sample}_{factor}-chipnexus-counts-midpoints.bedgraph".format(**wc) if wc.type=="experimental" else "coverage/sicounts/{sample}_{factor}-chipnexus-sicounts-midpoints.bedgraph".format(**wc)
    output:
        temp("diff_binding/{condition}-v-{control}/{sample}_{type}-{factor}-chipnexus-peakcounts.tsv")
    log: "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}-{sample}-{type}-{factor}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_peak_counts:
    input:
        lambda wc: ["diff_binding/{condition}-v-{control}/".format(**wc) + sample + "_{type}-{factor}-chipnexus-peakcounts.tsv".format(**wc) for sample in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_binding/{condition}-v-{control}/{condition}-v-{control}_allsamples-{type}-{factor}-chipnexus-allpeakcounts.tsv.gz"
    params:
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log: "logs/combine_peak_counts/combine_peak_counts_{condition}-v-{control}_{type}-{factor}.log"
    shell: """
        (bedtools unionbedg -i {input} -header -names {params.names} | bash scripts/cleanUnionbedg.sh | pigz -f > {output}) &> {log}
        """

rule differential_binding:
    input:
        expcounts = "diff_binding/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{factor}-chipnexus-allpeakcounts.tsv.gz",
        sicounts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{factor}-chipnexus-allpeakcounts.tsv.gz".format(**wc)
    output:
        results_all = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-all.tsv",
        results_up = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-up.tsv",
        results_down = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-down.tsv",
        results_unch = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-unchanged.tsv",
        normcounts = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-counts-sizefactornorm.tsv",
        rldcounts = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-counts-rlogtransformed.tsv",
        qcplots = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-qc-plots.svg"
    params:
        samples = lambda wc : get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc: [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["deseq"]["fdr"],
        lfc = log2(config["deseq"]["fold-change-threshold"]),
    script:
        "../scripts/call_diffbind_peaks.R"

#NOTE: coverage is smoothed using a Gaussian with bandwidth equal to the median
# fragment size across all samples as determined by MACS2 cross-correlation before
# finding the summit
rule diffbind_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-midpoints.bw", sample=get_samples("passing", wc.norm, wc.condition), norm=wc.norm, factor=FACTOR),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-midpoints.bw", sample=get_samples("passing", wc.norm, wc.control), norm=wc.norm, factor=FACTOR),
        diffbind_results = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-{direction}.tsv",
        tsv = lambda wc: expand("peakcalling/macs/{group}/{group}_experimental-{factor}-chipnexus_peaks.xls", factor=FACTOR, group=GROUPS),
    output:
        narrowpeak = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-{direction}.narrowpeak",
        summit_bed = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-{direction}-summits.bed",
    log: "logs/diffbind_results_to_narrowpeak/diffbind_results_to_narrowpeak-{condition}-v-{control}_{norm}-{direction}-{factor}.log"
    shell: """
        median_fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4 | sort -k1,1n | awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 2.0;}} }}' | xargs printf "%.*f\n" 0)
        (python scripts/diffbind_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -s $median_fragsize -d {input.diffbind_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """

rule summarise_db_results:
    input:
        total = "diff_binding/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-diffbind-results-all.tsv",
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
    script: "../scripts/plot_diffbind_summary.R"


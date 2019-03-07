#!/usr/bin/env python

localrules:
    map_counts_to_annotations,
    combine_annotation_counts

rule map_counts_to_annotations:
    input:
        bed = lambda wc: "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_{species}-{factor}-peaks.bed" if wc.annotation=="peaks" else config["differential_occupancy"]["annotations"][wc.annotation],
        bg = lambda wc: "coverage/counts/{sample}_{factor}-chipnexus-counts-midpoints.bedgraph".format(**wc) if wc.species=="experimental" else "coverage/sicounts/{sample}_{factor}-chipnexus-sicounts-midpoints.bedgraph".format(**wc)
    output:
        temp("diff_binding/{annotation}/{condition}-v-{control}/{sample}_{species}-{factor}-chipnexus-counts-{annotation}.tsv")
    log:
        "logs/map_counts_to_peaks/map_counts_to_peaks-{condition}-v-{control}-{sample}-{species}-{annotation}-{factor}.log"
    shell: """
        (bedtools map -a {input.bed} -b {input.bg} -c 4 -o sum > {output}) &> {log}
        """

rule combine_annotation_counts:
    input:
        lambda wc: ["diff_binding/{annotation}/{condition}-v-{control}/".format(**wc) + sample + "_{species}-{factor}-chipnexus-counts-{annotation}.tsv".format(**wc) for sample in get_samples("passing", "libsizenorm", [wc.control, wc.condition])]
    output:
        "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-{species}-{factor}-chipnexus-counts-{annotation}.tsv.gz"
    params:
        n = lambda wc: 7*len(get_samples("passing", "libsizenorm", [wc.control, wc.condition])),
        names = lambda wc: "\t".join(get_samples("passing", "libsizenorm", [wc.control, wc.condition]))
    log:
        "logs/combine_peak_counts/combine_peak_counts_{condition}-v-{control}_{species}-{annotation}-{factor}.log"
    shell: """
        (paste {input} | \
         cut -f$(paste -d, <(echo "1-6") <(seq -s, 7 7 {params.n})) | \
         cat <(echo -e "chrom\tstart\tend\tname\tscore\tstrand\t{params.names}" ) - > {output}) &> {log}
        """

rule differential_binding:
    input:
        exp_counts = "diff_binding/{annotation}/{condition}-v-{control}/{condition}-v-{control}_allsamples-experimental-{factor}-chipnexus-counts-{annotation}.tsv.gz",
        spike_counts = lambda wc: [] if wc.norm=="libsizenorm" else "diff_binding/peaks/{condition}-v-{control}/{condition}-v-{control}_allsamples-spikein-{factor}-chipnexus-counts-peaks.tsv.gz".format(**wc)
    output:
        results_all = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-all.tsv",
        results_up = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-up.tsv",
        results_down = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-down.tsv",
        results_unchanged = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-unchanged.tsv",
        counts_norm = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-counts-sizefactornorm.tsv",
        counts_rlog = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-counts-rlogtransform.tsv",
        qc_plots = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-qc-plots.svg"
    params:
        samples = lambda wc: get_samples("passing", wc.norm, [wc.control, wc.condition]),
        groups = lambda wc: [PASSING[x]["group"] for x in get_samples("passing", wc.norm, [wc.control, wc.condition])],
        alpha = config["differential_occupancy"]["fdr"],
        lfc = log2(config["differential_occupancy"]["fold_change_threshold"]),
    conda:
        "../envs/diff_bind.yaml"
    script:
        "../scripts/differential_binding.R"

#NOTE: coverage is smoothed using a Gaussian with bandwidth equal to the median
# fragment size across all passing samples as determined by MACS2 cross-correlation before
# finding the summit
rule diffbind_results_to_narrowpeak:
    input:
        condition_coverage = lambda wc: expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-midpoints.bw", sample=get_samples("passing", wc.norm, wc.condition), norm=wc.norm, factor=FACTOR),
        control_coverage = lambda wc: expand("coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-midpoints.bw", sample=get_samples("passing", wc.norm, wc.control), norm=wc.norm, factor=FACTOR),
        diffbind_results = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-{direction}.tsv",
        tsv = lambda wc: expand("peakcalling/sample_peaks/{sample}_experimental-{factor}-chipnexus_peaks.xls", sample=PASSING, factor=FACTOR),
    output:
        narrowpeak = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-{direction}.narrowpeak",
        summit_bed = "diff_binding/{annotation}/{condition}-v-{control}/{norm}/{condition}-v-{control}_{factor}-chipnexus-{norm}-{annotation}-diffbind-results-{direction}-summits.bed",
    conda:
        "../envs/diffbind_to_narrowpeak.yaml"
    log:
        "logs/diffbind_results_to_narrowpeak/diffbind_results_to_narrowpeak-{condition}-v-{control}_{norm}-{annotation}-{direction}-{factor}.log"
    shell: """
        median_fragsize=$(grep -e "^# d = " {input.tsv} | \
                cut -d ' ' -f4 | \
                sort -k1,1n | \
                awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 2.0;}} }}' | \
                xargs printf "%.*f\n" 0)
        (python scripts/diffbind_results_to_narrowpeak.py -i {input.condition_coverage} -j {input.control_coverage} -s $median_fragsize -d {input.diffbind_results} -n {output.narrowpeak} -b {output.summit_bed}) &> {log}
        """


#!/usr/bin/env python

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


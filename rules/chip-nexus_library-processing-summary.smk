#!/usr/bin/env python

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


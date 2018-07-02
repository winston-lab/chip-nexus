#!/usr/bin/env python

rule crosslink_coverage:
    input:
        lambda wc: "alignment/{sample}_{factor}-chipnexus-noPCRduplicates-".format(**wc) + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
    output:
        "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-{strand}.bedgraph",
    params:
        strand_symbol = lambda wc: {"plus": "+", "minus": "-"}.get(wc.strand)
    wildcard_constraints:
        counttype="counts|sicounts",
        strand="plus|minus"
    log: "logs/crosslink_coverage/crosslink_coverage-{sample}-{counttype}-{strand}-{factor}.log"
    shell: """
        (bedtools genomecov -bga -5 -strand {params.strand_symbol} -ibam {input} > {output}) &> {log}
        """

#extend reads to the median fragment size over all samples as
#called by MACS2 cross-correlation
rule protection_coverage:
    input:
        tsv = lambda wc: expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.xls", factor=FACTOR, group=GROUPS, species= ("experimental" if wc.counttype=="counts" else "spikein")),
        bam = lambda wc: "alignment/{sample}_{factor}-chipnexus-noPCRduplicates-".format(**wc) + ("experimental" if wc.counttype=="counts" else "spikein") + ".bam",
    output:
        "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-protection.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/genome_coverage/genome_coverage-{sample}-{counttype}-protection-{factor}.log"
    shell: """
        median_fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4 | sort -k1,1n | awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 2.0;}} }}' | xargs printf "%.*f\n" 0)
        (bedtools genomecov -bga -fs $median_fragsize -scale $(echo 1/$median_fragsize | bc -l) -ibam {input.bam} | LC_COLLATE=C sort -k1,1 -k2,2n > {output}) &> {log}
        """

#shift 5' ends of reads by half the median fragment size
#over all samples
rule midpoint_coverage:
    input:
        tsv = lambda wc: expand("peakcalling/macs/{group}/{group}_{species}-{factor}-chipnexus_peaks.xls", factor=FACTOR, group=GROUPS, species= ("experimental" if wc.counttype=="counts" else "spikein")),
        chrsizes = lambda wc: config["genome"]["chrsizes"] if wc.counttype=="counts" else config["genome"]["sichrsizes"],
        plus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-plus.bedgraph",
        minus = "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-minus.bedgraph"
    output:
        "coverage/{counttype}/{sample}_{factor}-chipnexus-{counttype}-midpoints.bedgraph",
    wildcard_constraints:
        counttype="counts|sicounts"
    log: "logs/genome_coverage/genome_coverage-{sample}-{counttype}-midpoints-{factor}.log"
    shell: """
        half_median_fragsize=$(grep -e "^# d = " {input.tsv} | cut -d ' ' -f4 | sort -k1,1n | awk '{{count[NR]=$1;}} END{{if (NR % 2) {{print count[(NR+1)/2]/2.0}} else {{print (count[(NR/2)] + count[(NR/2)+1]) / 4.0;}} }}' | xargs printf "%.*f\n" 0)
        (bedtools unionbedg -i <(bedtools shift -i {input.plus} -g {input.chrsizes} -s $half_median_fragsize) <(bedtools shift -i {input.minus} -g {input.chrsizes} -s -$half_median_fragsize) -g {input.chrsizes} -empty | awk 'BEGIN{{FS=OFS="\t"}}{{print $1, $2, $3, $4+$5}}' > {output}) &> {log}
        """

rule normalize_genome_coverage:
    input:
        counts = "coverage/counts/{sample}_{factor}-chipnexus-counts-{strand}.bedgraph",
        bam = lambda wc: "alignment/{sample}_{factor}-chipnexus-noPCRduplicates-".format(**wc) + ("experimental" if wc.norm=="libsizenorm" else "spikein") + ".bam",
    output:
        normalized = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bedgraph",
    params:
        scale_factor = lambda wc: config["spikein-pct"] if wc.norm=="spikenorm" else 1
    wildcard_constraints:
        norm="libsizenorm|spikenorm",
        strand="plus|minus|protection|midpoints"
    conda: "../envs/default.yaml"
    log: "logs/normalize_genome_coverage/normalize_genome_coverage-{sample}-{norm}-{strand}-{factor}.log"
    shell: """
        (awk -v norm_factor=$(samtools view -c {input.bam} | paste -d "" - <(echo "/({params.scale_factor}*1000000)") | bc -l) 'BEGIN{{FS=OFS="\t"}}{{$4=$4/norm_factor; print $0}}' {input.counts} > {output.normalized}) &> {log}
        """

rule make_stranded_bedgraph:
    input:
        plus = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-plus.bedgraph",
        minus = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-minus.bedgraph"
    output:
        sense = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-SENSE.bedgraph",
        antisense = "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-ANTISENSE.bedgraph"
    log : "logs/make_stranded_bedgraph/make_stranded_bedgraph-{sample}-{norm}-{factor}.log"
    shell: """
        (bash scripts/makeStrandedBedgraph.sh {input.plus} {input.minus}> {output.sense}) &> {log}
        (bash scripts/makeStrandedBedgraph.sh {input.minus} {input.plus}> {output.antisense}) &>> {log}
        """

def select_chromsizes(wc):
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
        chrsizes = select_chromsizes
    output:
        "coverage/{norm}/{sample}_{factor}-chipnexus-{norm}-{strand}.bw"
    log : "logs/bedgraph_to_bigwig/bedgraph_to_bigwig-{sample}-{norm}-{strand}-{factor}.log"
    wildcard_constraints:
        strand="plus|minus|midpoints|protection|SENSE|ANTISENSE"
    shell: """
        (bedGraphToBigWig {input.bg} {input.chrsizes} {output}) &> {log}
        """


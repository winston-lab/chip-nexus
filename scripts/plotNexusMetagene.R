library(tidyverse)
library(forcats)
library(viridis)

import = function(path){
    read_table2(path,
    	 col_names=c("group", "sample", "index", "position","cpm"),
    	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double())) %>%
    filter(cpm != "NA") %>%
    return()
}

plotmeta = function(intable.plus, intable.minus, refptlabel, factor, ylabel, samples_out, group_out){
    raw = intable.plus %>% import() %>% right_join(intable.minus %>% import, by=c("group","sample","index","position"))
    raw$cpm.y = -raw$cpm.y
    raw$sample = factor(raw$sample, ordered = TRUE)
    raw$group = factor(raw$group, ordered = TRUE)
    nindices = max(raw$index)
    nsamples = length(fct_unique(raw$sample))
    ngroups = length(fct_unique(raw$group))
    w = round((max(raw$position) - min(raw$position))*1000/148) 
    df = raw %>% group_by(group, sample, position) %>% summarise(pos.mean = mean(cpm.x, na.rm=TRUE), neg.mean = mean(cpm.y, na.rm=TRUE))
    
    metagene_base = ggplot(data = df, aes(x=position)) +
                    geom_vline(xintercept = 0, size=1, color="black", linetype="dashed") +
                    geom_ribbon(aes(ymax = pos.mean, ymin = neg.mean), fill="#781B6C") +
                    theme_minimal() +
                    xlab(paste("distance from", refptlabel, "(kb)")) +
                    ylab(paste("Average", factor, "ChIP-nexus coverage,\n over", nindices, ylabel)) +
                    theme(strip.text = element_text(size=12, face="bold"),
                          axis.text.x = element_text(size=12, face="bold"),
                          axis.text.y = element_text(size=12),
                          axis.title = element_text(size=12, face="bold"))
    
    metagene_group = metagene_base + facet_wrap(~group, ncol=ngroups)
    metagene_samples = metagene_base + facet_wrap(~sample, dir="v", ncol=ngroups) 
    ggsave(group_out, plot = metagene_group, height= 10, width = 5+2*w*ngroups, units = "cm")
    ggsave(samples_out, plot = metagene_samples, height= 10+5*(nsamples/ngroups), width = 5+2*w*ngroups, units = "cm")
}

plotmeta(intable.plus = snakemake@input[["plus"]],
             intable.minus = snakemake@input[["minus"]],
             #upstream = snakemake@params[["upstream"]],
             #downstream= snakemake@params[["downstream"]],
             refptlabel = snakemake@params[["refpointlabel"]],
             factor = snakemake@params[["factor"]],
             ylabel = snakemake@params[["ylabel"]],
             samples_out = snakemake@output[["meta_sample"]],
             group_out = snakemake@output[["meta_group"]])

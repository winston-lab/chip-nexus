library(tidyverse)
library(forcats)
library(viridis)

axislabel = function(up, dn, xlab){
    if (up < (1/3)*dn){
        return( c('', xlab, (dn/1000)))
    }
    else if (dn < (1/3)*up){
        return( c((-up/1000), xlab, ''))
    }
    else {
        return( c((-up/1000), xlab, (dn/1000)))
    }
}

plotheatmaps = function(intable, samplelist, upstream, dnstream, pct_cutoff, cluster, k, refptlab, factor, ylabel, cmap, samples_out, group_out){
    raw = read_tsv(intable, col_names=c("group", "sample", "index", "position","cpm")) %>%
            filter(sample %in% samplelist & !is.na(cpm)) 
    raw$sample = fct_inorder(raw$sample, ordered = TRUE)
    raw$group = fct_inorder(raw$group, ordered = TRUE)
    
    nindices = max(raw$index, na.rm=TRUE)
    nsamples = length(fct_unique(raw$sample))
    ngroups = length(fct_unique(raw$group))
    
    #clustering
    if (cluster){
        rr = raw %>% select(-group) %>% unite(cid, c(sample, position), sep="~") %>%
                spread(cid, cpm, fill=0) %>% select(-index)
        clust = kmeans(rr, centers = k)$cluster %>% as_tibble() %>% rename(cluster=value) %>%
                mutate_at(vars(cluster), as.factor) %>% mutate(og_index=row_number()) %>%
                arrange(cluster,og_index) %>% mutate(new_index=row_number())
        raw = raw %>% left_join(clust, by=c("index"="og_index")) %>% select(-index) %>% rename(index=new_index)
    } 
    
    #percentile cutoff for heatmap visualization
    cutoff = quantile(raw$cpm, probs=pct_cutoff, na.rm=TRUE)
    
    #plot heatmap facetted by sample and group
    heatmap_base = ggplot(data = raw %>% mutate_at(vars(cpm), funs(pmin(cutoff, .)))) +
      geom_raster(aes(x=position, y=index, fill=cpm)) +
      scale_y_reverse(name=paste(nindices, ylabel), expand=c(0.01, 0))+
      scale_x_continuous(breaks = c(-upstream/1000, 0, dnstream/1000),
                         labels= axislabel(up=upstream, dn=dnstream, xlab=refptlab), 
                         minor_breaks = scales::pretty_breaks(n=10),
                         name=paste("distance from", refptlab, "(kb)")) +
      scale_fill_viridis(option = cmap, na.value="FFFFFF00", name=paste(factor, 'ChIP-nexus signal'), guide=guide_colorbar(title.position="top", barwidth=15, barheight=1, title.hjust=0.5)) +
      theme_minimal() +
      theme(text = element_text(size=12, face="bold", color="black"),
              legend.position = "top",
              legend.title = element_text(size=12, face="bold", color="black"),
              legend.text = element_text(size=8, face="plain"),
              strip.text = element_text(size=12, face="bold", color="black"),
              axis.text.y = element_blank(),
              axis.text.x = element_text(size=12, face="bold", color="black", margin = unit(c(0,0,0,0),"cm")),
              panel.grid.major.x = element_line(color="black"),
              panel.grid.minor.x = element_line(color="grey80"),
              panel.grid.major.y = element_line(color="grey80"),
              panel.grid.minor.y = element_blank(),
              panel.spacing.x = unit(.5, "cm"))

    hmap.width = max(12, (.0008*(upstream+dnstream)+3.4)*ngroups)
    
    heatmap_samples = heatmap_base + facet_wrap(~sample, dir="v", ncol=ngroups)
    ggsave(samples_out, plot = heatmap_samples, height=(.0005*nindices+7.5), width = hmap.width, units = "cm", limitsize=FALSE)
    rm(heatmap_samples)
    heatmap_groups = heatmap_base + facet_wrap(~group, ncol=ngroups)
    ggsave(group_out, plot = heatmap_groups, height= .0009*nindices+11.5, width = hmap.width, units = "cm")
}

plotheatmaps(intable= snakemake@input[["matrix"]],
             samplelist = snakemake@params[["samplelist"]],
             upstream = snakemake@params[["upstream"]],
             dnstream= snakemake@params[["dnstream"]],
             pct_cutoff = snakemake@params[["pct_cutoff"]],
             cluster = snakemake@params[["cluster"]],
             k = snakemake@params[["nclust"]],
             refptlab = snakemake@params[["refpointlabel"]],
             factor = snakemake@params[["factor"]],
             ylabel = snakemake@params[["ylabel"]],
             cmap = snakemake@params[["heatmap_cmap"]],
             samples_out = snakemake@output[["heatmap_sample"]],
             group_out = snakemake@output[["heatmap_group"]])

library(tidyverse)
library(forcats)
#library(modeest)
library(gridExtra)

main = function(groups, factor, prefix, table.out, size.out, dist.out){
    dflist = list()
    groups = unique(groups)
    for (i in 1:length(groups)){
        g = groups[i]
        dflist[[g]][['all']] = read_tsv(paste0('peakcalling/macs/', groups[i], '-', prefix, '_peaks.narrowPeak'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10)
        dflist[[g]][['genic']] = read_tsv(paste0('peakcalling/macs/genic/', groups[i], '-exp-peaks-genic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11)
        dflist[[g]][['intragenic']] = read_tsv(paste0('peakcalling/macs/intragenic/', groups[i], '-exp-peaks-intragenic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10, gene.name=X11, peak.dist.to.ATG=X12)
        dflist[[g]][['intergenic']] = read_tsv(paste0('peakcalling/macs/intergenic/', groups[i], '-exp-peaks-intergenic.tsv'), col_names=FALSE) %>%
                                        select(chrom=X1, start=X2, end=X3, strand=X6, summit=X10)
    }

    ndf = tibble()
    for (i in 1:length(groups)){
        g = groups[i]
        dd = tibble(group=g,
                    all=nrow(dflist[[g]][['all']]), genic=nrow(dflist[[g]][['genic']]),
                    intragenic=nrow(dflist[[g]][['intragenic']]),
                    intergenic=nrow(dflist[[g]][['intergenic']]))
        ndf = ndf %>% bind_rows(dd)
    }

    binddfs = function(type){
        df = tibble()
        for (i in 1:length(groups)){
            g = groups[i]
            dd = dflist[[g]][[type]] %>% mutate(group=g)
            df = df %>% bind_rows(dd)
        }
        df$group = fct_inorder(df$group, ordered=TRUE)
        return(df)
    }

    alldf = binddfs('all')
    genicdf = binddfs('genic')
    intragenicdf = binddfs('intragenic')
    intergenicdf = binddfs('intergenic')

    sizedf = alldf %>% transmute(size=end-start, group=group, type='all') %>%
                bind_rows(genicdf %>% transmute(size=end-start, group=group, type='genic')) %>%
                bind_rows(intragenicdf %>% transmute(size=end-start, group=group, type='intragenic')) %>%
                bind_rows(intergenicdf %>% transmute(size=end-start, group=group, type='intergenic'))
    sizedf$type = fct_inorder(sizedf$type, ordered=TRUE)
    sizedf$group = fct_inorder(sizedf$group, ordered=TRUE)

    sizedf %>% group_by(group, type) %>% summarise(n= n()) %>% spread(key=type, value=n) %>% ungroup() %>%
        write_tsv(table.out)

    sizeannodf = sizedf %>% group_by(type, group, size) %>% mutate(freq = n()) %>%
                    group_by(type) %>% mutate(y = max(freq)) %>%
                    group_by(group, type) %>%
                    summarise(size = .7*quantile(sizedf$size, .9993), y = .5*unique(y), n = n())

    sizehist = ggplot() +
                geom_histogram(data = sizedf, aes(size), binwidth=1, fill="#08306b") +
                geom_text(data = sizeannodf, aes(x=size, y = y, label = n), hjust=1, size=4, fontface="bold") +
                facet_grid(type~group, scales="free_y", space="free_y", switch="y") +
                scale_x_continuous(limits = c(NA, quantile(sizedf$size, .9993)), name=paste(factor, "ChIP-nexus peak size (bp)")) +
                scale_y_continuous(position="right", name=NULL, breaks=scales::pretty_breaks(n=4)) +
                theme_light() +
                theme(text = element_text(size=12, color="black", face="bold"),
                      axis.text = element_text(size=10, color="black"),
                      strip.text = element_text(size=12, color="black", face="bold"),
                      strip.text.y = element_text(angle=180, hjust=1),
                      strip.background = element_blank(),
                      panel.grid.major = element_line(color="grey80"),
                      panel.grid.minor = element_line(color="grey80"))

    ggsave(size.out, plot = sizehist, width = length(groups)*6, height = 12, units = "cm", limitsize=FALSE)

    intraannodf = intragenicdf %>% select(group, dist = peak.dist.to.ATG) %>% mutate(x = .95*quantile(dist, .995)) %>% group_by(group) %>%
                #summarize(mode = mlv(dist, bw=100, method='parzen', kernel="gaussian")[['M']][[1]], median = median(dist), max=max(dist), x = first(x), n=n())
                summarize(median = median(dist), max=max(dist), x = first(x), n=n())
    intradistplot = ggplot() +
                    geom_histogram(data = intragenicdf, aes(peak.dist.to.ATG), binwidth=25, fill="#3f007d") +
                    geom_text(data = intraannodf,
                              aes(x=x, label = paste0("n= ", n, "\nmedian= ",median, "\nmax= ", max)), y=40,
                              hjust=1, size=3, fontface="bold") +
                    scale_x_continuous(limits = c(NA, quantile(intragenicdf$peak.dist.to.ATG, .995)),
                                       name=paste("distance from ATG to intragenic", factor, "ChIP-nexus peak (bp)"),
                                       minor_breaks = scales::pretty_breaks(n=10)) +
                    scale_y_continuous(name=NULL) +
                    facet_grid(~group) +
                    theme_light() +
                    theme(text = element_text(size=12, color="black", face="bold"),
                          axis.text = element_text(size=10, color="black"),
                          strip.text = element_text(size=12, color="black", face="bold"),
                          strip.text.y = element_text(angle=180, hjust=1),
                          strip.background = element_blank(),
                          panel.grid.major = element_line(color="grey80"),
                          panel.grid.minor = element_line(color="grey80"))

    ggsave(dist.out, plot = intradistplot, width = length(groups)*6, height=12, units="cm")
}

main(groups = snakemake@params[["groups"]],
     factor = snakemake@wildcards[["factor"]],
     prefix = snakemake@params[["prefix"]],
     table.out = snakemake@output[["table"]],
     size.out = snakemake@output[["size"]],
     dist.out = snakemake@output[["dist"]])

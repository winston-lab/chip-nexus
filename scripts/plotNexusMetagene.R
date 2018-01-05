library(psych)
library(tidyverse)
library(forcats)
library(ggthemes)

import = function(path){
    read_tsv(path,
    	 col_names=c("group", "sample", "index", "position","cpm"),
    	 col_types=cols(group=col_character(), sample=col_character(), index=col_integer(), position=col_double(), cpm=col_double())) %>%
    filter(cpm != "NA") %>%
    return()
}

format_xaxis_kb = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x))
}

format_xaxis_nt = function(refptlabel){
    function(x) if_else(x==0, refptlabel, as.character(x*1000))
}

label_xaxis= function(ggp, refptlabel, upstream, downstream){
    if(upstream>1000 | downstream>1000){
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_kb(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(kb)"),
                               limits = c(-upstream/1000, downstream/1000),
                               expand=c(0,0))
    } else{
        ggp = ggp +
            scale_x_continuous(breaks=scales::pretty_breaks(n=3),
                               labels=format_xaxis_nt(refptlabel=refptlabel),
                               name=paste("distance from", refptlabel, "(nt)"),
                               limits = c(-upstream/1000, downstream/1000),
                               expand=c(0,0))
    } 
    return(ggp)
}

apply_theme = function(ggp){
    ggp = ggp +
        theme_light() +
        theme(strip.text = element_text(size=12, face="bold", color="black"),
              strip.text.y = element_text(angle=0),
              axis.text.x = element_text(size=10, face="bold", color="black"),
              axis.text.y = element_text(size=10, color="black"),
              axis.title = element_text(size=12, face="bold"),
              panel.grid.major.x = element_line(color="grey60"),
              panel.grid.minor.x = element_line(color="grey80"),
              panel.grid.major.y = element_line(color="grey80"),
              panel.grid.minor.y = element_line(color="grey80"),
              panel.spacing.x = unit(0.5, "cm"),
              strip.background=element_blank(),
              plot.title = element_text(size=12, face="bold"),
              plot.subtitle = element_text(size=10))
    return(ggp)
}

stranded_meta = function(df, factor, nindices, ylabel, upstream, downstream, refptlabel){
    metagene_base = ggplot(data=df, aes(x=position)) +
                    geom_vline(xintercept=0, size=1, color="grey65") +
                    geom_col(aes(y=pos.mean), fill="#114477", color="#114477", alpha=.9, size=0.1) +
                    geom_col(aes(y=neg.mean), fill="#4477AA", color="#4477AA", alpha=.9, size=0.1) +
                    ggtitle(paste("mean", factor, "ChIP-nexus coverage"),
                            subtitle=paste(nindices, ylabel)) +
                    ylab("normalized counts")
    metagene_base %>% apply_theme() %>%
        label_xaxis(upstream=upstream, downstream=downstream, refptlabel=refptlabel) %>%
        return()
}

protection_meta = function(df, factor, nindices, ylabel, upstream, downstream, refptlabel){
    metagene_base = ggplot(data = df, aes(x=position)) +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                    fill="#114477", alpha=0.4, size=0) +
        geom_line(aes(y=mean), color="#114477") +
        scale_y_continuous(limits=c(0, NA), name="normalized counts") +
        ggtitle(paste("mean", factor, "protection"),
                subtitle = paste(nindices, ylabel))
    metagene_base %>% apply_theme() %>%
        label_xaxis(upstream=upstream, downstream=downstream, refptlabel=refptlabel) %>%
        return()
} 

main = function(intable.plus, intable.minus, intable.qfrags, samplelist, upstream, downstream,
                trim_pct, refptlabel, factor, ylabel, smeta_group_out, smeta_sample_out,
                pmeta_group_out, pmeta_sample_out, pmeta_goverlay_out,
                pmeta_soverlay_out, pmeta_soverlay_bygroup_out){
    stranded = intable.plus %>% import() %>%
        right_join(intable.minus %>% import,
                   by=c("group","sample","index","position")) %>% 
        filter(sample %in% samplelist) %>% 
        mutate_at(vars(cpm.y), funs(-.)) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
    
    repl_df = stranded %>% select(group, sample) %>% distinct() %>%
        group_by(group) %>% mutate(replicate=row_number()) %>% ungroup() %>% 
        select(-group)
    
    nindices = max(stranded$index, na.rm=TRUE)
    nsamples = length(fct_unique(stranded$sample))
    ngroups = length(fct_unique(stranded$group))
    
    sdf_group = stranded %>% group_by(group, position) %>%
        summarise(pos.mean = winsor.mean(cpm.x, trim=trim_pct, na.rm=TRUE),
                  neg.mean = winsor.mean(cpm.y, trim=trim_pct, na.rm=TRUE))
    smeta_group = stranded_meta(sdf_group, factor=factor, nindices=nindices,
                                ylabel=ylabel, upstream=upstream,
                                downstream=downstream, refptlabel=refptlabel) +
        facet_wrap(~group, ncol=ngroups)
    ggsave(smeta_group_out, plot=smeta_group, height=8, width=7*ngroups, units="cm")
    
    sdf_sample = stranded %>% left_join(repl_df, by="sample") %>%
        group_by(group, sample, replicate, position) %>%
        summarise(pos.mean = winsor.mean(cpm.x, trim=trim_pct, na.rm=TRUE),
                  neg.mean = winsor.mean(cpm.y, trim=trim_pct, na.rm=TRUE))
    smeta_sample = stranded_meta(sdf_sample, factor=factor, nindices=nindices, ylabel=ylabel, upstream=upstream,
                                downstream=downstream, refptlabel=refptlabel) +
        facet_grid(replicate~group)
    ggsave(smeta_sample_out, plot=smeta_sample, height=8+5*(nsamples/ngroups-1),
           width=7*ngroups, units="cm")
    
    protection = import(intable.qfrags) %>%
        filter(sample %in% samplelist) %>% 
        mutate_at(vars(sample, group), funs(fct_inorder(., ordered=TRUE)))
    
    pdf_group = protection %>% group_by(group, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct, na.rm=TRUE),
                  sd = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE),
                  sem = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE)/sqrt(n()))
    pmeta_group = protection_meta(pdf_group, factor=factor, nindices=nindices, ylabel=ylabel, upstream=upstream,
                                downstream=downstream, refptlabel=refptlabel) +
        facet_wrap(~group, ncol=ngroups)
    ggsave(pmeta_group_out, plot=pmeta_group, height=8, width=7*ngroups, units="cm")
    
    pdf_sample = protection %>% left_join(repl_df, by="sample") %>% 
        group_by(group, sample, replicate, position) %>%
        summarise(mean = winsor.mean(cpm, trim=trim_pct, na.rm=TRUE),
                  sd = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE),
                  sem = winsor.sd(cpm, trim=trim_pct, na.rm=TRUE)/sqrt(n()))
    pmeta_sample = protection_meta(pdf_sample, factor=factor, nindices=nindices, ylabel=ylabel, upstream=upstream,
                                downstream=downstream, refptlabel=refptlabel) +
        facet_grid(replicate~group)
    ggsave(pmeta_sample_out, plot=pmeta_sample, height=2+4.5*(max(repl_df$replicate)),
           width=7*ngroups, units="cm")
    
    pmeta_goverlay = ggplot(data=pdf_group, aes(x=position, color=group, fill=group)) +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                        alpha=0.3, size=0) +
        geom_line(aes(y=mean)) +
        scale_fill_ptol() + scale_color_ptol() +
        scale_y_continuous(limits=c(0, NA),
                           name="normalized coverage") +
        ggtitle(paste("mean", factor, "protection"),
                subtitle=paste(nindices, ylabel)) 
    pmeta_goverlay = pmeta_goverlay %>% apply_theme() %>%
        label_xaxis(upstream=upstream, downstream=downstream, refptlabel=refptlabel) +
        theme(legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold"))
    ggsave(pmeta_goverlay_out, plot=pmeta_goverlay, width=14, height=8, units="cm")
    
    pmeta_soverlay = ggplot(data=pdf_sample,
                            aes(x=position, group=sample, color=group, fill=group)) +
        geom_vline(xintercept = 0, size=1, color="grey65") +
        geom_ribbon(aes(ymax=mean+1.96*sem, ymin=mean-1.96*sem),
                        alpha=0.2, size=0) +
        geom_line(aes(y=mean)) +
        scale_fill_ptol() + scale_color_ptol() +
        scale_y_continuous(limits=c(0, NA),
                           name="normalized coverage") +
        ggtitle(paste("mean", factor, "protection"),
                subtitle=paste(nindices, ylabel)) 
    pmeta_soverlay = pmeta_soverlay %>% apply_theme() %>%
        label_xaxis(upstream=upstream, downstream=downstream, refptlabel=refptlabel) +
        theme(legend.title = element_blank(),
              legend.text = element_text(size=12, face="bold"))
    ggsave(pmeta_soverlay_out, plot=pmeta_soverlay, width=14, height=8, units="cm")
    
    pmeta_soverlay_bygroup = pmeta_soverlay +
        facet_grid(group~., switch="y") +
        theme(legend.position="none",
              strip.placement="outside",
              strip.text.y=element_text(angle=180, hjust=1))
    ggsave(pmeta_soverlay_bygroup_out, plot=pmeta_soverlay_bygroup,
           width=14, height=2+4.5*(ngroups), units="cm")
}

main(intable.plus = snakemake@input[["plus"]],
     intable.minus = snakemake@input[["minus"]],
     intable.qfrags = snakemake@input[["qfrags"]],
     samplelist = snakemake@params[["samplelist"]],
     upstream = snakemake@params[["upstream"]],
     downstream= snakemake@params[["dnstream"]],
     trim_pct = snakemake@params[["trim_pct"]],
     refptlabel = snakemake@params[["refpointlabel"]],
     factor = snakemake@params[["factor"]],
     ylabel = snakemake@params[["ylabel"]],
     smeta_group_out = snakemake@output[["smeta_group"]],
     smeta_sample_out = snakemake@output[["smeta_sample"]],
     pmeta_group_out = snakemake@output[["pmeta_group"]],
     pmeta_sample_out = snakemake@output[["pmeta_sample"]],
     pmeta_goverlay_out = snakemake@output[["pmeta_goverlay"]],
     pmeta_soverlay_out = snakemake@output[["pmeta_soverlay"]],
     pmeta_soverlay_bygroup_out = snakemake@output[["pmeta_soverlay_bygroup"]])

library(dplyr)
library(stringr)
library(ggplot2)
library(grid)

# load oldref results
old_ref_results = read.delim("data/ASM680v1_dura3_vs_drnr.tsv")
old_ref_results_plasmids = old_ref_results[grep("NC_002607.1",
                                                rownames(old_ref_results),
                                                invert = T),]

# load newref results
new_ref_results = read.delim("TK_dura3_vs_drnr.tsv")
new_ref_results_plasmids = new_ref_results[grep("SC_chr",
                                       rownames(new_ref_results),
                                       invert = T),]
# load genes in deletions
deleted_genes = read.delim("deleted-genes_drnr-minus-dura3.gff",
                           header=F)


# load orthogroups
orthogroups = read.delim("Orthogroups.tsv")


# 1) Extracting protein ID from Sig plasmidial genes (roldref)
#oldref_prot_ids = 
#        str_extract(rownames(old_ref_results),'DAC[0-9]*[.]?[1]?') %>%
#        na.omit()

# 1B) Extracting protein ID from deleted genes (roldref)
deleted_prot_ids = 
        str_extract(deleted_genes$V9,'DAC[0-9]*[.]?[1]?') %>%
        na.omit()

# 2) Grepping it with orthogoups
grep_results = orthogroups[grep(paste(deleted_prot_ids, collapse = "|"), 
                                orthogroups$PFEIFFER_annot),]

# 3) parsing result
df_deleted_orthogroups = 
        data.frame(Oldref = grep_results$PFEIFFER_annot, 
                   Newref = grep_results$SC_annot)

# 4) grep those on deleted)_orthogroups against Newref results
single_values = 
        df_deleted_orthogroups$Newref[grep(",", 
                                           df_deleted_orthogroups$Newref,
                                           invert = T)] 
single_values = single_values[single_values %>% str_detect("_")]

mult_values = 
        df_deleted_orthogroups$Newref[grep(",", 
                                           df_deleted_orthogroups$Newref,
                                           invert = F)]
Newref_deleted_orthogroups = mult_values %>% 
        str_split(", ") %>% 
        unlist() %>% 
        c(single_values)

# 5) check p-values of those who match
# df_to_plot must be all genes that has at least one ortholog

single_values =
        orthogroups$SC_annot[grep(",", orthogroups$SC_annot, invert = T)] 

single_values = single_values[single_values %>% str_detect("_")]

mult_values = 
        orthogroups$SC_annot[grep(",", orthogroups$SC_annot, invert = F)]

all_orthologs_in_newref_results = 
        mult_values %>% str_split(", ") %>% unlist() %>% c(single_values)


df_to_plot = new_ref_results_plasmids
df_to_plot$Status = "Presentes"

df_to_plot[grep(paste(Newref_deleted_orthogroups, collapse = "|"),
                              rownames(df_to_plot)),7] = "Deletados"

# doing the same for the oldref data
df_to_plot_new = old_ref_results_plasmids
df_to_plot_new$Status = "Presentes"
df_to_plot_new[grep(paste(deleted_prot_ids, collapse = "|"), 
                    rownames(df_to_plot_new)),7] = "Deletados"


df_to_plot_new$Reference = "ASM680v1"
df_to_plot$Reference = "TK"

col = c("#e64b35", "#4dbbd5")

rbind(df_to_plot_new, df_to_plot) %>%
        ggplot(aes(x=padj,col=Status, fill=Status)) +
        geom_density(alpha=0.2, size = 1.1) +
        theme_bw() + labs(y="Densidade", x = "p-valor ajustado das análises EGD") + 
        scale_fill_manual(values = col) + 
        scale_color_manual(values = col ) +
        theme(text = element_text(size=15), legend.position = c(0.92,0.85),
              legend.title = element_blank()) +
        facet_wrap(. ~ Reference)

# Writing file
rbind(df_to_plot_new, df_to_plot) %>%
        write.csv(file = "rnr_BothReferences_DGE.csv")
        

library(argparse)
library(BiocParallel)
library(MPRAnalyze)
library(tidyverse)

parser = ArgumentParser()
parser$add_argument('all_rep_counts', help = 'counts from all reps')
parser$add_argument('basename', help = 'output_basename')
args = parser$parse_args()

get_colannos = function(current_colnames){
    
    data.frame(col_ids = current_colnames) %>% 
        separate(col_ids, remove = FALSE, sep = "\\.", into = c('batch', 'condition')) %>%
        column_to_rownames('col_ids')
}

all_rep_counts = read_rds(args$all_rep_counts)

all_mpranalyse_input = all_rep_counts %>% 
    pivot_wider(id_cols = gBC, names_from = condition, names_sep = '.', values_from = colnames(select(., contains('count')))) %>%
    column_to_rownames('gBC') %>% 
    replace(is.na(.), 0) 

all_mpranalyse_DNA = all_mpranalyse_input %>% 
    select(contains('DNA')) %>% 
    rename_with(~ sub('DNA_count.', '', .x)) %>% 
    as.matrix()

all_mpranalyse_RNA = all_mpranalyse_input %>% 
    select(contains('RNA')) %>% 
    rename_with(~ sub('RNA_count.', '', .x)) %>% 
    as.matrix()

all_mpranalyse_colanno = get_colannos(colnames(all_mpranalyse_DNA)) %>% 
    mutate(batch = as.factor(batch),
           condition = as.factor(condition)) 

obj = MpraObject(dnaCounts = all_mpranalyse_DNA, rnaCounts = all_mpranalyse_RNA, 
              dnaAnnot = all_mpranalyse_colanno, rnaAnnot = all_mpranalyse_colanno,
              BPPARAM = MulticoreParam())

obj = estimateDepthFactors(obj, lib.factor = c("batch", "condition"),
                        which.lib = "both")

obj = analyzeQuantification(obj = obj, 
                            dnaDesign = ~ batch + condition,
                            rnaDesign = ~ condition)

mpranalyse_alphas = getAlpha(obj, by.factor = "condition")

write_rds(mpranalyse_alphas, paste0(args$basename, '_mpranalyze_alphas.rds'))
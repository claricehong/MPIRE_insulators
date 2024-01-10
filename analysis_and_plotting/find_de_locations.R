library(argparse)
library(BiocParallel)
library(MPRAnalyze)
library(tidyverse)

parser = ArgumentParser()
parser$add_argument('all_rep_counts', help = 'counts from all reps')
parser$add_argument('iBC', help = 'iBC to analyse')
parser$add_argument('basename', help = 'basename to add to output file')
parser$add_argument('-m', help = 'whether to compare mutants')
args = parser$parse_args()

get_colannos = function(current_colnames){
    
    data.frame(col_ids = current_colnames) %>% 
        separate(col_ids, remove = FALSE, sep = "\\.", into = c('batch', 'condition')) %>%
        column_to_rownames('col_ids')
}

format_mpranalyse = function(selected_iBC){
    
    if (args$m != 'none'){
    both = all_rep_counts %>%
        filter(condition %in% c(selected_iBC, args$m)) %>%
        pivot_wider(id_cols = gBC, names_from = condition, names_sep = '.', values_from = colnames(select(., contains('count')))) %>%     
        column_to_rownames('gBC') %>%
        drop_na()
        }
    else {
    both = all_rep_counts %>% 
        filter(condition == selected_iBC | condition == 'empty') %>% 
        pivot_wider(id_cols = gBC, names_from = condition, names_sep = '.', values_from = colnames(select(., contains('count')))) %>% 
        column_to_rownames('gBC') %>% 
        drop_na()
    }
    
    DNA = both %>% select(contains('DNA')) %>% rename_with(~ sub('DNA_count.', '', .x)) 
    
    RNA = both %>% select(contains('RNA')) %>% rename_with(~ sub('RNA_count.', '', .x)) 
    
    if (args$m != 'none'){
    colanno = get_colannos(colnames(DNA)) %>%
	mutate(batch = as.factor(batch),
	       condition = as.factor(condition)) %>%
        mutate(condition = factor(condition, levels = c(selected_iBC, args$m)))
    }
    else {
    colanno = get_colannos(colnames(DNA)) %>%
        mutate(batch = as.factor(batch),
              condition = as.factor(condition)) %>%
        mutate(condition = factor(condition, levels = c('empty', selected_iBC)))
    }

    list(DNA = as.matrix(DNA), colanno = colanno, RNA = as.matrix(RNA)) 
}

mpra_compare = function(current_dfs){
    
    obj = MpraObject(dnaCounts = current_dfs$DNA, rnaCounts = current_dfs$RNA, 
                  dnaAnnot = current_dfs$colanno, rnaAnnot = current_dfs$colanno,
                  BPPARAM = MulticoreParam())
    
    obj = estimateDepthFactors(obj, lib.factor = c("batch", "condition"),
                            which.lib = "both")
    
    obj = analyzeComparative(obj = obj, 
                           dnaDesign = ~ batch + condition, 
                           rnaDesign = ~ condition, 
                           reducedDesign = ~ 1)

    testLrt(obj)
}

all_rep_counts = read_rds(args$all_rep_counts) %>% 
    mutate(across(contains('count'), function(x){x+1}))

formatted = format_mpranalyse(args$iBC)
compared = mpra_compare(formatted)

write_rds(compared, paste(args$iBC, args$basename, 'mpranalyze_de.rds', sep = '_'))

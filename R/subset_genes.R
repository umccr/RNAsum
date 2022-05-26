##### Function used in the "sv_prioritize" function
subset_genes = function(genes, ind) {
  genes %>% str_split('&') %>% map(~ .[ind] %>% replace("", NA) %>% .[!is.na(.)]) %>% map_chr(~ ifelse(length(.) > 0, str_c(., collapse = '&'), ""))
}

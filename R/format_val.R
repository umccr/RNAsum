##### Function used in the "sv_prioritize" function
format_val = function(val, is_pct = F) {
  ifelse(!is.na(val),
         format(val,  digits = 1) %>% str_c(ifelse(is_pct, "%", "")), NA)
}

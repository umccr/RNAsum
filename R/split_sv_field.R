##### Function used in the "sv_prioritize" function
split_sv_field = function(.data, field, is_pct = F) {
  f_q = rlang::enquo(field)
  f_str = rlang::quo_name(f_q)
  f1_str = str_c(f_str, '1')
  f2_str = str_c(f_str, '2')
  f1_q = sym(f1_str)
  f2_q = sym(f2_str)
  .data %>%
    separate(!!f_q, c(f1_str, f2_str), ",") %>%
    dplyr::mutate(
      !!f1_q := as.double(!!f1_q) * ifelse(is_pct, 100, 1),
      !!f2_q := as.double(!!f2_q) * ifelse(is_pct, 100, 1),
      !!f_q  := (!!f1_q + ifelse(is.na(!!f2_q), !!f1_q, !!f2_q)) / 2,
      !!f_q  := format_val(!!f_q, is_pct),
      !!f1_q := format_val(!!f1_q, is_pct),
      !!f2_q := format_val(!!f2_q, is_pct)
    )
}

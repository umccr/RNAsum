# Not in 
"%!in%" <- function(x,table) {
    match(x,table, nomatch = 0) == 0
}


## From rapportools - https://github.com/Rapporter/rapportools/blob/f45730af9cbdf147cafbd7c030602bae2fe915d5/R/utils.R#L356
is.empty <- function(x, trim = TRUE, ...) {
    vgsub <- function(pattern, replacement, x, ...){
        for(i in 1:length(pattern))
            x <- gsub(pattern[i], replacement[i], x, ...)
        x
    }

    trim.space <- function(x, what = c('both', 'leading', 'trailing', 'none'), space.regex = '[:space:]', ...){
        if (missing(x))
            stop('nothing to trim spaces to =(')
        re <- switch(match.arg(what),
                     both     = sprintf('^[%s]+|[%s]+$', space.regex, space.regex),
                     leading  = sprintf('^[%s]+', space.regex),
                     trailing = sprintf('[%s]+$', space.regex),
                     none     = {
                         return (x)
                     })
        vgsub(re, '', x, ...)
    }

    if (length(x) <= 1) {
        if (is.null(x))
            return (TRUE)
        if (length(x) == 0)
            return (TRUE)
        if (is.na(x) || is.nan(x))
            return (TRUE)
        if (is.character(x) && nchar(ifelse(trim, trim.space(x), x)) == 0)
            return (TRUE)
        if (is.logical(x) && !isTRUE(x))
            return (TRUE)
        if (is.numeric(x) && x == 0)
            return (TRUE)
        return (FALSE)
    } else {
        return(sapply(x, is.empty, trim = trim, ...))
    }
}

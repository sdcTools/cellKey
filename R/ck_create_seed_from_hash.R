ck_create_seed_from_hash <- function(x) {
  hash <- sha1(x)
  hash <- unlist(strsplit(hash, ""))
  hash <- sapply(hash, function(x) {
    res <- suppressWarnings(as.numeric(x))
    ifelse(is.na(res), grep(x, letters), x)
  })

  n <- as.numeric(paste0("0.", paste(hash, collapse = "")))
  as.integer(n * .Machine$integer.max)
}

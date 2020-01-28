ck_create_seed_from_hash <- function(x) {
  hash <- digest::sha1(x)
  hash <- unlist(strsplit(hash, ""))
  hash <- sapply(hash, function(x) {
    res <- suppressWarnings(as.numeric(x))
    ifelse(is.na(res), grep(x, letters), x)
  })
  as.integer(as.numeric(paste0("0.", paste(hash, collapse = ""))) * .Machine$integer.max)
}

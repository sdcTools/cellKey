args = commandArgs(trailingOnly=TRUE)
if (nchar(args)>0) {
  for (p in strsplit(args,";")[[1]]) {
    message("Installation von Paket",shQuote(p),"\n")
    install.packages(p)
  }
} else {
  message("Es werden keine zusätzlichen Pakete manuell installiert!")
}
message("CHECK wird gestartet!")
library(devtools)
devtools::install(path = ".")
checkOutput <- devtools::check()
if(length(checkOutput$errors)!=0){
  stop("CHECK sind nicht erfolgreich durchgelaufen")
}else{
  if(length(checkOutput$warnings)>0){
    warning("Es wurden Warnings in dem Check gefunden:")
    print(checkOutput$warnings)
  }else{
    message("CHECK lief ohne Errors und Warnings durch.")
  }
}
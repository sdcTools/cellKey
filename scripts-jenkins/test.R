library(devtools)
if (dir.exists("tests") || dir.exists("test")) {
  message("Tests starten!")
  testOutput <- devtools::test()
  #TODO: testOutput weiter verarbeiten ?!?!
  message("Tests beendet!")
} else {
  message("Keine Tests definiert!")
}
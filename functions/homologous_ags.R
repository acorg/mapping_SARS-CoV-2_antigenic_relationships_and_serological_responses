
homologous_ags <- c(
  "D614G" = "D614G",
  "B.1.1.7" = "B.1.1.7",
  "P.1" = "P.1",
  "B.1.351" = "B.1.351",
  "B.1.526+E484K" = "B.1.526+E484K",
  "B.1.617.2" = "B.1.617.2",
  "B.1.637" = "B.1.429", # Closest in sequence
  "C.37" = "C.37",
  "BA.1" = "BA.1",
  "BA.2" = "BA.2",
  "2x mRNA-1273" = "D614G",
  "3x mRNA-1273 BD01" = "D614G",
  "3x mRNA-1273 BD29" = "D614G",
  "3x mRNA-1273 (6 month)" = "D614G",
  "2x mRNA-1273.351" = "B.1.351"
)

srHomologousLogTiters <- function(map) {
  logtiters <- logtiterTable(map)
  homologous_ags <- srHomologousAgs(map)
  lapply(seq_len(numSera(map)), function(i) {
    logtiters[homologous_ags[[i]], i]
  })
}


sr_group_labels_no_sera <- c(
  "D614G"                  = "D614G",
  "B.1.1.7"                = "B.1.1.7 (Alpha)",
  "P.1"                    = "P.1 (Gamma)",
  "B.1.351"                = "B.1.351 (Beta)",
  "B.1.526+E484K"          = "B.1.526+E484K (Iota)",
  "B.1.617.2"              = "B.1.617.2 (Delta)",
  "B.1.637"                = "B.1.637",
  "C.37"                   = "C.37 (Lambda)",
  "BA.1"                   = "BA.1 (Omicron)",
  "BA.2"                   = "BA.2 (Omicron)",
  "2x mRNA-1273"           = "4 weeks post 2x mRNA-1273",
  "3x mRNA-1273 BD01"      = ">3 months post 2x mRNA-1273",
  "3x mRNA-1273 BD29"      = "4 weeks post 3x mRNA-1273",
  "3x mRNA-1273 (6 month)" = ">3 months post 3x mRNA-1273",
  "2x mRNA-1273.351"       = "4 weeks post 2x mRNA-1273.351"
)

sr_group_labels_multiline <- c(
  "D614G"                  = "D614G",
  "B.1.1.7"                = "B.1.1.7",
  "P.1"                    = "P.1",
  "B.1.351"                = "B.1.351",
  "B.1.526+E484K"          = "B.1.526+E484K",
  "B.1.617.2"              = "B.1.617.2",
  "B.1.637"                = "B.1.637",
  "C.37"                   = "C.37",
  "BA.1"                   = "BA.1",
  "BA.2"                   = "BA.2",
  "2x mRNA-1273"           = "4 weeks post\n2x mRNA-1273",
  "3x mRNA-1273 BD01"      = ">3 months post\n2x mRNA-1273",
  "3x mRNA-1273 BD29"      = "4 weeks post\n3x mRNA-1273",
  "3x mRNA-1273 (6 month)" = ">3 months post\n3x mRNA-1273",
  "2x mRNA-1273.351"       = "4 weeks post\n2x mRNA-1273.351"
)

sr_group_labels_multiline_with_sera <- sr_group_labels_multiline
sr_group_labels_multiline_with_sera[1:10] <- paste(sr_group_labels_multiline_with_sera[1:10], "sera")

sr_group_labels <- paste(sr_group_labels_no_sera, "sera")
names(sr_group_labels) <- names(sr_group_labels_no_sera)

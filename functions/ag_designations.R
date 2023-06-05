
ag_designations <- c(
  "D614G"                  = "",
  "B.1.1.7"                = "Alpha",
  "B.1.1.7+E484K"          = "",
  "P.1"                    = "Gamma",
  "B.1.351"                = "Beta",
  "B.1.429"                = "Epsilon",
  "B.1.526+E484K"          = "Iota",
  "B.1.617.1"              = "Kappa",
  "B.1.617.2"              = "Delta",
  "B.1.617.2+K417N"        = "",
  "B.1.617.2 (AY.1)+K417N" = "",
  "B.1.617.2 (AY.2)+K417N" = "",
  "B.1.617.2 (AY.3)+E484Q" = "",
  "B.1.621"                = "Mu",
  "C.37"                   = "Lambda",
  "BA.1"                   = "Omicron",
  "BA.1.1"                 = "Omicron",
  "BA.2"                   = "Omicron",
  "BA.2.12.1"              = "Omicron",
  "BA.3"                   = "Omicron",
  "BA.4/BA.5"              = "Omicron"
)

add_ag_designation <- function(ag_names) {
 
  ag_name_designations <- ag_designations[ag_names]
  ag_names[ag_name_designations != ""] <- sprintf(
    "%s (%s)", 
    ag_names[ag_name_designations != ""],
    ag_name_designations[ag_name_designations != ""]
  )
  ag_names
   
}

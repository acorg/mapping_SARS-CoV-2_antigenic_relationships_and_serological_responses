
# Setup workspace
rm(list = ls())
library(tidyverse)
library(Racmacs)
set.seed(100)

source("functions/homologous_ags.R")
source("functions/adjust_titers_by_slope.R")
source("functions/set_pt_drawing_order.R")

# Read in the raw titers and sequence data
longtiters <- readr::read_csv("data/raw_data/titerdata.csv", col_types = cols(.default = "c"))
seq_info <- readxl::read_excel("data/raw_data/sequences.xlsx")

# Take average titer for each serum and antigen
longtiters %>%
  group_by(
    ag, SampleID, sr_group, sample_source
  ) %>%
  summarise(
    ID50_merged = paste(ID50, collapse = ", "),
    ID50 = Racmacs:::ac_merge_titers(ID50, options = RacMerge.options(sd_limit = NA, dilution_stepsize = 0, method = "likelihood")),
    .groups = "drop"
  ) -> longtiters

# Convert the data to a map file
longtiters %>%
  select(
    SampleID,
    sr_group,
    sample_source,
    ag,
    ID50
  ) %>%
  pivot_wider(
    names_from = "ag",
    values_from = "ID50",
    values_fill = "*"
  ) -> titertabledata

titertable <- as.matrix(titertabledata[,-(1:3)])
rownames(titertable) <- unname(unlist(titertabledata[,1]))

# Create map
map_full <- acmap(titer_table = t(titertable))
srGroups(map_full) <- factor(titertabledata$sr_group)
srExtra(map_full) <- titertabledata$sample_source
dilutionStepsize(map_full) <- 0

# Set antigen groups
ag_groups <- c(
  "D614G" = "wildtype",
  "B.1.1.7" = "wildtype",
  "B.1.1.7+E484K" = "wildtype",
  "P.1" = "wildtype",
  "B.1.351" = "wildtype",
  "B.1.429" = "wildtype",
  "B.1.526+E484K" = "wildtype",
  "B.1.617.1" = "wildtype",
  "B.1.617.2" = "wildtype",
  "B.1.617.2+K417N" = "wildtype",
  "B.1.617.2 (AY.1)+K417N" = "wildtype",
  "B.1.617.2 (AY.2)+K417N" = "wildtype",
  "B.1.617.2 (AY.3)+E484Q" = "wildtype",
  "B.1.621" = "wildtype",
  "C.37" = "wildtype",
  "BA.1" = "wildtype",
  "BA.1.1" = "wildtype",
  "BA.2" = "wildtype",
  "BA.2.12.1" = "wildtype",
  "BA.3" = "wildtype",
  "BA.4/BA.5" = "wildtype",
  "D614G+E484K" = "d614g_mutant",
  "D614G+E484Q" = "d614g_mutant",
  "D614G+L452R" = "d614g_mutant",
  "D614G+L452R+E484Q" = "d614g_mutant",
  "D614G+N501Y" = "d614g_mutant",
  "P.1+T417K" = "417_mutant",
  "B.1.429+K417N" = "417_mutant",
  "B.1.351+N417K" = "417_mutant",
  "B.1.617.1+K417N" = "417_mutant",
  "BA.1+A484K" = "BA1_mutant"
)

agGroups(map_full) <- factor(
  ag_groups[agNames(map_full)], 
  c("wildtype", "d614g_mutant", "417_mutant", "BA1_mutant")
)

# Reorder antigens
map_full <- orderAntigens(
  map_full,
  c(
    "D614G",
    "B.1.1.7",
    "B.1.1.7+E484K",
    "P.1",
    "B.1.351",
    "B.1.429",
    "B.1.526+E484K",
    "B.1.617.1",
    "B.1.617.2",
    "B.1.617.2+K417N",
    "B.1.617.2 (AY.1)+K417N",
    "B.1.617.2 (AY.2)+K417N",
    "B.1.617.2 (AY.3)+E484Q",
    "B.1.621",
    "C.37",
    "BA.1",
    "BA.1.1",
    "BA.2",
    "BA.2.12.1",
    "BA.3",
    "BA.4/BA.5",
    "D614G+E484K",
    "D614G+E484Q",
    "D614G+L452R",
    "D614G+L452R+E484Q",
    "D614G+N501Y",
    "P.1+T417K",
    "B.1.429+K417N",
    "B.1.351+N417K",
    "B.1.617.1+K417N",
    "BA.1+A484K"
  )
)

# Reorder serum groups
srGroups(map_full) <- factor(
  srGroups(map_full),
  c(
    "D614G",
    "B.1.1.7",
    "P.1",
    "B.1.351",
    "B.1.526+E484K",
    "B.1.617.2",
    "B.1.637",
    "C.37",
    "BA.1",
    "BA.2",
    "2x mRNA-1273",
    "3x mRNA-1273 BD01",
    "3x mRNA-1273 BD29",
    "3x mRNA-1273 (6 month)",
    "2x mRNA-1273.351"
  )
)

map_full <- orderSera(
  map_full,
  order(srGroups(map_full), srExtra(map_full))
)

# Style map
ag_colors <- c(
  "D614G"                  = "#393b79",
  "B.1.1.7"                = "#637939",
  "B.1.1.7+E484K"          = "#637939",
  "B.1.351"                = "#e7ba52",
  "B.1.617.2 (AY.1)+K417N" = "#d18652",
  "B.1.617.2 (AY.2)+K417N" = "#d18652",
  "B.1.617.2"              = "#d18652",
  "B.1.617.2 (AY.3)+E484Q" = "#d18652",
  "B.1.617.2+K417N"        = "#d18652",
  "B.1.429"                = "#9B9FD9",
  "P.1"                    = "#7b4173",
  "B.1.526+E484K"          = "#9ab370",
  "B.1.617.1"              = "#ad494a",
  "C.37"                   = "#79c9c9",
  "B.1.621"                = "#0096ad",
  "BA.1"                   = "#EF3737",
  "BA.1+A484K"             = "#EF3737",
  "BA.1.1"                 = "#EF3737",
  "BA.2"                   = "#EF3737",
  "BA.2.12.1"              = "#EF3737",
  "BA.3"                   = "#EF3737",
  "BA.4/BA.5"              = "#EF3737",
  # Mutants
  "D614G+E484K"            = "#08FDFF",
  "D614G+E484Q"            = "#0433FF",
  "D614G+L452R"            = "#AC88FF",
  "D614G+L452R+E484Q"      = "#ED82EE",
  "D614G+N501Y"            = "#26B300",
  "P.1+T417K"              = "#7b4173",
  "B.1.429+K417N"          = "#9B9FD9",
  "B.1.351+N417K"          = "#e7ba52",
  "B.1.617.1+K417N"        = "#ad494a"
)
agFill(map_full) <- ag_colors[agNames(map_full)]

sr_group_colors <- c(
  "2x mRNA-1273"           = "grey",
  "D614G"                  = "#333333",
  "B.1.1.7"                = "#637939",
  "B.1.351"                = "#e7ba52",
  "B.1.617.2"              = "#d18652",
  "B.1.617.2 (AY.12)"      = "#d18652",
  "B.1.617.2 (AY.7.1)"     = "#d18652",
  "P.1"                    = "#7b4173",
  "B.1.526+E484K"          = "#9ab370",
  "B.1.637"                = "#456C8C",
  "C.37"                   = "#79c9c9",
  "BA.1"                   = "#EF3737",
  "BA.2"                   = "#EF3737",
  "B.1.617.2 new"          = "#d18652",
  "2x mRNA-1273 new"       = "grey",
  "3x mRNA-1273"           = "grey",
  "2x mRNA-1273.351"       = colorspace::darken("#e7ba52"),
  "3x mRNA-1273 BD01"      = "grey",
  "3x mRNA-1273 BD29"      = "grey",
  "3x mRNA-1273 (6 month)" = "grey"
)
srOutline(map_full) <- sr_group_colors[as.character(srGroups(map_full))]

# Point sizes
srSize(map_full) <- 10
agSize(map_full) <- 18
agSize(map_full)[agNames(map_full) %in% c(
  agNames(map_full)[agGroups(map_full) != "wildtype"],
  "B.1.1.7+E484K",
  "B.1.617.2 (AY.1)+K417N",
  "B.1.617.2 (AY.2)+K417N",
  "B.1.617.2 (AY.3)+E484Q",
  "B.1.617.2+K417N",
  "BA.1+A484K",
  "BA.1.1",
  "BA.2",
  "BA.3",
  "BA.2.12.1",
  "BA.4/BA.5"
)] <- 12

# Add antigen and serum sequence info
base_seq <- 'MFVFLVLLPLVSSQCVNLTTRTQLPPAYTNSFTRGVYYPDKVFRSSVLHSTQDLFLPFFSNVTWFHAIHVSGTNGTKRFDNPVLPFNDGVYFASTEKSNIIRGWIFGTTLDSKTQSLLIVNNATNVVIKVCEFQFCNDPFLGVYYHKNNKSWMESEFRVYSSANNCTFEYVSQPFLMDLEGKQGNFKNLREFVFKNIDGYFKIYSKHTPINLVRDLPQGFSALEPLVDLPIGINITRFQTLLALHRSYLTPGDSSSGWTAGAAAYYVGYLQPRTFLLKYNENGTITDAVDCALDPLSETKCTLKSFTVEKGIYQTSNFRVQPTESIVRFPNITNLCPFGEVFNATRFASVYAWNRKRISNCVADYSVLYNSASFSTFKCYGVSPTKLNDLCFTNVYADSFVIRGDEVRQIAPGQTGKIADYNYKLPDDFTGCVIAWNSNNLDSKVGGNYNYLYRLFRKSNLKPFERDISTEIYQAGSTPCNGVEGFNCYFPLQSYGFQPTNGVGYQPYRVVVLSFELLHAPATVCGPKKSTNLVKNKCVNFNFNGLTGTGVLTESNKKFLPFQQFGRDIADTTDAVRDPQTLEILDITPCSFGGVSVITPGTNTSNQVAVLYQGVNCTEVPVAIHADQLTPTWRVYSTGSNVFQTRAGCLIGAEHVNNSYECDIPIGAGICASYQTQTNSPRRARSVASQSIIAYTMSLGAENSVAYSNNSIAIPTNFTISVTTEILPVSMTKTSVDCTMYICGDSTECSNLLLQYGSFCTQLNRALTGIAVEQDKNTQEVFAQVKQIYKTPPIKDFGGFNFSQILPDPSKPSKRSFIEDLLFNKVTLADAGFIKQYGDCLGDIAARDLICAQKFNGLTVLPPLLTDEMIAQYTSALLAGTITSGWTFGAGAALQIPFAMQMAYRFNGIGVTQNVLYENQKLIANQFNSAIGKIQDSLSSTASALGKLQDVVNQNAQALNTLVKQLSSNFGAISSVLNDILSRLDKVEAEVQIDRLITGRLQSLQTYVTQQLIRAAEIRASANLAATKMSECVLGQSKRVDFCGKGYHLMSFPQSAPHGVVFLHVTYVPAQEKNFTTAPAICHDGKAHFPREGVFVSNGTHWFVTQRNFYEPQIITTDNTFVSGNCDVVIGIVNNTVYDPLQPELDSFKEELDKYFKNHTSPDVDLGDISGINASVVNIQKEIDRLNEVAKNLNESLIDLQELGKYEQYIKWPWYIWLGFIAGLIAIVMVTIMLCCMTSCCSCLKGCCSCGSCCKFDEDDSEPVLKGVKLHYT'
base_seq <- strsplit(base_seq, "")[[1]]

ag_sequences <- matrix("-", numAntigens(map_full), length(base_seq), byrow = TRUE)
rownames(ag_sequences) <- agNames(map_full)
seq_info$`Pango Lineage`[seq_info$`Pango Lineage` == "B.1"] <- "D614G"
seq_info$`Pango Lineage`[seq_info$`Pango Lineage` == "B.1.617.2 (AY.3)+K417N"] <- "B.1.617.2+K417N"
seq_info$`Pango Lineage`[seq_info$`Pango Lineage` == "P.1.248"] <- "P.1"
seq_info$`Pango Lineage`[seq_info$`Pango Lineage` == "B.1.526"] <- "B.1.526+E484K"
seq_info$`Pango Lineage`[seq_info$`Pango Lineage` == "BA.1_A484K"] <- "BA.1+A484K"

for (x in seq_len(nrow(seq_info))) {
  
  variant    <- seq_info$`Pango Lineage`[x]
  aa_seq     <- base_seq
  aa_changes <- seq_info$`Mutations in Spike`[x]
  aa_changes <- trimws(strsplit(aa_changes, ",")[[1]])
  aa_changes <- gsub("in", "^", aa_changes, fixed = t)
  aa_from    <- substr(aa_changes, 1, 1)
  aa_pos     <- as.numeric(gsub("[^0-9]", "", aa_changes))
  aa_to      <- gsub(".*[0-9]", "", aa_changes)
  aa_inserts <- aa_from == "^"
  aa_orig    <- base_seq[aa_pos]
  aa_seq[aa_pos[!aa_inserts]] <- aa_to[!aa_inserts]
  aa_seq[aa_pos[aa_inserts]] <- paste0(aa_orig[aa_inserts], aa_to[aa_inserts])
  ag_sequences[variant,] <- aa_seq
  
}

# Add mutant info
ag_sequences["D614G+E484K",]       <- replace(ag_sequences["D614G",],     484, "K")
ag_sequences["D614G+E484Q",]       <- replace(ag_sequences["D614G",],     484, "Q")
ag_sequences["D614G+L452R",]       <- replace(ag_sequences["D614G",],     452, "R")
ag_sequences["D614G+L452R+E484Q",] <- replace(ag_sequences["D614G",],     c(452, 484), c("R", "Q"))
ag_sequences["D614G+N501Y",]       <- replace(ag_sequences["D614G",],     501, "Y")
ag_sequences["P.1+T417K",]         <- replace(ag_sequences["P.1",],       417, "K")
ag_sequences["B.1.429+K417N",]     <- replace(ag_sequences["B.1.429",],   417, "N")
ag_sequences["B.1.351+N417K",]     <- replace(ag_sequences["B.1.351",],   417, "K")
ag_sequences["B.1.617.1+K417N",]   <- replace(ag_sequences["B.1.617.1",], 417, "N")

# Set serum sequences
sr_sequences <- matrix(".", numSera(map_full), ncol(ag_sequences))
sr_sequences <- t(sr_sequences)
sr_sequences[,srGroups(map_full) == "2x mRNA-1273"]         <- ag_sequences["D614G",]
sr_sequences[,srGroups(map_full) == "D614G"]                <- ag_sequences["D614G",]
sr_sequences[,srGroups(map_full) == "B.1.1.7"]              <- ag_sequences["B.1.1.7",]
sr_sequences[,srGroups(map_full) == "P.1"]                  <- ag_sequences["P.1",]
sr_sequences[,srGroups(map_full) == "B.1.351"]              <- ag_sequences["B.1.351",]
sr_sequences[,srGroups(map_full) == "B.1.526+E484K"]        <- ag_sequences["B.1.526+E484K",]
sr_sequences[,srGroups(map_full) == "B.1.617.2"]            <- ag_sequences["B.1.617.2",]
sr_sequences[,srGroups(map_full) == "B.1.617.2 (AY.12)"]    <- ag_sequences["B.1.617.2",]
sr_sequences[,srGroups(map_full) == "B.1.617.2 (AY.7.1)"]   <- ag_sequences["B.1.617.2",]
# sr_sequences[,srGroups(map_full) == "B.1.637"]              <- ag_sequences["B.1.637",]
sr_sequences[,srGroups(map_full) == "C.37"]                 <- ag_sequences["C.37",]
sr_sequences[,srGroups(map_full) == "BA.1"]                 <- ag_sequences["BA.1",]
sr_sequences[,srGroups(map_full) == "BA.2"]                 <- ag_sequences["BA.2",]
sr_sequences[,srGroups(map_full) == "2x mRNA-1273 new"]     <- ag_sequences["D614G",]
sr_sequences[,srGroups(map_full) == "B.1.617.2 new"]        <- ag_sequences["B.1.617.2",]
sr_sequences[,srGroups(map_full) == "2x mRNA-1273.351"]     <- ag_sequences["B.1.351",]
sr_sequences[,srGroups(map_full) == "3x mRNA-1273"]         <- ag_sequences["D614G",]
sr_sequences <- t(sr_sequences)

# Add sequence info
agSequences(map_full) <- ag_sequences
srSequences(map_full) <- sr_sequences

# Set homologous antigens
srHomologousAgs(map_full) <- as.list(
  match(homologous_ags[as.character(srGroups(map_full))], agNames(map_full))
)

# Split maps by whether extra sera only titrated against D614G+N501Y and BA.1+A484K mutants are included
map_full_with_extras <- map_full
map_full <- subsetMap(
  map_full, 
  sera = !srExtra(map_full) %in% c(
    "COVE phase 3 (additional sera titrated against BA.1+A484K mutant)",
    "COVE phase 3 (additional sera titrated against D614G+N501Y mutant)"
  )
)

# Save the map
unlink("data/maps/", recursive = T)
dir.create("data/maps")
save.acmap(map_full, "data/maps/map_full.ace")
save.acmap(map_full_with_extras, "data/maps/map_full_with_extras.ace")


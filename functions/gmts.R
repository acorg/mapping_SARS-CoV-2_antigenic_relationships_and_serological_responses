
gmts_by_srgroup <- function(longinfo) {
  
  longinfo %>% group_by(
      sr_group,
      ag_name
    ) %>%
    summarise(
      logtiter_stat = list(
        titertools::gmt(
          titers = titer
        )
      )
    ) %>%
    mutate(
      logtiter = vapply(logtiter_stat, \(x) x["mean", "estimate"], numeric(1)),
      logtiter_upper = vapply(logtiter_stat, \(x) x["mean", "upper"], numeric(1)),
      logtiter_lower = vapply(logtiter_stat, \(x) x["mean", "lower"], numeric(1))
    )
  
}

srGroupGMTs <- function(map) {
  
  # Fetch titer table
  titer_table <- adjustedTiterTable(map)
  
  # Calculate gmts for each group and antigen
  sr_group_gmts <- lapply(
    unique(srGroups(map)),
    \(sr_group) {
      apply(
        titer_table[ , srGroups(map) == sr_group, drop = FALSE], 1, \(titers) {
          titertools::gmt(
            titers = titers,
            ci_method = "quap"
          )["mean", "estimate"]
        }
      )
    }
  )
  
  # Convert to a matrix and return the values
  sr_group_gmts <- do.call(cbind, sr_group_gmts)
  colnames(sr_group_gmts) <- as.character(unique(srGroups(map)))
  sr_group_gmts <- sr_group_gmts[ ,match(levels(srGroups(map)), colnames(sr_group_gmts))]
  sr_group_gmts
  
}


mapValues <- function(val_fn, name_fn) {
  function(map) {
    values <- val_fn(map)
    names(values) <- name_fn(map)
    values
  }
}

agFillValues <- mapValues(agFill, agNames)


# Set point drawing orders
set_pt_drawing_order <- function(map) {
  
  ag_plot_order <- c(
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
    "BA.1+A484K",
    "BA.2",
    "BA.2.12.1",
    "BA.3",
    "BA.4/BA.5",
    "BA.2.75",
    "BA.4.6",
    "BA.2.75.2",
    "BA.4_R346T",
    "BQ.1.1",
    "XBB.1",
    "D614G+E484K",
    "D614G+E484Q",
    "D614G+L452R",
    "D614G+L452R+E484Q",
    "D614G+N501Y",
    "P.1+T417K",
    "B.1.429+K417N",
    "B.1.351+N417K",
    "B.1.617.1+K417N"
  )
  
  ptDrawingOrder(map) <- c(
    seq_len(numSera(map)) + numAntigens(map),
    match(ag_plot_order[ag_plot_order %in% agNames(map)], agNames(map))
  )
  
  map
  
}


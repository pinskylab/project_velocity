# from http://stackoverflow.com/questions/9654244/multipage-lattice-panel-arrangement
# for drawing lattice plots across multiple pages

packet.panel.bycolumn <- function (layout, condlevels, page, row, column, skip) {
  dims <- sapply(condlevels, length)
  if(layout[2] != dims[2]) {
    stop("rows in layout must be equal to rows of second conditioning variable")
  }
  panels.per.row <- layout[1]
  panels.per.column <- layout[2]
  total.columns <- dims[1]
  panels.needed <- total.columns * panels.per.column
  panels.per.page <- layout[1] * layout[2]
  pages.needed <- ceiling(panels.needed / panels.per.page)
  empty.columns <- (panels.per.row - total.columns) %% panels.per.row
  panel.matrix <- rbind(matrix(1:panels.needed,ncol=panels.per.column),
                        matrix(NA, nrow=empty.columns, ncol=panels.per.column))
  panel.order <- as.vector(aperm(array(panel.matrix,
                                       dim=c(panels.per.row, pages.needed, panels.per.column)),
                                 c(1,3,2)))
  packet.order <- do.call(expand.grid, condlevels)[panel.order,]
  panel.number <- 1 + (page - 1) * panels.per.page + (row - 1) * panels.per.row + (column - 1)
  out <- as.numeric(packet.order[panel.number, ])
  if (any(is.na(out))) out <- NULL
  out
}
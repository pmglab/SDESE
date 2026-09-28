library(grid)

setwd("/public3/ly/SDESE_summary/aECS")

draw_legend <- function() {
  grid.newpage()
  txt <- function(label, x, bold = FALSE)
    grid.text(label, x, 0.51, just = "left",
              gp = gpar(fontfamily = "Arial",
                        fontsize = if (bold) 19.5 else 18.5,
                        fontface = if (bold) "bold" else "plain"))

  txt("Covariate Size (nSnpA)", 0.012, TRUE)
  grid.points(unit(0.286, "npc"), unit(0.50, "npc"), pch = 16,
              size = unit(4.0, "mm"),
              gp = gpar(col = "#E41A1C", fill = "#E41A1C"))
  txt("50", 0.312)
  grid.points(unit(0.363, "npc"), unit(0.50, "npc"), pch = 16,
              size = unit(4.0, "mm"),
              gp = gpar(col = "#377EB8", fill = "#377EB8"))
  txt("200", 0.389)
  txt("Method", 0.450, TRUE)
  grid.points(unit(0.560, "npc"), unit(0.50, "npc"), pch = 16,
              size = unit(4.0, "mm"))
  txt("Conditional aECS", 0.585)
  grid.points(unit(0.779, "npc"), unit(0.50, "npc"), pch = 17,
              size = unit(4.6, "mm"))
  txt("Conditional Z-Score", 0.801)
}

w <- 298 / 25.4
h <- 16 / 25.4
cairo_pdf("figures/legend.pdf", w, h, family = "Arial")
draw_legend()
dev.off()

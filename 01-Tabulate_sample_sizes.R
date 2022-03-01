library(stringr)
library(tidyr)
library(dplyr)

setwd("C:/Users/Quresh.Latif/files/projects/BLM/ARIM")
load("Data_compiled.RData")

rows <- c("HI", "LO", "BG", "Total")
cols <- c("Grid cells", "Cell X year", "Points", "Point X year")
out <- matrix(NA, nrow = length(rows), ncol = length(cols),
              dimnames = list(rows, cols))

strata <- ifelse(str_detect(pointXyears.list, "ARIM-HI"), "HI",
                 ifelse(str_detect(pointXyears.list, "ARIM-LO"), "LO", "BG"))
out[1:3, "Point X year"] <- tapply(strata, strata, length)[cols[1:3]]
out["Total", "Point X year"] <- length(pointXyears.list)

strata <- ifelse(str_detect(point.list, "ARIM-HI"), "HI",
                 ifelse(str_detect(point.list, "ARIM-LO"), "LO", "BG"))
out[1:3, "Points"] <- tapply(strata, strata, length)[cols[1:3]]
out["Total", "Points"] <- length(point.list)

gridXyear.list <- str_c(str_sub(pointXyears.list, 1, -9), str_sub(pointXyears.list, -5, -1)) %>% unique()
strata <- ifelse(str_detect(gridXyear.list, "ARIM-HI"), "HI",
                 ifelse(str_detect(gridXyear.list, "ARIM-LO"), "LO", "BG"))
out[1:3, "Cell X year"] <- tapply(strata, strata, length)[cols[1:3]]
out["Total", "Cell X year"] <- length(gridXyear.list)

strata <- ifelse(str_detect(grid.list, "ARIM-HI"), "HI",
                 ifelse(str_detect(grid.list, "ARIM-LO"), "LO", "BG"))
out[1:3, "Grid cells"] <- tapply(strata, strata, length)[cols[1:3]]
out["Total", "Grid cells"] <- length(grid.list)

write.csv(out, "Table_sample.csv", row.names = T)



library(data.table)
library(tidyverse)
library(limma)

clin <- fread("data/drug.naive_metadata.txt", header=TRUE) %>%
  mutate(DiseaseControl = 
           case_when(DISEASE == "HC" ~ 0,
                     DISEASE == "RA" ~ 1)
         )

edata_probes <- fread("data/gene.expression.matrix_drug.naive_raw.txt", 
                      header=TRUE) %>%
  mutate(ProbeID = sub(".*:","",ProbeID)) %>%
  dplyr::rename("SYMBOL" = "ProbeID")

edata <- avereps(edata_probes, ID = edata_probes$SYMBOL)

symbols <- edata[,1]
rownames(edata) <- symbols
edata <- edata[, -1]

edata <- edata[, clin$SAMPLE.ID]

write.csv(clin, "output/clin_drugnaive.csv", row.names=FALSE)
write.csv(edata, "output/edata_drugnaive.csv", row.names=TRUE)

# normalised data
edata_probes_norm <- fread("data/gene.expression.matrix_drug.naive_age.rin.removed.txt", 
                      header=TRUE) %>%
  mutate(ProbeID = sub(".*:","",ProbeID)) %>%
  dplyr::rename("SYMBOL" = "ProbeID")

edata_norm <- avereps(edata_probes_norm, ID = edata_probes_norm$SYMBOL)

symbols_norm <- edata_norm[,1]
rownames(edata_norm) <- symbols_norm
edata_norm <- edata_norm[, -1]

edata_norm <- edata_norm[, clin$SAMPLE.ID]

write.csv(clin, "output/clin_drugnaive_norm.csv", row.names=FALSE)
write.csv(edata_norm, "output/edata_drugnaive_norm.csv", row.names=TRUE)



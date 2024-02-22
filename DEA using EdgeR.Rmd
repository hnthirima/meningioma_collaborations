---
title: "Differential analysis between tumors in cluster A and cluster B"
output: html_notebook
author: Nayanga Thirimanne
---
Load packages
```{r}
library(tidyverse)
library(edgeR)
library(ggplot2)
library(ggrepel)
```

Load data
```{r}
#batch corrected raw counts
rawcounts1 <- readRDS("~/1298_combatseq_rawcounts.rds")

clusterID <- read_rds("~/dbscan_1298_cl9_uncl48.rds")
clusterID$cluster[clusterID$cluster == "1"] <- "A"
clusterID$cluster[clusterID$cluster == "2"] <- "B"
clusterID$cluster[clusterID$cluster == "3"] <- "C"
clusterID$cluster[clusterID$cluster == "4"] <- "D"
clusterID$cluster[clusterID$cluster == "5"] <- "E"
clusterID$cluster[clusterID$cluster == "6"] <- "F"
clusterID$cluster[clusterID$cluster == "7"] <- "G"
clusterID$cluster[clusterID$cluster == "8"] <- "H"
clusterID$cluster[clusterID$cluster == "9"] <- "I"

```

Generate sample dataframe
```{r}
#select samples in cluster A and B
clAandB <- clusterID %>%
            filter(clusterID$cluster == "A" | clusterID$cluster == "B")

rawcounts1 <- rawcounts1[, (colnames(rawcounts1) %in% clAandB$coordinate_ID)]

sample <- colnames(rawcounts1)
Adf <- data.frame(sample)

clA <- clAandB %>% filter(clAandB$cluster == "A")
clB <- clAandB %>% filter(clAandB$cluster == "B")

Adf$group[Adf$sample %in% clA$coordinate_ID] <- 'A'
Adf$group[Adf$sample %in% clB$coordinate_ID] <- 'B'

View(Adf)
```

EdgeR analysis
```{r}

#EdgeR DE analysis: GLM approach (quasi-likelihood F-tests)
#set group
group = factor(Adf$group, levels = c("A","B"))
cols_keep <- Adf$sample
rawcounts <- rawcounts1[, colnames(rawcounts1) %in% cols_keep]
y <- DGEList(counts = rawcounts, genes = rownames(rawcounts), samples = colnames(rawcounts), group = group)
#filter out the lowly expressed genes
#keep <- filterByExpr(y)
#y <- y[keep, , keep.lib.sizes=FALSE]
#recompute libsize after filtering
y$samples$lib.size <- colSums(y$counts)
y <- calcNormFactors(y)
#y$samples
#set design
design <- model.matrix(~0+group)
rownames(design) <- colnames(y)
#design
#estimate the NB disperson for the dataset
y <- estimateDisp(y, design, robust = TRUE)
y$common.dispersion
y$tagwise.dispersion
plotBCV(y) 
fit <- glmQLFit(y, design)
my.contrasts1 <- makeContrasts(AvsB=groupA-groupB, levels=design)
qlf <- glmQLFTest(fit, contrast = my.contrasts1) 
topTags(qlf)
W <- as.data.frame(topTags(qlf), n=nrow(y))
out <- topTags(qlf, n=Inf)
AvsB <- as.data.frame(out)
write.csv(AvsB, "~/DE_1298_clusterAvsB_unfilt_022024.csv")

```

Volcano plot
```{r}
DEoutput <- read.csv("~/DE_1298_clusterAvsB_unfilt_022024.csv")
DEoutput2 <- DEoutput %>%
  mutate(threshold = factor(case_when(logFC > 0.6 & FDR  < 0.05 ~ "Upregulated in cluster A",
                                      logFC < -0.6 & FDR < 0.05 ~ "Downregulated in cluster A",
                                      TRUE ~ "Non-diff.")))

volc <- ggplot(DEoutput2, aes(x = logFC, y=-log10(FDR))) + 
  geom_point(aes(color=DEoutput2$threshold), size = 0.8) +
  scale_color_manual(name = " ", values = c("Upregulated in cluster A" = "red", "Downregulated in cluster A"="blue", "Non-diff."="grey"))  +
  xlab("logFC") +
  ylab("-log10(FDR)") +
  theme_bw() +
  theme(legend.text = element_text(size = 20), text = element_text(size = 20),
        legend.position = "bottom",) + 
  xlim(-6,6) +
  labs(title = "DE analysis between cluster A and B")

DE_genes <- DEoutput2 %>%
  filter(DEoutput2$genes %in% c("CCN1", "CCN2", "AMOT", "AMOTL2",
                                "ANKRD1", "CITED2", "VGLL4", "NF2", "TEAD3", "TEAD4",
                                "YAP1", "IGF2", "IGF2BP1", "IGF2BP2", "IGFBP1", "IGFBP3"))

y_limits <- c(80, NA)
p <- volc + geom_text_repel(data = DE_genes, aes(label = genes),
                         box.padding = 4, max.overlaps = Inf, ylim = y_limits, force = 5) 
p
setwd("~/DEA/")
ggsave("VolcanoPlot_DE_clusterAvsB.pdf", height = 10, width = 10)
```


data <- read.csv("~/Downloads/SignificantGenes.csv", header=TRUE, stringsAsFactors=FALSE)
avg_il <- data$avera_il
avg_lps <- data$avera_lps

mtx <- matrix(1:150, nrow=3, ncol=50)
rownames(mtx) <- c("il4", "lps", "ratio")
colnames(mtx) <- data$external_gene_name[1:50]
mtx["il4",] <- avg_il[1:50]
mtx["lps",] <- avg_lps[1:50]

for(i in 1:50){
  mtx["ratio",i] <- avg_lps[i] / avg_il[i]
}



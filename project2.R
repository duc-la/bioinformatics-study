#Question 1
library("GEOquery")

eset<-getGEO("GSE19804", filename="GSE19804_series_matrix.txt.gz")

data <- exprs(eset)

pdata <- pData(eset)

#Question 2
control<-rownames(pdata[grep("Lung Normal",pdata$title),])

cancer<-rownames(pdata[grep("Lung Cancer",pdata$title),])

cancer.indices <- which(colnames(data) %in% cancer)
control.indices <- which(colnames(data) %in% control)

#T-score
TScore<-apply(data, 1, FUN = function(x, cancer, control) {
  t.test(x[cancer], x[control])$statistic
}, cancer = cancer, control = control)

#P-Value
PValue<-apply(data, 1, FUN = function(x, cancer, control) {
  t.test(x[cancer], x[control])$p.value
}, cancer = cancer, control = control)

#logFC
diff <- function(x, cancer, control) {
  mean(x[cancer]) - mean(x[control])
}

logFC <- apply(data, 1, FUN = diff, cancer = cancer, control = control)

geneID = rownames(data)

#Data frame
GSEData <- data.frame(geneID, PValue, TScore, logFC)

saveRDS(GSEData, file = "DE-results.rds")

#Question 3
library("ggplot2")
pdf("volcano.pdf")
p <- ggplot(data = GSEData, aes(x=logFC, y = -log10(PValue))) + geom_point(col = ifelse((abs(logFC)>1 & PValue<0.05), "red", "black"))
print(p)
dev.off()

#Question 4
t_OBSERVED <- data.frame(geneID, TScore)
saveRDS(t_OBSERVED, file = "t-observed.rds")

#Question 5
t_NULL_DISTRIBUTION <- NULL
for (i in 1:100){
  random_index <- sample(1:ncol(data))
  data_permuted <- data[, random_index]
  
  t_NULL <- apply(data_permuted, 1, FUN = function(x, cancer, control){
    t.test(x[cancer], x[control])$statistic
  }, cancer = cancer.indices, control=control.indices)
  
  t_NULL_DISTRIBUTION <- cbind(t_NULL_DISTRIBUTION, t_NULL)
}

pT <- NULL
for (i in 1:nrow(data)){
  t_OBSERVED_i <- t_OBSERVED$TScore[i]
  t_NULL_DIST_i <- abs(t_NULL_DISTRIBUTION[i, ]) > abs(t_OBSERVED_i)
  pT[i] <- sum(t_NULL_DIST_i)/100
}

saveRDS(pT, "p-empirical-t-score.rds")

#Question 6
e_Score<-apply(data, 1, FUN = function(x, cancer, control) {
  abs(mean(x[cancer]) - mean(x[control]))
}, cancer = cancer, control = control)

geneID = rownames(data)
e_OBSERVED <- data.frame(geneID, e_Score)

e_NULL_DISTRIBUTION <- NULL
for (i in 1:100){
  random_index <- sample(1:ncol(data))
  data_permuted <- data[, random_index]
  
  e_NULL <- apply(data_permuted, 1, FUN = function(x, cancer, control){
    abs(mean(x[cancer]) - mean(x[control]))
  }, cancer = cancer.indices, control=control.indices)
  
  e_NULL_DISTRIBUTION <- cbind(e_NULL_DISTRIBUTION, e_NULL)
}

pE <- NULL
for (i in 1:nrow(data)){
  e_OBSERVED_i <- e_OBSERVED$e_Score[i]
  e_NULL_DIST_i <- abs(e_NULL_DISTRIBUTION[i, ]) > abs(e_OBSERVED_i)
  pE[i] <- sum(e_NULL_DIST_i)/100
}

saveRDS(pE, "p-empirical-euclidean.rds")

#Question 7
pdf("hist-pT.pdf")
hist(pT,
     main= "Histogram of pT",
     xlab= "pT")
dev.off()

pdf("hist-pE.pdf")
hist(pE,
     main= "Histogram of pE",
     xlab= "pE")
dev.off()

#Question 8
cor(pT, pE, method = "pearson")

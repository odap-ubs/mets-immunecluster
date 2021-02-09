---
title: "Marker gene"
author: "García-Mulero S"
date: "3/30/2020"
output: html_document
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

## Identify gene to classify samples into ImmuneClusters
Load data 
```{r}
infoClin <- read.delim("data/mets_infoClin.txt", sep="\t", header=T, stringsAsFactors = F)
ex <- read.table("data/mets_expression_combat.txt", sep="\t", header=T, stringsAsFactors = F)
```

Generate training and testing data
```{r}
labels <- infoClin
labelsAll <- labels
# filter only high and low samples
labels <- labels[!labels$cluster=="Medium",]
head(labels); dim(labels)
table(labels$cluster)
# Take training data - 75 % of samples
head(labels)
high <- as.character(labels$GEO_ID[labels$cluster=="High"])
length(high)*0.75 # take 55 samples
set.seed(123)
high_train <- sample(high, 52, replace=F) # take high
low <- as.character(labels$GEO_ID[labels$cluster=="Low"])
length(low)*0.75 # take 85 samples
set.seed(123)
low_train <- sample(low, 85, replace=F) # take low
labels_train <- labels[labels$GEO_ID %in% c(high_train, low_train),]
head(labels_train); dim(labels_train) # 137
table(labels_train$cluster);
# Take the rest of samples - 25 %
labels_test <- labels[! (labels$GEO_ID %in% c(high_train, low_train)),]
head(labels_test); dim(labels_test) # 46
table(labels_test$cluster) 
```

Differential expression High vs Low
```{r}
### Limma High vs Low
library(limma)
# filter only training
ex2 <- ex[, colnames(ex) %in% labels_train$GEO_ID]
dim(ex2)
ex2 <- ex2[, order(colnames(ex2))]
infoClin2 <- infoClin[infoClin$GEO_ID %in% colnames(ex2), ]
ex2[1:5,1:5]
infoClin2[1:5,1:5]
dim(ex2); dim(infoClin2) # 137
ex2 <- ex2[, order(match(colnames(ex2), infoClin2$GEO_ID))]
identical(as.character(infoClin2$GEO_ID), as.character(colnames(ex2)))
# create design matrix
ImmunoClust <- as.factor(infoClin2$cluster)
design <- model.matrix(~0+ImmunoClust)
colnames(design) <- levels(ImmunoClust)
design
# groups comparison
fit <- lmFit(ex2, design)
names(fit)
# Extract the linear model fit for the contrasts
contrast.matrix <- makeContrasts(High-Low, levels=design)
contrast.matrix
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
# Take up/downregulated genes in two lists 
results <- decideTests(fit2, method = "separate",	adjust.method = "bonferroni", p.value = 0.01, lfc=2)
summary(results)
aux <- results@.Data
de <- as.data.frame(aux)
de$genes <- rownames(de)
DEGenes_High <- de$genes[de$`High - Low`==1]
DEGenes_Low <- de$genes[de$`High - Low`==-1]
# select genes for linear model
DiffExpr <- topTable(fit2, n=Inf, p.value = 0.01, lfc=2, sort.by="logFC", adjust.method = "bonferroni") 
DiffExpr <- DiffExpr[order(abs(DiffExpr$adj.P.Val), decreasing=FALSE), ]
dim(DiffExpr) # 43
write.table(DiffExpr, "data/DEG_genes_list_filtered.txt", sep="\t", quote=F, col.names=T, row.names=T)
```
Plot genes 
```{r}
top.genes <- rownames(DiffExpr)[1:10]
labels_genes <- merge(infoClin, t(ex[top.genes,]), by.x="GEO_ID", by.y="row.names")
labels_genes <- melt(labels_genes)
labels_genes$cluster <- factor(labels_genes$cluster, levels=c("Low", "Medium", "High"))
data <- labels_genes
### Function
plotGene <- function(gene, data){
              hline <- as.numeric(c(quantile(data$value[data$variable==gene], 0.25),
                                median(data$value[data$variable==gene]),
                                quantile(data$value[data$variable==gene], 0.75)))
              graph <- ggplot(data[data$variable==gene,], aes(x=cluster, y=value)) +
                      	geom_boxplot(color="black", outlier.shape = NA) +
                      	geom_hline(yintercept = hline, linetype="dashed") +
                        geom_jitter(position=position_jitter(0.1), size=0.7, 
                              aes(x=cluster, y=value, color = cluster)) +
                        labs(title=gene, x="", y = "log2(expression)") +
                          scale_color_manual(guide=FALSE, values=c("grey80", "grey60", "grey40")) +
                          scale_y_continuous(limits = c(4,12)) +
                          stat_compare_means(method="kruskal.test", label="p.format", size=3) +
		                theme_classic() +
						theme(legend.position="none",
								plot.title = element_text(hjust = 0.5), panel.spacing = unit(1, "lines"),
								axis.text.x = element_text(size=10),
					  			axis.text.y = element_text(size=8),
					  			strip.text.x = element_text(size = 10), 
					  			panel.grid.major = element_blank(), 
					  			panel.grid.minor = element_blank(),
								panel.background = element_rect(colour = "black", size=0.5),
								strip.background = element_rect(colour = "white", size=0.5))
              return(graph)
}
for (i in 1:length(top.genes)){
  gene <- top.genes[i]
  graph <- plotGene(gene, labels_genes)
  eval(parse(text=paste0("graph", i, "<- graph")))
}
graph1
require(gridExtra)
png("~/Desktop/eris-tmp/revision_mets/def2/Supp_Figure_9.tiff", 
    units="in", width=12, height=5, res=600)
grid.arrange(graph1, graph2, graph3, graph4, graph5, 
			 graph6, graph7, graph8, graph9, graph10, 
			 nrow = 2, ncol=5)
dev.off()
```


Parse data for model fitting
```{r}
## add info genes
genes <- rownames(DiffExpr)
# add info genes to labelsAll (the expression level is with all data)
exGenes <- ex[genes,]
rownames(exGenes) <- gsub( "\\-", "_", rownames(exGenes))
exGenes <- as.data.frame(t(exGenes))
dim(exGenes)
# merge
rownames(labelsAll) <- labelsAll$GEO_ID
labelsAll <- merge(labelsAll, exGenes , by="row.names")
# generate dichotomous
for (i in 14:ncol(labelsAll)){
   labelsAll[,i] <- ifelse(labelsAll[,i] > median(labelsAll[,i]), "High_ex", "Low_ex")
   labelsAll[,i] <- factor(labelsAll[,i], levels=c("Low_ex", "High_ex"))
}
# create data frame only with info for the model
labels_train2 <- labelsAll[labelsAll$GEO_ID %in% labels_train$GEO_ID, ]
train <- labels_train2[,c(9, 12:ncol(labels_train2))]
train$cluster <- factor(train$cluster)
table(train$cluster)
```
Decision tree model
```{r}
require(rpart); library(caret); library(rpart.plot)
## fit the model
set.seed(1)
caret.control = trainControl(method = "cv", 
	number = 5, repeats = 3, search="random")
fit = train(cluster ~ ., 
	data = train, 
	method = "rpart", 
	trControl = caret.control)
fit$finalModel 
```
Validation of the model (test data)
```{r}
## prepare testing data
labels_test2 <- labelsAll[labelsAll$GEO_ID %in% labels_test$GEO_ID, ]
## Validate model
pred <- predict(fit, newdata=labels_test2, method="prob")
names(pred) <- labels_test$GEO_ID
labels_test2$Predicted <- pred
# Add to labels
with(labels_test2, table(cluster, Predicted))
conf_matrix <- as.data.frame(with(labels_test2, table(cluster, Predicted)))
conf_matrix
# sensitivity
sens = (conf_matrix$Freq[1])/(conf_matrix$Freq[1]+conf_matrix$Freq[3])
#specificity
specif = (conf_matrix$Freq[4])/(conf_matrix$Freq[4]+conf_matrix$Freq[2])
omspec = 1 - specif
# calculate error
error = (conf_matrix$Freq[2]+conf_matrix$Freq[3])/sum(conf_matrix$Freq)
# accuracy
accuracy = (sens + specif) / sum(conf_matrix$Freq)
# AUC
library(pROC)
labels_test2$pred2 <- ifelse(labels_test2$Predicted=="High", 1, 0)
res.roc <- roc(labels_test2$cluster, labels_test2$pred2)
plot.roc(res.roc, print.auc = TRUE, main="")
# results
sens; 
specif; 
error; 
accuracy; 
res.roc$auc
```
Validation on dataset GSE51244
```{r}
infoClin <- read.table("data/clin_gse51244.1.txt", sep="\t", header=T, stringsAsFactors=F)
ex <- read.table("data/ex_gse51244.1.txt", sep="\t", header=T, stringsAsFactors=F)
load("data/pathways_for_gsva_mets.Rdata")

# Filter metastases
infoClin$MET_SITE <- gsub(" ", "", infoClin$MET_SITE)
infoClin <- infoClin[infoClin$TISSUE=="Metastasis",]
infoClin <- infoClin[grep("Liver|Lung", infoClin$MET_SITE),]
infoClin <- infoClin[infoClin$GEO_ID!="GSM1019407",]

# Look for the biomarker
ex <- ex[, colnames(ex) %in% infoClin$GEO_ID]

    ### Perform GSVA function ---------------------------------------
    library(GSVA)
    GSVA <- gsva(expr=as.matrix(ex), gset.idx.list=Pathways, min.sz=1, verbose=F)
    
    data <- t(GSVA)
    d <- dist(data, method = "euclidean") # distance matrix
    hcl1 <- hclust(d, method="ward.D2") 

    ### Create three groups of samples  
    ImmunoCluster <- cutree(hcl1, k=3)
    ImmunoCluster <- ImmunoCluster[order(factor(names(ImmunoCluster),
      levels=infoClin$GEO_ID))]
    identical(as.character(names(ImmunoCluster)), as.character(infoClin$GEO_ID))

    boxplot(as.numeric(ex["CD74",colnames(ex) %in% names(ImmunoCluster[ImmunoCluster==1])]),
    as.numeric(ex["CD74",colnames(ex) %in% names(ImmunoCluster[ImmunoCluster==2])]),
    as.numeric(ex["CD74",colnames(ex) %in% names(ImmunoCluster[ImmunoCluster==3])]))

    # add to infoClin
    infoClin$ImmunoCluster <- ImmunoCluster
    infoClin$cluster <- NA
    infoClin$cluster[infoClin$ImmunoCluster==1] = "Low"
    infoClin$cluster[infoClin$ImmunoCluster==2] = "Medium"
    infoClin$cluster[infoClin$ImmunoCluster==3] = "High"
    infoClin$cluster <- factor(infoClin$cluster, levels=c("Low", "Medium", "High"))


### ------------------------------------------------------------------------------------------
### CD74
### ------------------------------------------------------------------------------------------

# Create a function to validate in different datasets
require(reshape)
validateGene <- function(ex, Pathways, infoClin, gene){

  ### Look at prediction by gene expression 
  # add gene expression
  labels <- cbind(infoClin, as.data.frame(t(ex[gene,])))
  labels$geneEx <- NA
  labels$geneEx <- ifelse(labels[, gene] <= median(labels[,gene]), "Low" , "High")
  # remove medium 
  # labels <- labels[labels$cluster!="Medium",]
  labels$geneEx <- factor(labels$geneEx, levels=c("Low", "Medium", "High"))

  # Estimate error ---------------------------------------
  print(with(labels, table(geneEx, cluster)))
  print(with(labels, table(geneEx, MET_SITE)))

  conf_matrix <- as.data.frame(with(labels, table(geneEx, cluster)))
  error = (conf_matrix$Freq[2]+conf_matrix$Freq[5])/sum(conf_matrix$Freq)
  print(error)

  # generate data for plot
  labels$esp[labels$cluster=="High" & labels$geneEx=="Low"] = "False Negative"
  labels$esp[labels$cluster=="Low" & labels$geneEx=="High"] = "False Positive" 
  labels$esp[is.na(labels$esp)] = "NA"

  data <- melt(labels[,c("DATASET", "GEO_ID", gene, "geneEx", "cluster", "esp")])
  data$cluster <- factor(data$cluster, levels=c("Low", "Medium", "High"))
  return(data)

}


datCD74 <- validateGene(ex, Pathways, infoClin, "CD74")
conf_matrix <- as.data.frame(with(datCD74, table(cluster, geneEx)))

# sensitivity
sens = (conf_matrix$Freq[9])/(conf_matrix$Freq[9]+conf_matrix$Freq[7])

#specificity
specif = (conf_matrix$Freq[1])/(conf_matrix$Freq[1]+conf_matrix$Freq[2])

# calculate error
error = (conf_matrix$Freq[2]+conf_matrix$Freq[3]+conf_matrix$Freq[7]+conf_matrix$Freq[8])/sum(conf_matrix$Freq)

# accuracy
accuracy = (sens + specif) / sum(conf_matrix$Freq)

# AUC
res.roc <- roc(labels_test$cluster, labels_test$pred_prob)
plot.roc(res.roc, print.auc = TRUE, main="")

# results
sens; 
specif; 
error; 
accuracy; 
res.roc$auc


```


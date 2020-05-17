load("signatures.Rda")

library(GEOquery)
library(scatterplot3d)
library(genefilter)

strong_sig <- read.table("strong_signature.txt")[,1]
alpha <- 0.05

# rand_sig_0 <- as.list(c()) #contain no SPS
# rand_sig_5 <- as.list(c())
# rand_sig_10 <- as.list(c())
# rand_sig_20 <- as.list(c()) #contain all SPS
cancer.signatures <- as.list(c())

#non_SPS <- setdiff(as.vector(entrez_IDs[,1]),strong_sig)

for (i in 1:1000)
{
  #cancer.signatures[[i]] <- noquote(as.vector(sample(entrez_IDs[,1],83)))
  cancer.signatures[[i]] <- noquote(as.vector(sample(entrez_IDs[,1],20)))
  #just take from 5027 below 
  #rand_sig_0[[i]] <- sample(non_SPS,20)
  #rand_sig_5[[i]] <- c(sample(non_SPS,15), sample(strong_sig,5))
  #rand_sig_10[[i]] <- c(sample(non_SPS,10), sample(strong_sig,10))
  #rand_sig_20[[i]] <- sample(strong_sig,20)
}

sps_overlaps <- c()

for (i in 1:length(cancer.signatures))
{
  #sig_names <- append(sig_names,cancer.signatures[[i]]$name)
  sps_overlaps <- append(sps_overlaps, length(intersect(cancer.signatures[[i]], strong_sig)))
}  

table(sps_overlaps)

#GSE22544
GSE22544<- getGEO("GSE22544", GSEMatrix=T,AnnotGPL=T,getGPL=T)
disease_fact <- pData(phenoData(GSE22544[[1]]))[,11]
table(disease_fact)
GSE22544_eset <- GSE22544[[1]]
colnames(fData(GSE22544_eset))[4]
entrez_IDs_22544 <-(fData(GSE22544_eset)[4])

GSE22544_eset_exprs <- exprs(GSE22544_eset)

GSE22544_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_22544[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GSE22544_outcome <- append(GSE22544_outcome,0) 
  }
  else
  {
    data_mat <- GSE22544_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GSE22544_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GSE22544_class_pval <- append(GSE22544_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GSE22544_class_pval <- append(GSE22544_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GSE22544_outcome <- append(GSE22544_outcome, sum(as.numeric(GSE22544_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GSE22544_outcome <- append(GSE22544_outcome, sum(c(6,4,2)[GSE22544_class_pval[1:3] <= 0.05]))
    GSE22544_outcome <- append(GSE22544_outcome, min(GSE22544_class_pval))
  }
}


#GDS4061
GDS4061 <- getGEO("GDS4061", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS4061@dataTable@columns[["genotype/variation"]]
table(GDS4061@dataTable@columns[["genotype/variation"]])
GDS4061_eset <- GDS2eSet(GDS4061,do.log2=F)

entrez_IDs_4061 <-(fData(GDS4061_eset)[4])

GDS4061_eset_exprs <- exprs(GDS4061_eset)


GDS4061_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_4061[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS4061_outcome <- append(GDS4061_outcome,0) 
  }
  else
  {
    data_mat <- GDS4061_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS4061_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS4061_class_pval <- append(GDS4061_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS4061_class_pval <- append(GDS4061_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS4061_outcome <- append(GDS4061_outcome, sum(as.numeric(GDS4061_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS4061_outcome <- append(GDS4061_outcome, sum(c(6,4,2)[GDS4061_class_pval[1:3] <= 0.05]))
    GDS4061_outcome <- append(GDS4061_outcome, min(GDS4061_class_pval))
  }
}

#GDS4095
GDS4095 <- getGEO("GDS4095", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS4095@dataTable@columns[["genotype/variation"]]
table(GDS4095@dataTable@columns[["genotype/variation"]])
GDS4095_eset <- GDS2eSet(GDS4095,do.log2=F)
colnames(fData(GDS4095_eset))[4]
entrez_IDs_4095 <-(fData(GDS4095_eset)[4])

GDS4095_eset_exprs <- exprs(GDS4095_eset)


GDS4095_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_4095[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS4095_outcome <- append(GDS4095_outcome,0) 
  }
  else
  {
    data_mat <- GDS4095_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS4095_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS4095_class_pval <- append(GDS4095_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS4095_class_pval <- append(GDS4095_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS4095_outcome <- append(GDS4065_outcome, sum(as.numeric(GDS4065_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS4095_outcome <- append(GDS4065_outcome, sum(c(6,4,2)[GDS4065_class_pval[1:3] <= 0.05]))
    GDS4095_outcome <- append(GDS4095_outcome, min(GDS4095_class_pval))
  }
}

#GDS4092
GDS4092 <- getGEO("GDS4092", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS4092@dataTable@columns[["cell.line"]]
table(GDS4092@dataTable@columns[["cell.line"]])
GDS4092_eset <- GDS2eSet(GDS4092,do.log2=F)
colnames(fData(GDS4092_eset))[4]
entrez_IDs_4092 <-(fData(GDS4092_eset)[4])

GDS4092_eset_exprs <- exprs(GDS4092_eset)


GDS4092_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_4092[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS4092_outcome <- append(GDS4092_outcome,0) 
  }
  else
  {
    data_mat <- GDS4092_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS4092_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS4092_class_pval <- append(GDS4092_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS4092_class_pval <- append(GDS4092_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS4095_outcome <- append(GDS4065_outcome, sum(as.numeric(GDS4065_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS4095_outcome <- append(GDS4065_outcome, sum(c(6,4,2)[GDS4065_class_pval[1:3] <= 0.05]))
    GDS4092_outcome <- append(GDS4092_outcome, min(GDS4092_class_pval))
  }
}

#GDS2250
GDS2250 <- getGEO("GDS2250", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS2250@dataTable@columns[["disease.state"]]
table(GDS2250@dataTable@columns[["disease.state"]])
GDS2250_eset <- GDS2eSet(GDS2250,do.log2=F)
colnames(fData(GDS2250_eset))[4]
entrez_IDs_2250 <-(fData(GDS2250_eset)[4])

GDS2250_eset_exprs <- exprs(GDS2250_eset)


GDS2250_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_2250[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS2250_outcome <- append(GDS2250_outcome,0) 
  }
  else
  {
    data_mat <- GDS2250_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS2250_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS2250_class_pval <- append(GDS2250_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS2250_class_pval <- append(GDS2250_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(as.numeric(GDS2250_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(c(6,4,2)[GDS2250_class_pval[1:3] <= 0.05]))
    GDS2250_outcome <- append(GDS2250_outcome, min(GDS2250_class_pval))
  }
}

#GDS4818
GDS4818 <- getGEO("GDS4818", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS4818@dataTable@columns[["genotype/variation"]]
table(GDS4818@dataTable@columns[["genotype/variation"]])
GDS4818_eset <- GDS2eSet(GDS4818,do.log2=F)
colnames(fData(GDS4818_eset))[4]
entrez_IDs_4818 <-(fData(GDS4818_eset)[4])

GDS4818_eset_exprs <- exprs(GDS4818_eset)


GDS4818_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_4818[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS4818_outcome <- append(GDS4818_outcome,0) 
  }
  else
  {
    data_mat <- GDS4818_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS4818_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS4818_class_pval <- append(GDS4818_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS4818_class_pval <- append(GDS4818_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(as.numeric(GDS2250_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(c(6,4,2)[GDS2250_class_pval[1:3] <= 0.05]))
    GDS4818_outcome <- append(GDS4818_outcome, min(GDS4818_class_pval))
  }
}

#GDS2760
GDS2760 <- getGEO("GDS2760", GSEMatrix=TRUE,AnnotGPL=FALSE,getGPL=TRUE)
disease_fact <- GDS2760@dataTable@columns[["protocol"]]
table(GDS2760@dataTable@columns[["protocol"]])
GDS2760_eset <- GDS2eSet(GDS2760,do.log2=F)
colnames(fData(GDS2760_eset))[4]
entrez_IDs_2760 <-(fData(GDS2760_eset)[4])

GDS2760_eset_exprs <- exprs(GDS2760_eset)


GDS2760_outcome <- c()
for (i in 1:length(cancer.signatures))
{
  print(i)
  positions <- as.numeric(entrez_IDs_2760[,1] %in% cancer.signatures[[i]])
  
  if (sum(positions) == 0)
  {
    GDS2760_outcome <- append(GDS2760_outcome,0) 
  }
  else
  {
    data_mat <- GDS2760_eset_exprs
    data_mat <- data_mat[which(positions==1),]
    pca_norm <- prcomp(t(data_mat),scale=T, center=T) 
    
    GDS2760_class_pval <- c()
    
    if (ncol(pca_norm$x) < 10)
    {
      for (j in 1:ncol(pca_norm$x))
      {
        GDS2760_class_pval <- append(GDS2760_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }
      
    }
    else
    {
      for (j in 1:10)
      {
        GDS2760_class_pval <- append(GDS2760_class_pval,kruskal.test(pca_norm$x[,j]~disease_fact)$p.value)
      }  
    }  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(as.numeric(GDS2250_class_pval[1:3] <= 0.05)))  #if at least 1 of the 3 PCs is significant  
    #GDS2250_outcome <- append(GDS2250_outcome, sum(c(6,4,2)[GDS2250_class_pval[1:3] <= 0.05]))
    GDS2760_outcome <- append(GDS2760_outcome, min(GDS2760_class_pval))
  }
}

outcome_test <- cbind(paste('R_',1:1000, sep=""), sps_overlaps, GSE22544_outcome, GDS4061_outcome, GDS4095_outcome, GDS4092_outcome, GDS2250_outcome, GDS4818_outcome, GDS2760_outcome)
#outcome_test <- cbind(paste('R_',1:1000, sep=""), sps_overlaps, GSE22544_outcome, GDS4061_outcome, GDS4095_outcome, GDS4092_outcome, GDS4818_outcome, GDS2760_outcome)

write.table(outcome_test, file="random_lowestPCpval_48_sig_7GDS_test.txt", quote=F, row.names=F, col.names=T, sep="\t")
#write.table(outcome_test, file="random_lowestPCpval_48_sig_6GDS_test.txt", quote=F, row.names=F, col.names=T, sep="\t")


results_test <- read.table("random_lowestPCpval_48_sig_7GDS_test.txt", header=T,sep="\t")
min_SPS_test <- read.table("SPS_PC1to10_lowest_pval_test.txt", sep="\t", header=F)[,2]


across_dataset_count_test <- c()
for (i in 1:nrow(results_test))
{
  across_dataset_count_test <- append(across_dataset_count_test, sum(results_test[i,3:9] <= min_SPS))
  #across_dataset_count <- append(across_dataset_count, sum(results[i,3:5] <= min_SPS))
  #across_dataset_count <- append(across_dataset_count, sum(results[i,3:4] <= min_SPS))
}

table(across_dataset_count_test < 1)
#hist(across_dataset_count)
expected_dist <- rbinom(1000, 7, 0.46)
expected_dist <- rbinom(1000, 6, 0.46)
#expected_dist <- rbinom(1000, 3, 0.46)

#the 48 signatures
known_sig <- read.table("lowest_pval_PC1to10_48_sig_7GDS.txt", header=T,sep="\t")


across_known_sig <- c()
for (i in 1:nrow(known_sig))
{
  across_known_sig <- append(across_known_sig, sum(known_sig[i,3:9] <= 0.05))
  #across_known_sig <- append(across_known_sig, sum(known_sig[i,c(4:5,8)] <= 0.05))
}  


pdf("expected_distribution3_test.pdf")
#hist(across_dataset_count, xlim=c(0,7), main=NULL, xlab="Significance in x datasets", prob=T, col="blue")
hist(expected_dist, col="red", ylim=c(0,1000), xlim=c(0,7), main=NULL)
#hist(expected_dist, col="red", ylim=c(0,1000), xlim=c(0,3), main=NULL)
hist(across_dataset_count_test, col="blue", add=T)
hist(across_known_sig, col="yellow", add=T, lty=3)
box()
# abline(v=7, col="pink", lty=2)
dev.off()

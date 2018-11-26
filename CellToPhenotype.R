## Code for CellToPhenotype predictor
## Authors: Nishanth Ulhas Nair and Avinash Das

library(glmnet)

# normalization function 
qnorm.array <- function(mat)
{
  mat.back = mat 
  mat = mat[!is.na(mat)]
  mat = rank(mat,  rank, ties.method = "average");
  mat = qnorm(mat / (length(mat)+1));
  mat = mat/sd(mat)
  mat.back[!is.na(mat)] = mat 
  mat.back
}


load("store_data_for_celltophenotype.RData") ## saved data required for running the code
## actual_mig -- experimentally measured migration values in 43 cell lines
## actual_prol -- experimentally measured migration values in 46 cell lines
## tcga_genes_uppercase -- gene names in the TCGA data
## cellline_genes_uppercase -- gene names in the cell line data
## expression_for_migration -- cell line normalized expression data for 43 cell lines used for migration model
## expression_for_proliferation -- cell line normalized expression data for 46 cell lines used for migration model
## match.gene.tcga -- 2448 genes selected from survival analysis from METABRIC dataset
## mRNA.norm -- normalized gene expression data for 1043 breast cancer TCGA patients


pheno.inx.tcga = match(match.gene.tcga, cellline_genes_uppercase)
tcga.inx = match(match.gene.tcga, tcga_genes_uppercase)

no_of_loops = 50 # no. of iterations

tcga.motility_store = matrix(0,no_of_loops,dim(mRNA.norm)[2])
tcga.GR.high_store = matrix(0,no_of_loops,dim(mRNA.norm)[2])
for (i in 1:no_of_loops)
{
  print(i)
  cv.fit.motility.tcga = cv.glmnet(x=t(expression_for_migration[pheno.inx.tcga,]), y=actual_mig, alpha=1,  family="gaussian",nfold=5, standardize = F) # LASSO regression for migration model
  cv.fit.GR.high.tcga = cv.glmnet(x=t(expression_for_proliferation[pheno.inx.tcga,]), y=actual_prol, alpha=1,  family="gaussian",nfold=5, standardize = F) # LASSO regression for proliferation model
  tcga.motility =  qnorm.array(predict(cv.fit.motility.tcga, newx=t(mRNA.norm[tcga.inx,]), s="lambda.min"))
  tcga.GR.high =  qnorm.array(predict(cv.fit.GR.high.tcga, newx=t(mRNA.norm[tcga.inx,]), s="lambda.min"))
  tcga.motility_store[i,] = tcga.motility
  tcga.GR.high_store[i,] = tcga.GR.high
}
predicted_migration = apply(tcga.motility_store,2,median) ## predicted migration for TCGA breast cancer patients
predicted_proliferation = apply(tcga.GR.high_store,2,median) ## predicted proliferation for TCGA breast cancer patients



# Ancestry_Voting_for_Sea_Urchins
These scripts were created to show how we were able to run Ancestry Voting in The Wray-McClay lab to identify developmental trajectories in Sea Urchins. We are extrimely grateful with Jonas Brandenburg who shared his orginal scripts with us.


This is an standalone R script to implement ancestry voting, as in Qiu et al (2022) > https://www.nature.com/articles/s41588-022-01018-x

The following pipelines will allow you to build a developmental cladogram like the following in Sea Urchins

![image](https://user-images.githubusercontent.com/5439367/206006293-870215a2-83bd-409f-aca8-1d9e795a05aa.png)


To run using SLURM in an HPC, run the following job launcher:



```bash
#! /bin/bash -l
#SBATCH -J Anc2
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ab620@duke.edu
#SBATCH --mem 50G
cd /data/wraycompute/alejo/singlecell/ancestry_voting2
Rscript ancestry.R
```

But first, make sure you have an R script called ancestry.R with the following code. 

note: Remember to modify the hard code to fit your experiment


```R
#code adapted from https://github.com/ChengxiangQiu/tome_code
#author of RScript: Jonas Brandenburg, implemented and adapted to sea urchin data by Alejandro Berrio 

#packages
library(Seurat)
library(tidyverse)
library(monocle3)
require(lattice)
library(FNN)
library(viridis)
library(reshape)
library(doParallel)


#Create a directory location to save data
inDir=""
outDir <- paste0(inDir,  "AncestryVoting4/")
dir.create(outDir)

#functions 
createLineage_Knn <- function(emb, pd, reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5){
  
  print(dim(emb))
  if(!"Reclustered" %in% names(pd) | !"Stage" %in% names(pd)) {print("Error: no Anno or day in pd")}
  if(sum(rownames(pd)!=rownames(emb))!=0) {print("Error: rownames are not matched")}
  pd$state = pd$Reclustered
  
  res = list()
  
  rep_i = 1
  
  while(rep_i < (replication_times+1)){
    print(rep_i) #to monitor progress. 
    sampling_index = sample(1:nrow(pd),round(nrow(pd)*(1-removing_cells_ratio)))
    
    emb_sub = emb[sampling_index,]
    pd_sub = pd[sampling_index,]
    
    irlba_pca_res_1 <- emb_sub[as.vector(pd_sub$Stage)==anc,]
    irlba_pca_res_2 <- emb_sub[as.vector(pd_sub$Stage)==off,]
    pd_sub1 <- pd_sub[pd_sub$Stage == anc,]
    pd_sub2 <- pd_sub[pd_sub$Stage == off,]
    
    pre_state_min = min(table(as.vector(pd_sub1$state)))
    
    if (pre_state_min < k_neigh & pre_state_min >= 3){
      k_neigh = pre_state_min
      print(k_neigh)
    }
    
    if (pre_state_min < 3){
      next
    }
    
    neighbors <- get.knnx(irlba_pca_res_1, irlba_pca_res_2, k = k_neigh)$nn.index
    
    tmp1 <- matrix(NA,nrow(neighbors),ncol(neighbors))
    for(i in 1:k_neigh){
      tmp1[,i] <- as.vector(pd_sub1$state)[neighbors[,i]]
    }
    state1 <- names(table(as.vector(pd_sub1$state)))
    state2 <- names(table(as.vector(pd_sub2$state)))
    
    tmp2 <- matrix(NA,length(state2),length(state1))
    for(i in 1:length(state2)){
      x <- c(tmp1[as.vector(pd_sub2$state)==state2[i],])
      for(j in 1:length(state1)){
        tmp2[i,j] <- sum(x==state1[j])
      }
    }
    tmp2 <- tmp2/apply(tmp2,1,sum)
    tmp2 <- data.frame(tmp2)
    row.names(tmp2) = state2
    names(tmp2) = state1
    
    res[[rep_i]] = tmp2
    
    rep_i = rep_i + 1
    
  }
  
  return(res)
}

load("dataset_single_cell.Rda", verbose=T)


control_2to24$group <- factor(control_2to24$group,  levels = c("2_HPF", "3_HPF","4_HPF","5_HPF","6_HPF", "7_HPF","8_HPF",
                                                               "9_HPF", "10_HPF" , "11_HPF" , "12_HPF" ,
                                                               "13_HPF", "14_HPF" ,  "15_HPF", "16_HPF" ,
                                                               "18_HPF", "20_HPF", "24_HPF"))

control_2to24$Stage <- control_2to24$group  




#load Data, should be all in one merged seurat object. 
#mergedRNA <- readRDS(file = paste0(inDir, "../mergedRNA.rds"))
mergedRNA <- control_2to24
#reorder developmental time
mergedRNA$Stage<- factor(mergedRNA$Stage, levels = c("2_HPF", "3_HPF","4_HPF","5_HPF","6_HPF",
                                                     "7_HPF","8_HPF",
                                                     "9_HPF", "10_HPF" ,
                                                     "11_HPF" , "12_HPF" ,
                                                     "13_HPF", "14_HPF" ,
                                                     "15_HPF", "16_HPF" ,
                                                     "18_HPF", "20_HPF", "24_HPF" ))




Idents(mergedRNA) 
mergedRNA@active.ident

# here you can assigne names for each cluster in the Seurat object

mergedRNA <- RenameIdents(mergedRNA, `0` = "xxx0", `1` = "xxx1", `2` = "xxx2",
                                  `3` = "xxx3", `4` = "xxx4", `5` = "xxx5", `6` = "xxx6", `7` = "xxx7", `8` = "xxx8", `9` = "xxx9",
                                  `10` = "xxx10", `11` = "xxx11", `12` = "xxx12", `13` = "xxx13", `14` = "xxx14", `15` = "xxx15", `16` = "xxx16")


mergedRNA$Reclustered <- Idents(mergedRNA)

mergedRNA@meta.data #likely phenodata of integrated data 
head(mergedRNA[[]])





ancestry.mappings.list <- list()
ancestries<- list()

#registerDoParallel(56)
#now run ancestry voting on dataset. 
for (i in 2:length(levels(mergedRNA$Stage))){ 
  anc <- levels(mergedRNA$Stage)[i-1]
  off <- levels(mergedRNA$Stage)[i]

  ancestor<- subset(mergedRNA, Stage == anc)
  offspring<- subset(mergedRNA, Stage == off)
  if(dim(ancestor)[2] < 100 | dim(offspring)[2] < 100) {
    next
  } else {
    
    ancestry.list <- list(ancestor, offspring)
    ancestry.anchors <- FindIntegrationAnchors(object.list = ancestry.list, dims = 1:30)
    ancestry.integrated <- IntegrateData(anchorset = ancestry.anchors, normalization.method = 'LogNormalize') ##for RNA
    ancestry.integrated<- ScaleData(ancestry.integrated)
    ancestry.integrated<- RunPCA(ancestry.integrated)
    
    #next calculate w subsampling 
    pd<- ancestry.integrated@meta.data # phenodata of integrated data 
    emb<- Embeddings(object = ancestry.integrated[["pca"]]) # embedding of integrated data 
    
    ancestry.list.tmp <- createLineage_Knn(emb, pd,  reduction="umap", replication_times=500, removing_cells_ratio=0.2, k_neigh = 5)
    
    #calculate median
    for(j in 1:nrow(ancestry.list.tmp[[1]])){
      for(k in 1:ncol(ancestry.list.tmp[[1]])){
        if(k == 1){
          tmp<- median(unlist(lapply(ancestry.list.tmp, function(x) x[j,k])))
          vec <- tmp
        } else {
          tmp<- median(unlist(lapply(ancestry.list.tmp, function(x) x[j,k])))
          vec <- c(vec, tmp)
        }
      }
      if(j == 1){
        mat <- vec
      } else {
        mat <- cbind(mat, vec)
      }
    }
    rownames(mat) <- colnames(ancestry.list.tmp[[1]])
    colnames(mat) <- rownames(ancestry.list.tmp[[1]])
    ancestries[[i-1]] <- mat
    names(ancestries)[i-1] <- off
  }
print(paste0("Current ancestry size: ",ancestries))
print(paste0("Current ancestry column: ",i))

}

#write.csv(ancestries, file = paste0(outDir, "ancestry.csv"), quote = F,  row.names = F, col.names = F)
save(ancestries,file="ancestries.Rda")


#good cutoff would be > 0.3!
#transform to csv, with following meta: 
##Parent, Children, Celltype, GermLayer, Support, Stage
filtered.list <- list()
for(i in 1:9){
  tmp <- reshape::melt(ancestries[i])
  colnames(tmp)<- c('Parent', 'Child', 'Support', 'Stage')
  tmp <- subset(tmp, Support> 0.3)
  filtered.list [[i]]<- tmp
}
filtered.list<- do.call(rbind, filtered.list)
write.csv(filtered.list, file = paste0(outDir, "filtered.ancestry.csv"), quote = F, row.names = F, col.names = F)

#plot heatmaps 
pdf(paste0(outDir, "AncestriesHeatmap.pdf"))
lapply(ancestries, function(x) heatmap(t(x), scale= "none", Rowv = NA, Colv = NA))
dev.off()

```




Once this job has been completed, there should be two files in AncestryVoting4 directory:
AncestriesHeatmap.pdf
and,
filtered.ancestry.csv

Which will have all the cell trajectories in a csv table that includes child to parent information:


|Parent|Child|Support|Stage|
| --------- | --------- | --------- | --------- |
|xxx15|xxx15|0.645792483660131|3_HPF|
|xxx8|xxx15|0.354207516339869|3_HPF|
|xxx16|xxx16|0.663450834879406|3_HPF|
|xxx8|xxx16|0.336126343496381|3_HPF|
|xxx8|xxx8|0.997613365155131|3_HPF|
|xxx15|xxx15|1|4_HPF|
|xxx16|xxx16|0.975455147302218|4_HPF|
|xxx8|xxx3|0.7|4_HPF|
|xxx8|xxx1|1|5_HPF|


if (1<0){
  install.packages('foreach')
  install.packages('doParallel')
  install.packages('flashClust')
  install.packages('dplyr')
  install.packages('ff')
  install.packages('WGCNA')
  install.packages('stringr')
  install.packages('reshape2')
  install.packages('caret')
  install.packages('blob') 
}
if(1>0){
  library(WGCNA)
  library(stringr)
  library(reshape2)
  library(caret)
  #library(dismay)
  library(foreach)
  library(doParallel)
  library(flashClust)
  library(dplyr)
  library(ff)
  bigcor <- function(
    x, 
    y = NULL,
    fun = c("cor", "cov"), 
    size = 2000, 
    verbose = TRUE, 
    ...)
  {
    fun <- match.arg(fun)
    if (fun == "cor") FUN <- cor else FUN <- cov
    if (fun == "cor") STR <- "Correlation" else STR <- "Covariance" 
    if (!is.null(y) & NROW(x) != NROW(y)) stop("'x' and 'y' must have compatible dimensions!")
    
    NCOL <- ncol(x)
    if (!is.null(y)) YCOL <- NCOL(y)
    
    ## calculate remainder, largest 'size'-divisible integer and block size
    REST <- NCOL %% size
    LARGE <- NCOL - REST  
    NBLOCKS <- NCOL %/% size
    
    ## preallocate square matrix of dimension
    ## ncol(x) in 'ff' single format
    if (is.null(y)) resMAT <- ff(vmode = "double", dim = c(NCOL, NCOL))  
    else resMAT <- ff(vmode = "double", dim = c(NCOL, YCOL))
    
    ## split column numbers into 'nblocks' groups + remaining block
    GROUP <- rep(1:NBLOCKS, each = size)
    if (REST > 0) GROUP <- c(GROUP, rep(NBLOCKS + 1, REST))
    SPLIT <- split(1:NCOL, GROUP)
    
    ## create all unique combinations of blocks
    COMBS <- expand.grid(1:length(SPLIT), 1:length(SPLIT))
    COMBS <- t(apply(COMBS, 1, sort))
    COMBS <- unique(COMBS)  
    if (!is.null(y)) COMBS <- cbind(1:length(SPLIT), rep(1, length(SPLIT)))
    
    ## initiate time counter
    timeINIT <- proc.time() 
    
    ## iterate through each block combination, calculate correlation matrix
    ## between blocks and store them in the preallocated matrix on both
    ## symmetric sides of the diagonal
    for (i in 1:nrow(COMBS)) {
      COMB <- COMBS[i, ]    
      G1 <- SPLIT[[COMB[1]]]
      G2 <- SPLIT[[COMB[2]]]    
      
      ## if y = NULL
      if (is.null(y)) {
        if (verbose) cat(sprintf("#%d: %s of Block %s and Block %s (%s x %s) ... ", i, STR,  COMB[1],
                                 COMB[2], length(G1),  length(G2)))      
        RES <- FUN(x[, G1], x[, G2], ...)
        resMAT[G1, G2] <- RES
        resMAT[G2, G1] <- t(RES) 
      } else ## if y = smaller matrix or vector  
      {
        if (verbose) cat(sprintf("#%d: %s of Block %s and 'y' (%s x %s) ... ", i, STR,  COMB[1],
                                 length(G1),  YCOL))    
        RES <- FUN(x[, G1], y, ...)
        resMAT[G1, ] <- RES             
      }
      
      if (verbose) {
        timeNOW <- proc.time() - timeINIT
        cat(timeNOW[3], "s\n")
      }
      
      gc()
    } 
    
    return(resMAT)
  }
}


setwd('C:/Users/user/Desktop/test/scanpy')
foldername='pbmc3k'
Clusteringmethod='Leiden'
clustering_size=read.csv('/',foldername,'/',Clusteringmethod,'_clustering_size.csv', check.names=FALSE)


for(i in 1:nrow(clustering_size))
{
  cat("round ",i,"\n")
  clusterName <- clustering_size[i,1]# = 第i個cluster的檔案名稱
  clusterCellNumber <- clustering_size[i,2]# = 第i個cluster的cell數量
  dir.create(paste("./",foldername,"/",clusterName,sep=""))#建立該cluster的資料夾
  dir.create(paste("./",foldername,"/",clusterName,"/modules",sep=""))#建立儲存該cluster module gene的資料夾
  ## SAVE!
  cat(as.character(clusterName),'\t',clusterCellNumber,'\n',sep=" ",file=(paste("./",foldername,"/",clusterName,"/info.txt",sep="")),append = T)
  
  #input RNAseq的檔案，一個column是一個基因，一個row是一個cell
  options(stringsAsFactors = FALSE)
  Data = read.csv(paste("./",foldername,"/",clusterName,".csv",sep=""), check.names=FALSE)
  rownames(Data) <- c(Data[,1])
  Data <- Data[,-1]#將barcode移轉到rownames
  #查看資料
  Data[1:10,1:10]
  ## SAVE!
  #cat('Before MAD, data dimension:',dim(Data),'\n',sep=" ",file=(paste("./",foldername,"/",clusterName,"/info.txt",sep="")),append = T)
  
  
  ##這邊做一次PCC
  
  #Pcc = bigcor(as.matrix(Data))
  #dim(Pcc)
  #Pcc
  #write.table(Pcc,paste("./",foldername,"/",clusterName,"/",clusterName,"AfterMAD_Pcc.csv",sep=""),row.names=T, na = "NA")    
  #Calculate the p-values of the pearson correlation
  #corPvalueFisher(cor, nSamples, twoSided = TRUE)
  #Pcc_pvalue = corPvalueFisher(Pcc, nrow(Data), twoSided = TRUE)
  #dim(Pcc_pvalue)
  #Pcc_pvalue[1:4,1:4]
  
  #取絕對中位數偏差(MAD)前75%的基因，並且MAD至少要大於0.01
  dataExpr <- as.data.frame(t(Data))
  m.mad <- apply(dataExpr,1,mad)
  dataExprVar <- dataExpr[which(m.mad > 
                                  max(quantile(m.mad, probs=seq(0, 1, 0.25))[2],0.01)),]
  #因為MAD那步轉了一次矩陣，所以再轉一次
  dataExpr <- as.data.frame(t(dataExprVar))
  ###ataExpr <- as.data.frame(t(Data))
  ##再做一次PCC
  
  #Pcc = bigcor(as.matrix(dataExpr))
  #dim(Pcc)
  #Pcc[1:4,1:4]
  #write.table(Pcc,paste("./",foldername,"/",clusterName,"/",clusterName,"AfterMAD_Pcc.csv",sep=""),row.names=T, na = "NA")
  #Calculate the p-values of the pearson correlation
  #corPvalueFisher(cor, nSamples, twoSided = TRUE)
  #Pcc_pvalue = corPvalueFisher(Pcc, nrow(dataExpr), twoSided = TRUE)
  #dim(Pcc_pvalue)
  #Pcc_pvalue[1:4,1:4]
  
  #篩選掉變異數為0或太多0的基因
  gsg = goodSamplesGenes(dataExpr, verbose = 3)
  #印出被篩選掉的基因
  if (!gsg$allOK){
    # Optionally, print the gene and sample names that were removed:
    if (sum(!gsg$goodGenes)>0)
      printFlush(paste("Removing genes:", 
                       paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
    if (sum(!gsg$goodSamples)>0) 
      printFlush(paste("Removing samples:", 
                       paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
    # Remove the offending genes and samples from the data:
    dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
  }
  
  nGenes = ncol(dataExpr)
  nSamples = nrow(dataExpr)
  
  
  ## SAVE!
  cat('After MAD, data dimension:',dim(dataExpr),'\n',sep=" ",file=(paste("./",foldername,"/",clusterName,"/info.txt",sep="")),append = T)
  
  #計算軟閥值(要讓數值幾次方)
  powers = c(c(1:10), seq(from = 12, to=30, by=2))
  sft = pickSoftThreshold(dataExpr, powerVector=powers, verbose=5)
  #樹狀圖
  ## SAVE!
  sampleTree = hclust(dist(dataExpr), method = "average")
  
  pdf(paste("./",foldername,"/",clusterName,"/sampleTree.pdf",sep=""),width = 20, height = 20)
  plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="")
  dev.off()
  ########
  
  
  ## SAVE!
  pdf(paste("./",foldername,"/",clusterName,"/scaleFree.pdf",sep=""),width = 20, height = 20)
  par(mfrow = c(1,2))
  cex1 = 0.9
  # 橫軸是Soft threshold (power)，縱軸是無標度網絡的評估參數，數值越高，
  # 網絡越符合無標度特徵 (non-scale) 
  
  plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       xlab="Soft Threshold (power)",
       ylab="Scale Free Topology Model Fit,signed R^2",type="n",
       main = paste("Scale independence"))
  text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
       labels=powers,cex=cex1,col="red")
  # 篩選標準。 R-square=0.85
  abline(h=0.85,col="red")
  
  # Soft threshold與平均連通性 
  plot(sft$fitIndices[,1], sft$fitIndices[,5],
       xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
       main = paste("Mean connectivity"))
  text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, 
       cex=cex1, col="red")
  dev.off()
  #########
  
  
  #自動選擇軟閾值數值
  ## SAVE!
  power = sft$powerEstimate
  cat('Power:',power,'\n',sep=' ',file=(paste("./",foldername,"/",clusterName,"/info.txt",sep="")),append = T)
  
  
  
  
  #自動建立網路跟模塊
  net = blockwiseModules(dataExpr, power = power, maxBlockSize = nGenes,
                         TOMType = 'unsigned', minModuleSize = 30,
                         reassignThreshold = 0, mergeCutHeight = 0.25,
                         numericLabels = TRUE, pamRespectsDendro = FALSE,
                         saveTOMs=TRUE, 
                         loadTOMs=TRUE,
                         verbose = 3)
  #模塊數跟基因數 
  ## SAVE!
  write.table(table(net$colors),file=(paste("./",foldername,"/",clusterName,"/info.txt",sep="")),append = T,row.names=F)
  
  
  ## 灰色的為未分類到模塊的基因。
  # Convert labels to colors for plotting
  moduleLabels = net$colors
  moduleColors = labels2colors(moduleLabels)
  # Plot the dendrogram and the module colors underneath
  # 如果對結果不滿意，還可以recutBlockwiseTrees，節省計算時間 
  ## SAVE!
  pdf(paste("./",foldername,"/",clusterName,"/moduleDendrogram.pdf",sep=""),width = 20, height = 20)
  plotDendroAndColors(net$dendrograms[[1]], moduleColors[net$blockGenes[[1]]],
                      "Module colors",
                      dendroLabels = FALSE, hang = 0.03,
                      addGuide = TRUE, guideHang = 0.05, lwd = 2, font=2)
  dev.off()
  
  #可以不用跑，也可以跑，不會太久
  # module eigengene, 可以繪製線圖，作為每個模塊的基因表達趨勢的展示 
  MEs = net$MEs
  if(ncol(MEs)>2)
  {
    #防呆,Module過少就不會執行,若module只有一個的話dendrogram畫不出來
    ### 不需要重新計算，改下列名字就好
    ### 官方教程是重新計算的，起始可以不用這麼麻煩
    MEs_col = MEs
    colnames(MEs_col) = paste0("ME", labels2colors(as.numeric(str_replace_all(colnames(MEs),"ME",""))))
    MEs_col = orderMEs(MEs_col)
    
    # 根據基因間表達量進行聚類所得到的各模塊間的相關性圖
    # marDendro/marHeatmap 設置下、左、上、右的邊距
    ## SAVE!
    if(ncol(MEs)>3)
    {
      pdf(paste("./",foldername,"/",clusterName,"/moduleEigengene.pdf",sep=""),width = 20, height = 20)
      plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap and dendrogram",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = T, 
                            xLabelsAngle = 90)
      dev.off()
    }else
    {
      pdf(paste("./",foldername,"/",clusterName,"/moduleEigengene.pdf",sep=""),width = 20, height = 20)
      plotEigengeneNetworks(MEs_col, "Eigengene adjacency heatmap",
                            marDendro = c(3,3,2,4),
                            marHeatmap = c(3,4,2,2), plotDendrograms = F, 
                            xLabelsAngle = 90)
      dev.off()
    }
  }
  
  #這步可以不用跑，很花時間，有需要圖再跑
  # 如果採用分步計算，或設置的blocksize>=總基因數，直接load計算好的TOM結果
  # 否則需要再計算一遍，比較耗費時間
  # TOM = TOMsimilarityFromExpr(dataExpr, power=power, corType=corType, networkType=type) 
  load(net$TOMFiles[1], verbose=T)
  
  ## Loading objects:
  ##   TOM
  
  TOM <- as.matrix(TOM)
  
  dissTOM = 1-TOM
  # Transform dissTOM with a power to make moderately strong 
  # connections more visible in the heatmap
  plotTOM = dissTOM^7
  # Set diagonal to NA for a nicer plot
  diag(plotTOM) = NA
  # Call the plot function
  ## SAVE!
  # 這一部分特別耗時，行列同時做層級聚類 
  pdf(paste("./",foldername,"/",clusterName,"/TOM.pdf",sep=""),width = 20, height = 20)
  TOMplot(plotTOM, net$dendrograms, moduleColors, 
          main = "Network heatmap plot, all genes")
  dev.off()
  
  
  
  #這步可以先跳過，之後有需要再回來跑，這步是可以設定閥值output出網路圖的步驟
  probes = colnames(dataExpr)
  dimnames(TOM) <- list(probes, probes)
  workdir <- getwd()
  setwd(paste(workdir,"/",foldername,"/",clusterName,"/modules",sep=""))
  
  # Export the network into edge and node list files Cytoscape can read
  # threshold 默認為0.5, 可以根據自己的需要調整，也可以都導出後在
  # cytoscape中再調整 
  cyt = exportNetworkToCytoscape(TOM,
                                 edgeFile = paste(clusterName,"edges.txt",sep=""),
                                 nodeFile = paste(clusterName,"nodes.txt",sep=""),
                                 weighted = TRUE, threshold = 0,
                                 nodeNames = probes, nodeAttr = moduleColors)
  setwd(paste(workdir))
  
  
  
  
  #####################
  #preservation######
  #####################
  
  
  if(ncol(MEs)>1)
  {
    ###new code
    for(j in 1:nrow(clustering_size))
    {
      cat("round ",i,',',j,"\n")
      if(j!=i)
      {
        #主要cluster對所有其他的cluster, 在其他cluster中用下面的test和train
        cat("other cluster\n")
        test = read.csv(paste('./',foldername,'/Leiden_cluster_',j-1,'.csv',sep=""), check.names=FALSE)
        test = data.matrix(test)
        train = dataExpr
        ######
      }else
      {
        #這步是把原本的data分一部分出來當test，所以數字要依cell數量去改(我這邊是242 cells)
        #samople數量/4
        #now 438 cells
        #438=109*4+2
        #m=c(rep(c(1:4),109))  "109*4"
        #n=c(1,2)   "+2"
        cat("same cluster\n")
        train_dataset = as.data.frame(dataExpr)
        m = c(rep(c(1:4),(clusterCellNumber%/%4)))
        n = c(1:(clusterCellNumber%%4))
        if(clusterCellNumber%%4!=0){
          m = c(m,n)
        }
        train_dataset = cbind(m,train_dataset)
        colnames(train_dataset)[1] = "group"
        inTrain = createDataPartition(y=train_dataset$group,p=0.25,list=FALSE)
        train = train_dataset[inTrain,-1]
        test = train_dataset[-inTrain,-1]
        
      }
      #train = as.data.frame(Data)
      
      gsg = goodSamplesGenes(test, verbose = 3);
      gsg$allOK
      if (!gsg$allOK)
      {
        # Optionally, print the gene and sample names that were removed:
        if (sum(!gsg$goodGenes)>0) 
          printFlush(paste("Removing genes:", paste(names(test)[!gsg$goodGenes], collapse = ", ")));
        if (sum(!gsg$goodSamples)>0) 
          printFlush(paste("Removing samples:", paste(rownames(test)[!gsg$goodSamples], collapse = ", ")));
        # Remove the offending genes and samples from the data:
        test = test[gsg$goodSamples, gsg$goodGenes]
      }
      
      setLabels = c("Train", "Test")
      multiExpr = list(Train = list(data = train), Test = list(data = test))
      multiColor = list(Train = moduleColors)
      cat("round ",i,',',j,"mp \n")
      #preservation，npermutations 200就夠了，每次跑出來數字不一樣是正常的，但是不會差太多
      mp = modulePreservation(multiExpr, multiColor,
                              referenceNetworks = 1,
                              nPermutations = 200,
                              randomSeed = 1,
                              quickCor = 0,
                              verbose = 3)
      
      ref=1
      test=2
      modColors = rownames(mp$preservation$observed[[ref]][[test]])
      moduleSizes = mp$preservation$Z[[ref]][[test]][, 1]
      
      #去掉"grey", "gold"
      #%in%不在這裡的基因
      plotMods = !(modColors %in% c("grey", "gold"));
      text = modColors[plotMods]
      plotData = cbind(mp$preservation$observed[[ref]][[test]][, 2], mp$preservation$Z[[ref]][[test]][, 2])
      
      #這邊Zsummary.pres那欄就是Zsummary score
      statsObs = cbind(mp$quality$observed[[ref]][[test]][, -1], mp$preservation$observed[[ref]][[test]][, -1])
      statsZ = cbind(mp$quality$Z[[ref]][[test]][, -1], mp$preservation$Z[[ref]][[test]][, -1])
      #儲存module的preservation Zscore table
      ## SAVE!
      #此處是cluster的module在自己的cluster中的preservation
      write.csv(data.frame(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],
                                 signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2))),
                paste("./",foldername,"/",clusterName,"/preservation_c",i-1,"-c",j-1,"_200permutation_Zscore.csv",sep=""))
      #####
      ##收集module在每個cluster的preservation Zscore
      #####
      if (j == 1){
        mpZscore<-as.data.frame(statsZ["Zsummary.pres"])
        #mpZscore 就是 module preservation Z score
      }
      else
        mpZscore <- cbind(mpZscore,as.data.frame(statsZ["Zsummary.pres"]))
      
      ###顯示 z score的 code
      #0要改成可變動變數,指出cluster幾,c0=cluster0,迴圈index:i要-1
      # Compare preservation to quality:
      #print(cbind(statsObs[, c("medianRank.pres", "medianRank.qual")],signif(statsZ[, c("Zsummary.pres", "Zsummary.qual")], 2)) )
      ###
      
      ##preservation可視化，一般來說只需要改檔名路徑
      mains = c("Preservation Median rank", "Preservation Zsummary")
      #檔名
      ## SAVE!
      pdf(paste("./",foldername,"/",clusterName,"/preservation_c",i-1,"-c",j-1,"_200permutation_Zscore.pdf",sep=""),width = 20, height = 10)
      ##一行兩列
      par(mfrow = c(1,2))
      ##到四邊的距離
      par(mar = c(4.5,4.5,2.5,1))
      for (p in 1:2)
      {
        min = min(plotData[, p], na.rm = TRUE);
        max = max(plotData[, p], na.rm = TRUE);
        # Adjust ploting ranges appropriately
        if (p==2)
        {
          if (min > -max/10) min = -max/10
          ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
        } else
          ylim = c(min - 0.1 * (max-min), max + 0.1 * (max-min))
        #bg 顏色 pch 圓圈的種類
        plot(moduleSizes[plotMods], plotData[plotMods, p], col = 1, bg = modColors[plotMods], pch = 21,
             main = mains[p],
             ##圓圈的大小
             cex = 2.4,
             ylab = mains[p], xlab = "Module size", log = "x",
             ylim = ylim,
             xlim = c(10, 2000), cex.lab = 1.2, cex.axis = 1.2, cex.main =1.4, font.lab=2)
        ##貼上標籤
        labelPoints(moduleSizes[plotMods], plotData[plotMods, p], text, cex = 1, offs = 0.08, font.lab=2);
        # For Zsummary, add threshold lines
        if (p==2)
        {
          abline(h=0)
          abline(h=2, col = "blue", lty = 2)
          abline(h=10, col = "darkgreen", lty = 2)
        }
      }
      box(which = "plot", col = "black", lwd = 5)
      dev.off()
      
    }
    ##儲存module在每個cluster 的preservation Zscore
    ## SAVE!
    mpZscore <- t(mpZscore)
    rownames(mpZscore) <- c(0:(nrow(clustering_size)-1))
    write.csv(mpZscore,paste("./",foldername,"/",clusterName,"/module_preservation_Zscore.csv",sep = ""))
    
    #輸出每個module的基因，在程式一開頭已經創好資料夾,在該cluster的modlue資料夾中操作
    ## SAVE!
    for(j in 1:length(text))
    {
      y=t(dataExpr)[which(moduleColors==text[j]),]
      write.csv(y,paste("./",foldername,"/",clusterName,"/modules/",text[j],".csv",sep = ""),quote=F)
    }
  }else
  {
    grey <- c(rep(0,nrow(clustering_size)))
    gold <- c(rep(0,nrow(clustering_size)))
    blankdataframe <- data.frame(gold, grey)
    rownames(blankdataframe) <- paste(c(0:(nrow(blankdataframe)-1)))
    ##要插入一個空白的preservation z score
    ##column只有grey 和gold, row =cluster數
    #blankdataframe=
    write.csv(blankdataframe,paste("./",foldername,"/",clusterName,"/module_preservation_Zscore.csv",sep = ""))         
    cat('no module.\n')
  }
}
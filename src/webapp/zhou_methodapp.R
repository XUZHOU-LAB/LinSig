library(DT)
library(ggplot2)
library(ggpubr)
library(dplyr)
library(BiocManager)
library(circlize)
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
library(RColorBrewer)
library(VennDiagram)
library(ggvenn)

#TODO
#optional qvalues
#multiple replicates
#clearer, more concise UI
# colnames error message

ui = fluidPage(
  sidebarPanel(
    fileInput("df", "Upload Counts, q-values",
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    
    checkboxInput("lowess", "LOWESS Norm", value=FALSE),
    
    numericInput("pseudo", "Add Pseudo Count",
                 value=1,
                 min=0.1,
                 step=0.5),
    numericInput("cntThres","Count Threshold",
                 value=10,
                 min=0,
                 step=0.1),
    #  actionButton("normalise","Compute Normalized Ratios"),
    hr(),
    sliderInput("lfcThres", "LFC Threshold (LFC>X & qval<0.05)",
                value=0.585, min=0, max=3, step=0.001),
    actionButton("deconvolute", "Deconvolute Signals"),
    hr(),
    sliderInput("H0Thres", "Null-Hypothesis Test FC Threshold",
                value=1.5, min=1, max=4, step=0.01),
    sliderInput("R2Thres", "R2 Threshold",
                value=0.8, min=0.1, max=1, step=0.01),
    actionButton("compFDR", "Compute False Discovery Rate"),
    numericInput("RDFsize", "Random Dataset Size",
                 value=40000, min=10000, step=1),
    sliderInput("FCFDR", "FC Selection for FDR",
                value=1, min=1, max=5, step=0.01),
    downloadButton('downloadData', 'Download Data')
  ),
  
  
  mainPanel(
    tabsetPanel(type="tabs",
                tabPanel("Data Table", dataTableOutput("stats"), textOutput("files"),
                         plotOutput("venn")),
                tabPanel("R2~Beta Plots", plotOutput("betaR2")),
                tabPanel("FDR",
                         fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"), plotOutput("betaR2_FDR"), plotOutput("FDRdistribution"))
                         ),
                         #plotOutput("betaR2_FDR"),
                         #plotOutput("FDRdistribution"),
                         verbatimTextOutput("cFDR"),
                         verbatimTextOutput("pval")),
                tabPanel("Volcano Plots",
                         selectInput("plottype", "Choose a plot:",
                                     choices = c("1", "2", "3")),
                         plotOutput('volcanoPlot', click='plot_click',
                                    brush='plot_brush'),
                         tableOutput('clickedPoints')),
                tabPanel("Heatmap",
                         actionButton("show_heatmap", "Generate Heatmap"),
                         htmlOutput("heatmap_output"))
    )
  )
)

################################################################################

geomean <- function(x){
  exp(mean(log(x)))
}

MedianNorm <- function(Data, CountThres=10, pseudo=1){
  
  Subdata = Data[rowMeans(Data)>=CountThres,]     # Subset rows with average count higher than 10
  Subdata = Subdata + pseudo                       # Add Pseudo count of 1
  
  t = apply(Subdata,1, geomean)         # Calculate Geometric mean
  # Divide expression of every element per row to the geometric mean of that row
  # Every 8 elements of Data need to be divided by an element of t
  
  DataRatio = log2(Subdata / t)
  T_med = apply(DataRatio,2,median) # median per column of array
  T_med = 2^T_med
  #print(T_med)
  C <- t(t(Subdata) / T_med)
  return(C)
}

RNAseqLowess <- function(LogIntensity, LogRatios){
  
  Ynorm = loess(LogRatios~LogIntensity,
                span=0.05,degree=1,family="gaussian",
                iterations=1,surface="direct")
  Ratiosnorm = LogRatios - fitted(Ynorm)
  
  
  
  return(Ratiosnorm)
}
# Function for normalization with LOWESS
normalize <- function(inputdf, cntThres, pseudo, lowess=TRUE){
  
  print("Applying Normalisation")
  cntMat <- as.matrix(inputdf[,1:8])
  DataPseudo <- MedianNorm(cntMat, CountThres=cntThres, pseudo=pseudo)
  
  LogRatios <- matrix(nrow=length(DataPseudo[,1]), ncol=10) #empty matrix
  LogRatios[,1] <- log2(DataPseudo[,3] / DataPseudo[,1])
  LogRatios[,2] <- log2(DataPseudo[,7] / DataPseudo[,5])
  LogRatios[,3] <- log2(DataPseudo[,5] / DataPseudo[,1])
  LogRatios[,4] <- log2(DataPseudo[,7] / DataPseudo[,3])
  LogRatios[,5] <- log2(DataPseudo[,7] / DataPseudo[,1])
  LogRatios[,6] <- log2(DataPseudo[,4] / DataPseudo[,2])
  LogRatios[,7] <- log2(DataPseudo[,8] / DataPseudo[,6])
  LogRatios[,8] <- log2(DataPseudo[,6] / DataPseudo[,2])
  LogRatios[,9] <- log2(DataPseudo[,8] / DataPseudo[,4])
  LogRatios[,10] <-log2(DataPseudo[,8] / DataPseudo[,2])
  
  LogIntensity <- matrix(nrow=length(DataPseudo[,1]), ncol=10) #empty matrix
  LogIntensity[,1] <- log2(DataPseudo[,3] * DataPseudo[,1])
  LogIntensity[,2] <- log2(DataPseudo[,7] * DataPseudo[,5])
  LogIntensity[,3] <- log2(DataPseudo[,5] * DataPseudo[,1])
  LogIntensity[,4] <- log2(DataPseudo[,7] * DataPseudo[,3])
  LogIntensity[,5] <- log2(DataPseudo[,7] * DataPseudo[,1])
  LogIntensity[,6] <- log2(DataPseudo[,4] * DataPseudo[,2])
  LogIntensity[,7] <- log2(DataPseudo[,8] * DataPseudo[,6])
  LogIntensity[,8] <- log2(DataPseudo[,6] * DataPseudo[,2])
  LogIntensity[,9] <- log2(DataPseudo[,8] * DataPseudo[,4])
  LogIntensity[,10] <-log2(DataPseudo[,8] * DataPseudo[,2])
  
  NumRatios <- length(LogRatios[1,])
  
  if(lowess==FALSE){
    rownames(LogRatios) <- rownames(DataPseudo)
    return(LogRatios)
  }
  
  withProgress(message="Computing Normalized Ratios", value=0,{
    
    Ratios <- matrix(ncol=10, nrow=length(LogRatios[,1]))
    for (i in 1:NumRatios){
      Ratios[,i] <- RNAseqLowess(LogIntensity[,i], LogRatios[,i])
      incProgress(1/NumRatios)
    }
  })
  rownames(Ratios) <- rownames(DataPseudo)
  return(Ratios)
  
}

deconvolutef <- function(ratiosDF, countDF){
  #Ratios <- normRatios()[firstFilter(),]
  Ratios <- ratiosDF
  strucvec <- c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,1)
  strucmat <- matrix(strucvec, ncol=3, byrow=T)
  X <- strucmat
  
  A = length(Ratios[,1]) # of gene
  B = matrix(0, nrow=A,ncol=4);
  
  #get column data
  cc <- strsplit(colnames(countDF)[1:8], "_")
  cols <- unique(unlist(cc)[2*(1:length(cc))-1])[-1]
  colnames(X) <- c(cols[2], cols[1], cols[3])
  
  
  rsquare <- vector()
  ngen <- length(Ratios[,1])
  covB_l <- NULL
  covB <- data.frame(cov_int=double(ngen), cov_x1=double(ngen), cov_x2=double(ngen), cov_x3=double(ngen))
  
  #adds progress bar
  withProgress(message="Computing Model Statistics", value=0,{
    
    print(ngen)
    for (i in 1:ngen){
      fit <- lm(Ratios[i,] ~ X)
      covB[i,] <- diag(vcov(fit))
      rsquare[i] <- summary(fit)$r.squared
      incProgress(1/ngen)
    }
  })
  
  multiplefit <- lm(t(Ratios) ~ X)
  B <- t(multiplefit$coefficients)
  Resid <- multiplefit$residuals
  
  modelStats <- as.data.frame(cbind(B, covB))
  
  colnames(modelStats) <- c("int", cols[2], cols[1], cols[3], "p_int",
                            paste0("cov_", cols[2]), paste0("cov_", cols[1]), paste0("cov_", cols[3]))
  modelStats$R2 <- rsquare
  modelStats <- round(modelStats, digits=7)
  return(modelStats)
}

#compute P value with threshold
calcPvalue <- function(Betas, covB, Threshold){
  ttest_stat <- (abs(Betas) - log2(Threshold)) / sqrt(covB)
  ttest_stat <- data.frame(ttest_stat)
  Pvalue = 1 - apply(ttest_stat, 2, pt, df=4)
  return(Pvalue)
}

NearestNeighbours <- function(x,y,i){ #x:betas, y:FDR, i:numbertotest
  loc <- which.min(abs(x-abs(i)))
  return(y[loc])
}

################################################################################

server = function(input,output, session){
  
  options(shiny.maxRequestSize=20*1024^2) #20MB max file size
  inputfile <- reactive({
    req(input$df)
    df <- read.csv(input$df$datapath, header=T,row.names=1)
    return(df)
  })
  
  
  normRatios <- eventReactive(input$deconvolute, {
    return(normalize(inputfile(), cntThres=input$cntThres, pseudo=input$pseudo,
                     lowess=input$lowess))
  })
  # Function that filters based on q-value and Fold Change
  firstFilter <- reactive({
    cntMat <- as.matrix(inputfile()[,1:8])
    
    DataPseudo <- MedianNorm(cntMat, CountThres=input$cntThres, pseudo=input$pseudo)
    F_AvsC = log2((DataPseudo[,3] + DataPseudo[,4]) / (DataPseudo[,1] + DataPseudo[,2]))
    F_ABvsB = log2((DataPseudo[,7] + DataPseudo[,8]) / (DataPseudo[,5] + DataPseudo[,6]))
    F_ABvsA = log2((DataPseudo[,7] + DataPseudo[,8]) / (DataPseudo[,3] + DataPseudo[,4]))
    F_BvsC = log2((DataPseudo[,5] + DataPseudo[,6]) / (DataPseudo[,1] + DataPseudo[,2]))
    
    Fold <- data.frame(F_AvsC, F_ABvsB, F_ABvsA, F_BvsC)
    
    qvals <- inputfile()[,9:12]
    countthresholdFilter <- rownames(DataPseudo)
    qvalDF <- qvals[rownames(qvals) %in% countthresholdFilter,]
    
    SigIDX = rowSums(abs(Fold[,1:4]) >= input$lfcThres & qvalDF <= 0.05)>0
    
    return(SigIDX)
  })
  
  deconvolute <- eventReactive(input$deconvolute, {
    deconvolutef(normRatios()[firstFilter(),], inputfile())
  })
  
  #returns Pvalues based on LFC and COVs
  Pvalues_real <- reactive({
    real_stats <- deconvolute()
    p <- calcPvalue(real_stats[,1:4], real_stats[,5:8], input$H0Thres)
    print(paste("real", sum(rowSums(p<=0.05)>0)))
    print(sum(rowSums(p>0.5)>0))
    colnames(p) <- paste0("p_", colnames(p))
    return(p)
    
  })
  
  sigReal <- reactive({
    Pvalues <- Pvalues_real()
    real_stats <- deconvolute()
    r2 <- real_stats$R2
    FCFDR <- log2(input$FCFDR)
    A <- (abs(real_stats[,2])>FCFDR & r2>input$R2Thres & Pvalues[,2]<0.05)
    B <- (abs(real_stats[,3])>FCFDR & r2>input$R2Thres & Pvalues[,3]<0.05)
    AB <- (abs(real_stats[,4])>FCFDR & r2>input$R2Thres & Pvalues[,4]<0.05)
    print(paste("real A|B|AB:", sum(A|B|AB)))
    print(paste("A B AB", sum(A), sum(B), sum(AB)))
    
    return(list(A, B, AB))
  })
  
  #returns LFC distribution of deconvoluted random genes
  compFalseDisc <- eventReactive(input$compFDR, {
    #normalize countdata
    cntMat <- as.matrix(inputfile()[,1:8])
    normedCounts <- MedianNorm(cntMat, CountThres=input$cntThres, pseudo=input$pseudo)
    
    #compute mean/sd per condition
    mean_per_condition <- matrix(ncol=4,nrow=nrow(normedCounts))
    mean_per_condition[,1] <- rowMeans(normedCounts[, c(1,2)])
    mean_per_condition[,2] <- rowMeans(normedCounts[, c(3,4)])
    mean_per_condition[,3] <- rowMeans(normedCounts[, c(5,6)])
    mean_per_condition[,4] <- rowMeans(normedCounts[, c(7,8)])
    
    sd_per_condition <- matrix(ncol=4,nrow=nrow(normedCounts))
    sd_per_condition[,1] <- sqrt(rowSums((rowMeans(normedCounts[,c(1,2)]) - normedCounts[,c(1,2)])**2))
    sd_per_condition[,2] <- sqrt(rowSums((rowMeans(normedCounts[,c(3,4)]) - normedCounts[,c(3,4)])**2))
    sd_per_condition[,3] <- sqrt(rowSums((rowMeans(normedCounts[,c(5,6)]) - normedCounts[,c(5,6)])**2))
    sd_per_condition[,4] <- sqrt(rowSums((rowMeans(normedCounts[,c(7,8)]) - normedCounts[,c(7,8)])**2))
    
    mucov<-data.frame(mu=array(mean_per_condition), cov=array(sd_per_condition/mean_per_condition))
    lmucov <- log(mucov)
    olmucov <- lmucov[order(lmucov$mu),]
    binmucov <- olmucov %>% mutate(mu_bin = cut(mu, breaks=20))
    covsList <- c()
    sizeRDF <- input$RDFsize
    #sample cov between each bin
    for (i in unique(binmucov$mu_bin)){
      nsamp <- sizeRDF
      bin_size <- nrow(binmucov[binmucov$mu_bin==i,])
      newsamplesize <- round((bin_size/nrow(binmucov))*nsamp)
      covs <- sample(binmucov[binmucov$mu_bin==i,]$cov, newsamplesize, replace=T)
      covsList <- c(covsList, covs)
    }
    
    #compute distribution of gene means
    mustat <-MASS::fitdistr(rowMeans(normedCounts), "lognormal")
    new_row_means <- rlnorm(sizeRDF, mustat[[1]][1], mustat[[1]][2])
    onrm <- new_row_means[order(new_row_means)]
    sampledMUCOV <- data.frame(cbind(log(onrm), covsList))
    colnames(sampledMUCOV) <- c("mu", "cov")
    #generate new counts
    drawreps <- function(mucovr, rep=4*2){
      rnorm(n=rep, mean=exp(mucovr[1]), sd=exp(mucovr[1]+mucovr[2]))
    }
    RandomDF <- t(round(apply(sampledMUCOV,1, FUN=drawreps), 1))
    RandomDF[RandomDF<0] <- 1 # convert negative counts to 1
    
    RandomDF <- data.frame(RandomDF)
    colnames(RandomDF) <- c("c_A","c_B", "A_A","A_B", "B_A", "B_B", "AB_A", "AB_B")
    
    normalizedRandomDF <- normalize(RandomDF, cntThres=input$cntThres, pseudo=input$pseudo, lowess=input$lowess)
    mstats <- deconvolutef(normalizedRandomDF, RandomDF)
    
    mstats[,5:8] <- calcPvalue(mstats[,1:4], mstats[5:8], 1)
    colnames(mstats)[5:8] <- c("p_int", "p_B", "p_A", "p_AB")
    
    mstats$SIGB <- mstats[,6]<0.05 & mstats$R2>0.8 # 0.7 for adjusted R2, 0.8 for Multiple Rsquared
    mstats$SIGA <- mstats[,7]<0.05 & mstats$R2>0.8
    mstats$SIGAB <-mstats[,8]<0.05 & mstats$R2>0.8
    FDRA <- mstats[mstats$SIGA,]
    FDRB <- mstats[mstats$SIGB,]
    FDRAB<- mstats[mstats$SIGAB,]
    
    FDRab <- FDRb <- FDRa <- c()
    for (i in 1:length(seq(0,4,0.001))){
      FDRab[i] <- mean(abs(FDRAB$AB)>(seq(0,4,0.001)[i]))
    }
    for (i in 1:length(seq(0,4,0.001))){
      FDRb[i] <- mean(abs(FDRB$B)>(seq(0,4,0.001)[i]))
    }
    for (i in 1:length(seq(0,4,0.001))){
      FDRa[i] <- mean(abs(FDRA$A)>(seq(0,4,0.001)[i]))
    }
    xs <- seq(0,4,0.001)
    FDRDF <- data.frame("Xs"=xs, "FDRa"=FDRa,"FDRb"=FDRb, "FDRab"=FDRab)
    print(paste("5% FDR B_threshold A:", seq(0,4,0.001)[which.min(abs(FDRa-0.05))]))
    print(paste("5% FDR B_threshold B:", seq(0,4,0.001)[which.min(abs(FDRb-0.05))]))
    print(paste("5% FDR B_threshold AB:",seq(0,4,0.001)[which.min(abs(FDRab-0.05))]))
    
    return(list(sampledMUCOV, FDRa, FDRb, FDRab))
  })
  
  #returns df with LFC, COV and FDR values
  joinFDRandGenes <- reactive({
    #after deconvolute button
    decoDF <- deconvolute()
    sigReal <- sigReal()
    boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
    decoDF[,c(6,7,8)] <- Pvalues_real()[,c(2,3,4)]
    colnames(decoDF)[6:8] <- colnames(Pvalues_real())[2:4]
    
    #after FDR computation
    FDRa_dist <- compFalseDisc()[[2]]
    FDRb_dist <- compFalseDisc()[[3]]
    FDRab_dist<- compFalseDisc()[[4]]
    xs <- seq(0,4,0.001)
    FDR_A <- sapply(decoDF[,2], NearestNeighbours, x=xs, y=FDRa_dist)
    FDR_B <- sapply(decoDF[,3], NearestNeighbours, x=xs, y=FDRb_dist)
    FDR_AB<- sapply(decoDF[,4], NearestNeighbours, x=xs, y=FDRab_dist)
    decoDF[,10] <- FDR_A
    decoDF[,11] <- FDR_B
    decoDF[,12] <- FDR_AB
    colnames(decoDF)[10:12] <- paste0("FDR_", colnames(decoDF[,2:4]))
    filteredDF <- round(decoDF[boolfilt,c(2,3,4,6,7,8,9,10,11,12)], digits=6)
    
    return(filteredDF)
  })
  
  output$venn <- renderPlot({
    df <- sigReal()
    # use data frame as input
    mstats <- deconvolute()
    cols <- colnames(mstats)[2:4]
    
    M <-tibble(value=1:length(df[[1]]), 'A'= df[[1]],
               'B'= df[[2]],
               'AB'= df[[3]])
    names(M)[2:4] <- cols
    # create Venn diagram and display all sets
    ggvenn(M, fill_color=c("blue","red", "purple"), fill_alpha=0.25)
  })
  
  
  output$betaR2 <- renderPlot({
    mstats <- deconvolute()
    sigReal <- sigReal()
    cols <- colnames(mstats)[2:4]
    R2 <- mstats$R2
    
    plotFun <- function(dat, x, y, sigs){
      ggplot(dat, aes(x=.data[[x]], y=.data[[y]])) +
        geom_point() +
        geom_point(data=dat[sigs,], color="blue")+
        coord_cartesian(xlim=c(-5,5))
    }
    
    A <- plotFun(dat=mstats, x=cols[1], y="R2", sigs=sigReal[[1]])
    B <- plotFun(dat=mstats, x=cols[2], y="R2", sigs=sigReal[[2]])
    AB <-plotFun(dat=mstats, x=cols[3], y="R2", sigs=sigReal[[3]])
    
    plot <- ggarrange(A, B, AB, ncol=3)
    annotate_figure(plot, top=text_grob("Betas after LFC>X and qval<0.05 Selection"))
    plot
  })
  
  # RANDOM
  output$betaR2_FDR <- renderPlot({
    mucovdf <- compFalseDisc()[[1]]
    plot(mucovdf$mu, mucovdf$cov, pch=".", asp=1,
         xlab="replicate mean", ylab="replicate variation",
         main="Mean ~ Variation relation between sample replicates")
  }, width=400, height=400)
  
  output$FDRdistribution <- renderPlot({
    FDRb <- compFalseDisc()[[2]]
    FDRa <- compFalseDisc()[[3]]
    FDRab <-compFalseDisc()[[4]]
    
    plot(seq(0,4,0.001), FDRab, log='x', type='l', col='green',lwd=1.5,
         xlim=c(0.001,2), ylim=c(0,0.99),
         xlab="Fold change", ylab="Probability of FDR",
         main="Distribution of False Discoveries over Fold Changes")
    lines(seq(0,4,0.001), FDRb, log='x', type='l', col='orange',lwd=1.5)
    lines(seq(0,4,0.001), FDRa, log='x', type='l', col='blue',lwd=1.5)
    abline(h=0.05, col="red", lwd=2)
    
    #print(head(joinFDRandGenes()))
    
  },width=400, height=400)
  
  output$files <- renderText({
    paste("Number of genes after Count Threshold:", length(normRatios()[,1]))
  })
  
  
  output$volcanoPlot <- renderPlot({
    mstats <- deconvolute()
    Pvals <- Pvalues_real()
    ptype <- as.character(input$plottype)
    
    X=mstats[,(as.integer(ptype)+1)]
    Y=Pvals[,(as.integer(ptype)+1)]
    
    ggplot()+
      geom_point(aes(x=X, y=Y))+
      scale_y_continuous(trans="log10")+
      xlab(colnames(mstats[as.integer(ptype)+1])) + ylab("log P-value")+
      ggtitle(paste0("Volcano Plot of ", colnames(mstats[as.integer(ptype)+1])), " Induced Genes")+
      theme(plot.title = element_text(size = 20, face = "bold"))
    
  })
  
  clicked <- reactive({
    mstats <- deconvolute()
    Pvals <- Pvalues_real()
    
    ptype <- as.character(input$plottype)
    
    df <- data.frame(cbind(mstats[,2:4], Pvals[,2:4]))
    
    ggdf <- data.frame(X=mstats[,(as.integer(ptype)+1)], Y=Pvals[,(as.integer(ptype)+1)])
    ggdf <- cbind(ggdf, df)
    brushedPoints(ggdf, input$plot_brush)
  })
  
  output$clickedPoints <- renderTable({
    clicked()[,3:8]
  }, rownames=T, digits=4)
  
  
  output$stats <- renderDT({
    sigReal <- sigReal()
    boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
    df <- deconvolute()
    df[,c(6,7,8)] <- Pvalues_real()[,c(2,3,4)]
    colnames(df)[6:8] <- colnames(Pvalues_real())[2:4]
    round(df[boolfilt,c(2,3,4,6,7,8,9)], digits=6)
    
  })
  
  HDF <- reactive({ #heatmap df function
    df <- deconvolute()
    p <- Pvalues_real()
    
    assign_genes <- function(df, B_thr=1, R_thr=0.8){
      print(colnames(df))
      SIG_genes <- ((abs(df[,2]) > B_thr & p[,2] < 0.05) |
                      (abs(df[,3]) > B_thr & p[,3] < 0.05) |
                      (abs(df[,4]) > B_thr & p[,4] < 0.05)) & df$R2>R_thr
      print(paste("Significant Genes: ", sum(SIG_genes)))
      df[,6:8] <- p[,2:4]
      sigs <- df[SIG_genes,]
      sts <- sign(sigs[,2:4])*(sigs[,6:8]<0.05)
      sts[sts==-1] <- 2 # positive reg 1, negative reg 2, no reg, 0
      clus <- rowSums(t(t(sts)*c(1,3,9))) # generate 26 unique cluster labels based on regulation
      c_order <- c(1,2,3,6,21,22,19,15,17,11,9,12,10,13,18,24,20,26,4,5,7,8,14,23,16,25)
      
      clus_inc <- which(table(clus)>30)
      cc_order <- paste0("clus\n", c_order[c_order %in% clus_inc])
      hmdf <- sigs[(clus %in% clus_inc), 2:4]
      clus.split <- factor(paste0("clus\n", clus[clus %in% clus_inc]),
                           levels=cc_order)
      col_fun = colorRamp2(c(-2,0, 2), c("blue","white", "red"))
      
      heatmap_obj <- Heatmap(as.matrix(hmdf), split=clus.split, col=col_fun,
                             cluster_row_slices = F, cluster_columns = F)
      
      
      return(heatmap_obj)
      
    }
    
    return(assign_genes(df, R_thr=input$R2Thres))
    
  })
  
  outputDF <- reactive({
    odf <- deconvolute()
    return(colnames(odf))
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("datatable", ".csv", sep = "")
    },
    content = function(file) {
      # add try catch with try:df+FDR , catch:df
      
      #sigReal <- sigReal()
      #boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
      #df <- deconvolute()
      #df[,c(6,7,8)] <- Pvalues_real()[,c(2,3,4)]
      #colnames(df)[6:8] <- colnames(Pvalues_real())[2:4]
      #outcsv <- df[boolfilt,c(2,3,4,6,7,8,9)]
      
      outcsv <- joinFDRandGenes()
      
      write.csv(outcsv, file, row.names = T)
    }
  )
  
  output$coln <- renderText({outputDF()})
  
  observeEvent(input$norm, {print("apply norm")})
  observeEvent(input$deconvolute, {print("Deconvolute")})
  #observeEvent(input$lfcThres)
  observe(firstFilter())
  observeEvent(input$H0Thres, {print(input$H0Thres)})
  observe(outputDF())
  observeEvent(input$R2Thres, {print(input$R2Thres)})
  observeEvent(input$plottype, {print(input$plottype)})
  observeEvent(input$compFDR, {print("Compute FDR")})
  observeEvent(input$show_heatmap, {
    ht1 <- HDF()
    InteractiveComplexHeatmapModal(input,output, session, ht1)
  })
}

shinyApp(ui, server)

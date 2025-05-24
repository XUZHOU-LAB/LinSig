library(DT)
library(ggplot2)
library(dplyr)
library(BiocManager)
library(circlize) # for colorRamp2 function
library(ComplexHeatmap)
library(InteractiveComplexHeatmap)
#library(RColorBrewer)
library(ggvenn)
library(patchwork)

# GSEA libraries
library(AnnotationDbi)
library("org.Mm.eg.db")
library(clusterProfiler)
library(markdown)


source("./linsig_helper_functions.R")
source("./gen_synthetic_data_helpers.R")

# for deploying paste in R console:
#library(BiocManager)
#options(repos = BiocManager::repositories())

#TODO
# colnames error message -> column name function to apply everywhere
# fix up heatmap functions repetition (improve speed?)
# remove LOWESS option (or keep?)
# Use ratios.fit function instead of custom deconvoluteFunc
# Add GSEA download results + speed up

#TODO:
# gneratate cluster funcite vervangen door gneratecluster1 - niet essentieel
# dotplot size aanpassen? - niet essentieel
# outputCSV werkend maken met en zonder FDR - nodig.
# Background Genes? - niet nodig. standaard is genoeg

# MULTIPLE REPLICATES (difficult)


###################################################
########## U S E R   I N T E R F A C E ############
###################################################

ui = fluidPage(
  sidebarPanel(
    fileInput("df", "Upload Counts, q-values",
              accept = c("text/csv",
                         "text/comma-separated-values,text/plain",
                         ".csv")),
    numericInput("reps", "Number of Replicates",
                 value = 2, min = 2, step = 1),
    checkboxInput("lowess", "LOWESS Norm", value = FALSE), # Kept here as per general controls
    hr(),
    sliderInput("lfcThres", "LogFoldChange Threshold",
                value = 0.585, min = 0, max = 3, step = 0.001),
    actionButton("deconvolute", "Deconvolute Signals"),
    hr(),
    actionButton("compFDR", "Compute False Discovery Rate"),
    hr(), # Added for visual separation before download
    downloadButton('downloadData', 'Download Data')
  ),
  
  
  mainPanel(
    tabsetPanel(type = "tabs",
                tabPanel("Data Table", 
                         dataTableOutput("stats"), 
                         plotOutput("venn")),
                tabPanel("R2~Beta Plots", plotOutput("betaR2")),
                tabPanel("FDR",
                         fluidRow(
                           splitLayout(cellWidths = c("50%", "50%"), 
                                       plotOutput("betaR2_FDR"), 
                                       plotOutput("FDRdistribution"))
                         ),
                         verbatimTextOutput("cFDR"),
                         verbatimTextOutput("pval")),
                tabPanel("Volcano Plots",
                         selectInput("plottype", "Choose a term for Volcano Plot:",
                                     choices = c("1", "2", "3")),
                         plotOutput('volcanoPlot', click = 'plot_click',
                                    brush = 'plot_brush'),
                         tableOutput('clickedPoints')),
                tabPanel("Heatmap",
                         actionButton("show_heatmap", "Generate Heatmap"),
                         htmlOutput("heatmap_output")),
                tabPanel("Enrichment Analysis",
                         numericInput("minimumGenesinClus", "Minimum Genes in a Cluster for GSEA",
                                      value = 20, min = 1, max = 1000, step = 1),
                         checkboxGroupInput("clusterChoice", "Choose Clusters for GSEA", choices = NA),
                         actionButton("actionButtonEnrich", "Start GSEA"),
                         plotOutput("enrichPlot")
                ),
                tabPanel("Advanced Parameters", # Maybe move some of them back to the main page?
                         numericInput("pseudo", "Add Pseudo Count",
                                      value = 1,
                                      min = 0.1,
                                      step = 0.5),
                         numericInput("cntThres", "Count Threshold",
                                      value = 10,
                                      min = 0,
                                      step = 0.1),
                         hr(),
                         h4("Threshold Parameters for Analysis"),
                         sliderInput("lfcThres", "LogFoldChange Threshold",
                                     value = 0.585, min = 0, max = 3, step = 0.001),
                         sliderInput("H0Thres", "Null-Hypo Test FC Threshold",
                                     value = 1.5, min = 1, max = 4, step = 0.01),
                         sliderInput("R2Thres", "R2 Threshold",
                                     value = 0.8, min = 0.1, max = 1, step = 0.01),
                         hr(),
                         h4("FDR Calculation Parameters"),
                         numericInput("RDFsize", "Random Dataset Size for FDR",
                                      value = 20000, min = 10000, step = 1),
                         sliderInput("FCFDR", "FC Selection for FDR Significance",
                                     value = 1, min = 1, max = 5, step = 0.01)
                ),
                tabPanel("Instructions",
                         downloadButton("downloadExample", "Download Example Dataset"),
                         includeMarkdown("./instructions.html")
                )
    )
  )
)


##############################################################
####################### S E R V E R ##########################
##############################################################

server = function(input,output, session){
  
  options(shiny.maxRequestSize=20*1024^2) #20MB max file size
  
  inputfile <- reactive({
    req(input$df)
    df <- read.csv(input$df$datapath, header = T, row.names = 1)
    return(df)
  })
  
  # Function that filters based on q-value and Fold Change
  firstFilter <- reactive({
    inputdf <- inputfile()
    
    DataPseudo <- MedianNorm(base::as.matrix(inputdf[,1:8]), countThres=input$cntThres, pseudo=input$pseudo)
    
    F_AvsC = log2((DataPseudo[,3] + DataPseudo[,4]) / (DataPseudo[,1] + DataPseudo[,2]))
    F_ABvsB = log2((DataPseudo[,7] + DataPseudo[,8]) / (DataPseudo[,5] + DataPseudo[,6]))
    F_ABvsA = log2((DataPseudo[,7] + DataPseudo[,8]) / (DataPseudo[,3] + DataPseudo[,4]))
    F_BvsC = log2((DataPseudo[,5] + DataPseudo[,6]) / (DataPseudo[,1] + DataPseudo[,2]))
    
    Fold <- data.frame(F_AvsC, F_ABvsB, F_ABvsA, F_BvsC)
    
    if (ncol(inputdf) != input$reps*4){ # optional q-values filtering, if no q-values provided assume q_val=0
      inputdf[,9:12] <- 0
    }
    
    qvals <- inputdf[,9:12]
    countthresholdFilter <- rownames(DataPseudo)
    qvalDF <- qvals[rownames(qvals) %in% countthresholdFilter,]
    
    # grab all genes that have a higher foldchange than lfcThres and have a significant q-value for at least one of the conditions.
    SignificantGenesIDX = rowSums(abs(Fold[,1:4]) >= input$lfcThres & qvalDF <= 0.05)>0
    
    return(SignificantGenesIDX)
  })
  

  deconvolute <- eventReactive(input$deconvolute, {
    sigGenes <- firstFilter()
    deconvoluteFunction(inputfile(), input$cntThres,
                        n_rep=input$reps, input$H0Thres,
                        pseudo=input$pseudo, lowess=input$lowess)
  })

  
  sigReal <- reactive({
    real_stats <- deconvolute()
    A <- (abs(real_stats[,2]) > log2(input$FCFDR) &
            real_stats$R2 > input$R2Thres & real_stats[,6] < 0.05)
    B <- (abs(real_stats[,3]) > log2(input$FCFDR) &
            real_stats$R2 > input$R2Thres & real_stats[,7] < 0.05)
    AB <- (abs(real_stats[,4]) > log2(input$FCFDR) &
             real_stats$R2 > input$R2Thres & real_stats[,8]< 0.05)
    
    print(paste("real A|B|AB:", sum(A|B|AB)))
    print(paste("A B AB", sum(A), sum(B), sum(AB)))
    
    return(list("A"=A, "B"=B, "AB"=AB))
  })
  
  #returns LFC distribution of deconvoluted random genes
  compFalseDisc <- eventReactive(input$compFDR, {
    #normalize countdata
    cntMat <- as.matrix(inputfile()[,1:8])
    
    # generate dataframe with mean row ~ coefficient of variation (between replicates)
    sampledMuCoV <- generate_mucov_df(cntMat, size=input$RDFsize, nrep=input$reps,
                                      countThres=input$cntThres)
    
    # draw new counts from normal distribution using mu (row mean) and CoV
    RandomDF <- t(apply(sampledMuCoV, 1, FUN=drawreps, nrep=input$reps))
    
    colnames(RandomDF) <- c("c_A","c_B", "A_A","A_B", "B_A", "B_B", "AB_A", "AB_B")

    modelStats <- deconvoluteFunction(RandomDF, input$cntThres,
                                      n_rep=input$reps, H0_threshold=1,
                                      pseudo=input$pseudo, lowess=input$lowess) # why hardcoded 1 (=0) here?
    # I guess because we want to test for any False Discovered genes (so LFC>0) and not just False Discoveries that are above e.g. FC 1.5
    
    FDRA <- modelStats[modelStats[,7] < 0.05 & modelStats$R2 > 0.8,] # make these parameters move with regular model parameters?
    FDRB <- modelStats[modelStats[,6] < 0.05 & modelStats$R2 > 0.8,] # 0.7 for adjusted R2, 0.8 for Multiple Rsquared
    FDRAB<- modelStats[modelStats[,8] < 0.05 & modelStats$R2 > 0.8,]
   
    xs <- seq(0,4,0.001) # TODO: make generic parameter
    FDRa <- sapply(xs, function(t) mean(abs(FDRA$A) > t)) # for each LFC check if above threshold 0-4
    FDRb <- sapply(xs, function(t) mean(abs(FDRB$B) > t)) # will generate a for each threshold a percentage of genes above it
    FDRab <- sapply(xs, function(t) mean(abs(FDRAB$AB) > t)) # where this threshold is 5%, its the accepted FDR value
    
    print(paste("5% FDR LFC threshold A:", NearestNeighbours(FDRa, xs, 0.05)))
    print(paste("5% FDR LFC threshold B:", NearestNeighbours(FDRb, xs, 0.05)))
    print(paste("5% FDR LFC threshold AB:", NearestNeighbours(FDRab, xs, 0.05)))
    
    return(list(sampledMuCoV, FDRa, FDRb, FDRab))
  })
  
  #returns DF with LFC, COV and FDR values
  joinFDRandGenes <- reactive({
    #after deconvolute button
    decoDF <- deconvolute()
    sigReal <- sigReal()
    boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
    
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
    colnames(decoDF)[10:12] <- paste0("FDR_", colnames(decoDF[,2:4])) # col_func
    filteredDF <- round(decoDF[boolfilt,c(1,2,3,4,6,7,8,9,10,11,12)], digits=6)
    
    return(filteredDF)
  })
  
  # VENN DIAGRAM
  output$venn <- renderPlot({
    df <- sigReal()
    # use data frame as input
    modelStats <- deconvolute() # col_func

    M <-tibble('A' = df[[1]],
               'B' = df[[2]],
               'AB'= df[[3]])
    names(M) <- colnames(modelStats)[2:4]
    # create Venn diagram and display all sets
    ggvenn(M, fill_color=c("blue","red", "purple"), fill_alpha=0.25)
  })
  
  # Beta vs R2 Plot 
  output$betaR2 <- renderPlot({
    modelStats <- deconvolute()
    sigReal <- sigReal()
    cols <- colnames(modelStats)[2:4] # Typically A, B, AB col_func
    R2 <- modelStats$R2
    
    # Generate individual plots
    A_plot  <- plotLFC_R2(dat = modelStats, x = cols[1], y = "R2", sigs = sigReal[[1]], xlab=paste("LFC", cols[1]))
    B_plot  <- plotLFC_R2(dat = modelStats, x = cols[2], y = "R2", sigs = sigReal[[2]], xlab=paste("LFC", cols[2]))
    AB_plot <- plotLFC_R2(dat = modelStats, x = cols[3], y = "R2", sigs = sigReal[[3]], xlab=paste("LFC", cols[3]))
    
    # Combine using patchwork
    combined_plot <- (A_plot | B_plot | AB_plot) +
      plot_annotation(
        title = "Betas after LFC > X and qval < 0.05 Selection",
        theme = theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold"))
      )
    
    combined_plot
  })
  
  # Mean ~ Variation Plot of Random Data
  output$betaR2_FDR <- renderPlot({
    mucovdf <- compFalseDisc()[[1]]
    plot(mucovdf$mu, mucovdf$cov, pch=".", asp=1,
         xlab="replicate mean", ylab="replicate variation",
         main="Mean ~ Variation relation between sample replicates")
  }, width=400, height=400)
  
  # FDR Cumulative Distribution Plot 
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
    
  },width=400, height=400)
  

  # Volcano Plots
  output$volcanoPlot <- renderPlot({
    modelStats <- deconvolute()
    ptype <- as.character(input$plottype)
    
    X=modelStats[,(as.integer(ptype)+1)] # lfc values
    Y=modelStats[,(as.integer(ptype)+1+4)] # pvalues
    
    ggplot()+
      geom_point(aes(x=X, y=Y))+
      scale_y_continuous(trans="log10")+
      xlab(colnames(modelStats[as.integer(ptype)+1])) + ylab("log P-value")+
      ggtitle(paste0("Volcano Plot of ", colnames(modelStats[as.integer(ptype)+1])), " Induced Genes")+
      theme(plot.title = element_text(size = 20, face = "bold"))
    
  })
  
  # Click function for Volcano Plot
  clicked <- reactive({
    modelStats <- deconvolute()
    ptype <- as.character(input$plottype)
    df <- data.frame(cbind(modelStats[,2:4], modelStats[,6:8]))
    
    X=modelStats[,(as.integer(ptype)+1)]   # columns 2,3,4 for LFC
    Y=modelStats[,(as.integer(ptype)+1+4)] # columns 6,7,8 for P-values
    
    ggdf <- data.frame(x=X,y=Y) 
    ggdf <- cbind(ggdf, df)
    brushedPoints(ggdf, input$plot_brush)
  })
  
  # Render table of selected genes in Volcano Plot
  output$clickedPoints <- renderTable({
    clicked()[,3:8]
  }, rownames=T, digits=4)
  
  # Render Table of all Genes
  output$stats <- renderDT({
    df <- getSignificantOutputTable()
    round(df, digits=4)
    df$cluster <- generateClusters1()
    df[,c(1:4,6:10)]
  })
  
  getSignificantOutputTable <- reactive({ # add these stats to the full table?
    sigReal <- sigReal()
    boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
    df <- deconvolute()
    df <- df[boolfilt,c(1,2,3,4,5,6,7,8,9)]
  })
  
  generateClusters <- reactive({
  
    minGenesInCluster = input$minimumGenesinClus

    sigs <- getSignificantOutputTable()
    
    sts <- sign(sigs[,2:4])*(sigs[,6:8]<0.05)
    sts[sts==-1] <- 2 # positive reg 1, negative reg 2, no reg, 0
    clus <- rowSums(t(t(sts)*c(1,3,9))) # generate 26 unique cluster labels based on regulation 1:LPS, 3:pH, 9:pHLPS
    cluster_order <- c(1,2,3,6,21,22,19,15,17,11,9,12,10,13,18,24,20,26,4,5,7,8,14,23,16,25,0)

    clustersToInclude <- which(table(clus)>minGenesInCluster) # parameter for minimum amount of genes in a cluster
    clustersToInclude <- names(clustersToInclude)
    
    clusterDF <- sigs[(clus %in% clustersToInclude), 2:4]
    geneClusters <- clus[clus %in% clustersToInclude]
    
    modelTerms <- colnames(sigs[,2:4]) #col_func
    clusterCode = c("T2+", "T2-", "T1+", "T1-", "T1+T3-", "T1+T2+T3-", "T2+T3-",
                    "T1-T3+", "T1-T2-T3+", "T2-T3+", "T3+", "T1+T3+", "T2+T3+",
                    "T1+T2+T3+", "T3-", "T1-T3-", "T2-T3-", "T1-T2-T3-", "T1+T2+",
                    "T1+T2-", "T1-T2+", "T1-T2-", "T1+T2-T3+", "T1+T2-T3-",
                    "T1-T2+T3+", "T1-T2+T3-","no regulation")
    clusterCode <- gsub("T1", modelTerms[2], clusterCode)
    clusterCode <- gsub("T2", modelTerms[1], clusterCode)
    clusterCode <- gsub("T3", modelTerms[3], clusterCode)
    
    clusterNames <- data.frame(clusterID = cluster_order,
                               clusterCode = clusterCode)
    geneClusters <- clusterNames$clusterCode[match(geneClusters, clusterNames$clusterID)]
    clusterDF$cluster <- geneClusters
    return(clusterDF)
  })
  

  generateClusters1 <- reactive({
    decoDF <- getSignificantOutputTable()
    sts <- sign(decoDF[,2:4])*(decoDF[,6:8]<0.05)

    sts[sts==-1] <- 2 # positive reg 1, negative reg 2, no reg, 0
    clus <- rowSums(t(t(sts)*c(1,3,9))) # generate 26 unique cluster labels based on regulation 1:LPS, 3:pH, 9:pHLPS
    
    modelTerms <- colnames(decoDF[,2:4]) # replace pH / LPS etc with model terms
    
    cluster_order <- c(1,2,3,6,21,22,19,15,17,11,9,12,10,13,18,24,20,26,4,5,7,8,14,23,16,25,0)
    clusterCode = c("T2+", "T2-", "T1+", "T1-", "T1+T3-", "T1+T2+T3-", "T2+T3-",
                    "T1-T3+", "T1-T2-T3+", "T2-T3+", "T3+", "T1+T3+", "T2+T3+",
                    "T1+T2+T3+", "T3-", "T1-T3-", "T2-T3-", "T1-T2-T3-", "T1+T2+",
                    "T1+T2-", "T1-T2+", "T1-T2-", "T1+T2-T3+", "T1+T2-T3-",
                    "T1-T2+T3+", "T1-T2+T3-","no regulation")
    clusterCode <- gsub("T1", modelTerms[2], clusterCode)
    clusterCode <- gsub("T2", modelTerms[1], clusterCode)
    clusterCode <- gsub("T3", modelTerms[3], clusterCode)
    
    clusterNames <- data.frame(clusterID = cluster_order,
                               clusterCode = clusterCode)
    geneClusters <- clusterNames$clusterCode[match(clus, clusterNames$clusterID)]

    return(geneClusters)
  })

  
  observeEvent(generateClusters(), {
    choices <- unique(generateClusters()$cluster)
    updateCheckboxGroupInput(inputId = "clusterChoice", choices = choices) 
  })
  
  #Render Heatmap Function
  HDF <- reactive({ 
    df <- deconvolute()
    p <- df[,5:8]
    
    assign_genes <- function(df, B_thr=0.585, R_thr=0.8){
      minimumGenesInClus <- 20
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
      
       
      clustersToInclude <- which(table(clus)>minimumGenesInClus) # parameter for minimum amount of genes in a cluster
      clustersToInclude <- names(clustersToInclude)
      
      print(clustersToInclude)
      
      
      
      hmdf <- sigs[(clus %in% clustersToInclude), 2:4]
      
      clusterCode = c("T2+", "T2-", "T1+", "T1-", "T1+T3-", "T1+T2+T3-", "T2+T3-",
                      "T1-T3+", "T1-T2-T3+", "T2-T3+", "T3+", "T1+T3+", "T2+T3+",
                      "T1+T2+T3+", "T3-", "T1-T3-", "T2-T3-", "T1-T2-T3-", "T1+T2+",
                      "T1+T2-", "T1-T2+", "T1-T2-", "T1+T2-T3+", "T1+T2-T3-",
                      "T1-T2+T3+", "T1-T2+T3-")
      
      modelTerms <- colnames(df)[2:4]
      clusterCode <- gsub("T1", modelTerms[2], clusterCode)
      clusterCode <- gsub("T2", modelTerms[1], clusterCode)
      clusterCode <- gsub("T3", modelTerms[3], clusterCode)
      
      clusterNames <- data.frame(clusterID = c_order,
                                 clusterCode = clusterCode)

      cc_order <- c_order[c_order %in% clustersToInclude]
      clus.split <- factor(clus[clus %in% clustersToInclude],
                           levels=cc_order)
      
      geneClusters <- clusterNames$clusterCode[match(clus.split, clusterNames$clusterID)]
      geneClusters_ord <- factor(geneClusters, level=clusterCode)
      
      col_fun = colorRamp2(c(-2,0, 2), c("blue","white", "red"))
      
      
      heatmap_obj <- Heatmap(as.matrix(hmdf), 
                             split = geneClusters_ord, 
                             col = col_fun,
                             cluster_row_slices = F,
                             cluster_columns = F,
                             show_row_dend = F,
                             heatmap_legend_param = list(title = "LFC"))
      return(heatmap_obj)
    }
    heatmap_obj <- assign_genes(df, R_thr=input$R2Thres)
    return(heatmap_obj)
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste("datatable", ".csv", sep = "")
    },
    content = function(file) {
      outcsv <- tryCatch( 
        {
          joinFDRandGenes()
        },
        error = function(e) {
          decoDF <- deconvolute()
          sigReal <- sigReal()
          boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])
          filteredDF <- round(decoDF[boolfilt,c(1,2,3,4,6,7,8,9)], digits=6)
          return(filteredDF)
          
        }
      )
      
      write.csv(outcsv, file, row.names = T)
    }
  )
  
  startGSEA <- eventReactive(input$actionButtonEnrich, {
      clusterDT <- generateClusters()
      clusterDT <- clusterDT[clusterDT$cluster %in% input$clusterChoice,]
      
      clusterDT$entrez <- mapIds(org.Mm.eg.db, keys = rownames(clusterDT),
                                 column = "ENTREZID", keytype = "SYMBOL")
      clusterDT <- na.omit(clusterDT)

      clusterList <- list()
      for (cluster in unique(clusterDT$cluster)){
        entrezInCluster <- clusterDT[clusterDT$cluster==cluster,]$entrez
        cluster <- as.character(cluster)
        clusterList[[cluster]] <- entrezInCluster
      }

      ck <- compareCluster(geneCluster = clusterList, 
                           fun = enrichGO, 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP")
      
      return(dotplot(ck,
                     label_format = 125))
  })
  
  output$enrichPlot <- renderPlot({
    startGSEA()
  })
  
  observeEvent(input$norm, {print("apply norm")})
  observeEvent(input$deconvolute, {print("Deconvolute")})
  observe(firstFilter())
  observeEvent(input$H0Thres, {print(input$H0Thres)})
  observeEvent(input$R2Thres, {print(input$R2Thres)})
  observeEvent(input$plottype, {print(input$plottype)})
  observeEvent(input$compFDR, {print("Compute FDR")})
  
  observeEvent(input$show_heatmap, {
    ht1 <- HDF()
    InteractiveComplexHeatmapWidget(input,output, session, ht1,
                                    output_id = "heatmap_output")
  })

  
  output$downloadExample <- downloadHandler(
    filename = function() {
      paste("IL6IL10_reduced_dataset", ".csv", sep = "")
    },
    content = function(file) {

      exampleData <- read.csv("./IL6IL10_reduced_df.csv",
                              header = T,
                              row.names = 1)
      
      write.csv(exampleData, file, row.names = T)
    }
  )

}

shinyApp(ui, server)


#clus#	LPS	pH	pHLPS code
#1	    0	  1 	 0    pH+
#2	    0	 -1 	 0    pH-
#3	    1	  0 	 0    LPS+
#6	   -1	  0 	 0    LPS-
#21	    1	  0 	-1    LPS+pHLPS-
#22	    1	  1 	-1    LPS+pH+pHLPS-
#19	    0	  1 	-1    pH+pHLPS-
#15	   -1	  0 	 1    LPS-pHLPS+
#17	   -1	 -1 	 1    LPS-pH-pHLPS+
#11	    0	 -1 	 1    pH-pHLPS+
#9	    0	  0 	 1    pHLPS+
#12	    1	  0 	 1    LPS+pHLPS+
#10	    0	  1 	 1    pH+pHLPS+
#13	    1	 -1 	 1    LPS+pH+pHLPS+
#18	    0	  0 	-1    pHLPS-
#24	   -1	  0 	-1    LPS-pHLPS-
#20	    0	 -1 	-1    pH-pHLPS-
#26	   -1	 -1 	-1    LPS-pH-pHLPS-
#4	    1	  1 	 0    LPS+pH+
#5	    1	 -1 	 0    LPS+pH-
#7	   -1	  1 	 0    LPS-pH+
#8	   -1	 -1 	 0    LPS-pH-
#14	    1	 -1 	 1    LPS+pH-pHLPS+
#23	    1	 -1 	-1    LPS+pH-pHLPS-
#16	   -1	  1 	 1    LPS-pH+pHLPS+
#25	   -1	  1 	-1    LPS-pH+pHLPS-
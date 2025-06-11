# for deploying paste in R console:
#library(BiocManager)
#options(repos = BiocManager::repositories())

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


#TODO: high priority
# Add GSEA download results + speed up
#   Add GSEA input data dataset? Just for verification purposes
# Add Cluster assignment in download
# FirstFilter? Optional? Give suggestion for LFC threshold?

#TODO: normal priority
# outputCSV werkend maken met en zonder FDR - nodig. (add FDR function to analyze model since its quick enough)
# Background Genes? - niet nodig. standaard is genoeg

# TODO: low priority
# remove LOWESS option (or keep?)
# MULTIPLE REPLICATES (difficult)
# dotplot size aanpassen? - niet essentieel
# Flow chart where genes are discarded and how many? first cnt>10, then lfc>1.5, then R2/pval, then FDR


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
                         downloadButton("downloadEnrichResults", "Download Enrichment Results"),
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
                                      value = 40000, min = 10000, step = 1),
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
    
    DataPseudo <- MedianNorm(base::as.matrix(inputdf[,1:8]), count_threshold=input$cntThres, pseudo=input$pseudo)
    
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
    sigGenesIDX <- firstFilter() # implement maybe as optional since we also have the FDR filter?
    # TODO: what to do with firstFilter? now its not being used...
    deconvoluteFunction(inputfile(), count_threshold=input$cntThres,
                        n_rep=input$reps, H0_threshold=input$H0Thres,
                        pseudo=input$pseudo, lowess=input$lowess,
                        beta_threshold=input$FCFDR, r2_threshold=input$R2Thres,
                        sig_index=sigGenesIDX)
  })
  
  
  #returns LFC distribution of deconvoluted random genes
  compFalseDisc <- reactive({
    #normalize countdata
    cntMat <- as.matrix(inputfile()[,1:8])
    
    # generate dataframe with mean row ~ coefficient of variation (between replicates)
    sampledMuCoV <- generate_mucov_df(cntMat, size=input$RDFsize, nrep=input$reps,
                                      count_threshold=input$cntThres)
    
    # draw new counts from normal distribution using mu (row mean) and CoV
    RandomDF <- t(apply(sampledMuCoV, 1, FUN=drawreps, nrep=input$reps))
    
    colnames(RandomDF) <- c("c_A","c_B", "A_A","A_B", "B_A", "B_B", "AB_A", "AB_B")

    modelStats <- deconvoluteFunction(RandomDF, input$cntThres,
                                      n_rep=input$reps, H0_threshold=1, # hardcoded here because we want to test for any False Discovered genes (so LFC>0) and not just False Discoveries that are above e.g. FC 1.5
                                      pseudo=input$pseudo, lowess=input$lowess,
                                      beta_threshold=1, # also set to LFC=0 because I want all genes in the FDR dataset in order to make a good distribution of the random genes.
                                      r2_threshold=input$R2Thres) #
   
    FDRA <- modelStats[modelStats$isSignificantA==T,]
    FDRB <- modelStats[modelStats$isSignificantB==T,]
    FDRAB <- modelStats[modelStats$isSignificantAB==T,]
    
    FDRa <- sapply(params$lfc_thresholds, function(t) mean(abs(FDRA$A) > t)) # for each LFC check if above threshold 0-4
    FDRb <- sapply(params$lfc_thresholds, function(t) mean(abs(FDRB$B) > t)) # will generate a for each threshold a percentage of genes above it
    FDRab <- sapply(params$lfc_thresholds, function(t) mean(abs(FDRAB$AB) > t)) # where this threshold is 5%, its the accepted FDR value
    
    print(paste("5% FDR LFC threshold A:", NearestNeighbours(FDRa, params$lfc_thresholds, 0.05)))
    print(paste("5% FDR LFC threshold B:", NearestNeighbours(FDRb, params$lfc_thresholds, 0.05)))
    print(paste("5% FDR LFC threshold AB:", NearestNeighbours(FDRab, params$lfc_thresholds, 0.05)))
    
    return(list(sampledMuCoV, FDRa, FDRb, FDRab))
  })
  
  #returns DF with LFC, COV and FDR values
  joinFDRandGenes <- reactive({
    #after deconvolute button
    decoDF <- deconvolute()
    boolfilt <- decoDF$isSignificantA | decoDF$isSignificantB | decoDF$isSignificantAB
    
    fdr_df <- compFalseDisc()
    
    #after FDR computation
    FDRa_dist <- fdr_df[[2]]
    FDRb_dist <- fdr_df[[3]]
    FDRab_dist<- fdr_df[[4]]
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
    df <- deconvolute()

    M <-tibble('A' = df$isSignificantA,
               'B' = df$isSignificantB,
               'AB'= df$isSignificantAB)
    names(M) <- colnames(df)[2:4]
    # create Venn diagram and display all sets
    ggvenn(M, fill_color=c("blue","red", "purple"), fill_alpha=0.25)
  })
  
  # Beta vs R2 Plot 
  output$betaR2 <- renderPlot({
    modelStats <- deconvolute()
    cols <- colnames(modelStats)[2:4] # Typically A, B, AB col_func
    R2 <- modelStats$R2
    
    # Generate individual plots
    A_plot  <- plotLFC_R2(dat = modelStats, x = cols[1], y = "R2", sigs = modelStats$isSignificantA, xlab=paste("LFC", cols[1]))
    B_plot  <- plotLFC_R2(dat = modelStats, x = cols[2], y = "R2", sigs = modelStats$isSignificantB, xlab=paste("LFC", cols[2]))
    AB_plot <- plotLFC_R2(dat = modelStats, x = cols[3], y = "R2", sigs = modelStats$isSignificantAB, xlab=paste("LFC", cols[3]))
    
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
    
    plot(params$lfc_thresholds, FDRab, log='x', type='l', col='green',lwd=1.5,
         xlim=c(0.001,2), ylim=c(0,0.99),
         xlab="Fold change", ylab="Probability of FDR",
         main="Distribution of False Discoveries over Fold Changes")
    lines(params$lfc_thresholds, FDRb, log='x', type='l', col='orange',lwd=1.5)
    lines(params$lfc_thresholds, FDRa, log='x', type='l', col='blue',lwd=1.5)
    abline(h=0.05, col="red", lwd=2, lty=2)
    
  },width=400, height=400)
  

  # Volcano Plots
  output$volcanoPlot <- renderPlot({
    modelStats <- deconvolute()
    ptype <- as.character(input$plottype)

    ggplot()+
      geom_point(aes(x=modelStats[,(as.integer(ptype)+1)], # lfc values
                     y=modelStats[,(as.integer(ptype)+1+4)] # pvalues
                     ))+
      scale_y_continuous(trans="log10")+
      xlab(colnames(modelStats[as.integer(ptype)+1])) + ylab("log P-value")+ # do something
      ggtitle(paste0("Volcano Plot of ", colnames(modelStats[as.integer(ptype)+1])), " Induced Genes")+
      theme(plot.title = element_text(size = 20, face = "bold"))
    
  })
  
  # Click function for Volcano Plot
  clicked <- reactive({
    modelStats <- deconvolute()
    ptype <- as.character(input$plottype)
    df <- data.frame(cbind(modelStats[,2:4], modelStats[,6:8]))

    ggdf <- data.frame(x=modelStats[,(as.integer(ptype)+1)],   # columns 2,3,4 for LFC
                       y=modelStats[,(as.integer(ptype)+1+4)] # columns 6,7,8 for P-values
                       ) 
    ggdf <- cbind(ggdf, df)
    brushedPoints(ggdf, input$plot_brush, xvar = "x", yvar = "y")
  })
  
  # Render table of selected genes in Volcano Plot
  output$clickedPoints <- renderTable({
    clicked()[,3:8]
  }, rownames=T, digits=4)
  
  
  getSignificantOutputTable <- reactive({ # add these stats to the full table?
    df <- deconvolute()
    boolfilt <- df$isSignificantA | df$isSignificantB | df$isSignificantAB
    df <- df[boolfilt,c(1:9)]
  })
  
  # Render Table of all Genes
  output$stats <- renderDT({
    df <- getSignificantOutputTable()
    df <- round(df, digits=4)
    df$cluster <- generateClusters1()
    df[,c(2:4,6:10)] # 
  })
  
  
  HDF <- reactive({ 
    df <- deconvolute()
    p <- df[, 5:8]
    
    heatmap_obj <- assign_genes(df, R_thr = input$R2Thres)
    return(heatmap_obj)
  })
  
  generateClusters <- reactive({
    sigs <- getSignificantOutputTable()
    minGenes <- input$minimumGenesinClus

    clus <- get_cluster_labels(sigs)
    
    included <- names(which(table(clus) > minGenes))
    clusterDF <- sigs[clus %in% included, 2:4]
    geneClusters <- clus[clus %in% included]
    
    modelTerms <- colnames(sigs)[2:4]
    clusterNames <- get_cluster_code_mapping(modelTerms)
    clusterDF$cluster <- clusterNames$clusterCode[match(geneClusters, clusterNames$clusterID)]
    
    return(clusterDF)
  })
  
  generateClusters1 <- reactive({
    df <- getSignificantOutputTable()

    clus <- get_cluster_labels(df)
    
    modelTerms <- colnames(df)[2:4]
    clusterNames <- get_cluster_code_mapping(modelTerms)
    geneClusters <- clusterNames$clusterCode[match(clus, clusterNames$clusterID)]
    
    return(geneClusters)
  })
  
  observeEvent(generateClusters(), {
    choices <- unique(generateClusters()$cluster)
    updateCheckboxGroupInput(inputId = "clusterChoice", choices = choices) 
  })
  
  output$downloadData <- downloadHandler(
    filename = function() {
      paste0("model_stats_", Sys.Date(), ".csv")
    },
    content = function(file) {
      tryCatch({
        df <- joinFDRandGenes()
        clusters <- generateClusters1()
        
        # Check row alignment
        if (nrow(df) == length(clusters)) {
          df$cluster <- clusters
        } else {
          warning("Row count mismatch between data and clusters")
          df$cluster <- NA
        }
        
        write.csv(df, file, row.names = TRUE)
      }, error = function(e) {
        # Optional fallback if needed
        message("Download error: ", e$message)
        showNotification("Error while preparing download. Please check inputs.", type = "error")
      })
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

      ck <- compareCluster(geneCluster = clusterList, # this function is slow (100s)
                           fun = enrichGO, 
                           OrgDb = org.Mm.eg.db, 
                           ont = "BP")
      
      enrichResult(ck)  # Store the result
      
      return(dotplot(ck,
                     label_format = 125))
  })
  enrichResult <- reactiveVal(NULL)
  
  
  output$enrichPlot <- renderPlot({
    startGSEA()
  })
  
  output$downloadEnrichResults <- downloadHandler(
    filename = function() {
      paste0("GO_enrichment_results_", Sys.Date(), ".csv")
    },
    content = function(file) {
      res <- enrichResult()
      if (is.null(res)) {
        write.csv(data.frame(Message = "No enrichment results available."), file, row.names = FALSE)
      } else {
        write.csv(as.data.frame(res), file, row.names = FALSE)
      }
    }
  )
  
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
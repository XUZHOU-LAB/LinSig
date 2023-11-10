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
    fileInput("rdf", "Upload Random Counts",
              accept=c("text/csv", ".csv")),

    checkboxInput("lowess", "LOWESS Norm", value=FALSE),

    numericInput("pseudo", "Add Pseudo Count",value=1, min=0.1, step=0.5),
    numericInput("cntThres", "Count Threshold", value=10, min=0, step=0.1),
    actionButton("normalise","Compute Normalized Ratios"),
    hr(),
    sliderInput("lfcThres", "LFC Threshold (LFC>X & qval<0.05)", value=0.585, min=0, max=3, step=0.001),
    actionButton("deconvolute", "Deconvolute Signals"),
    hr(),
    sliderInput("H0Thres", "Hypothesis Test FC Threshold", value=1.5, min=1, max=4, step=0.01),
    sliderInput("R2Thres", "R2 Threshold", value=0.8, min=0.1, max=1, step=0.01),
    actionButton("compFDR", "Compute False Discovery Rate"),
    sliderInput("FCFDR", "FC Selection for FDR", value=1, min=1, max=5, step=0.01),
    downloadButton('downloadData', 'Download Data')
    ),


    mainPanel(
      tabsetPanel(type="tabs",
                  tabPanel("Data Table", dataTableOutput("stats"), textOutput("files"), plotOutput("venn")),
                  tabPanel("R2~Beta Plots", plotOutput("betaR2")),
                  tabPanel("FDR", plotOutput("betaR2_FDR"),
                           verbatimTextOutput("cFDR"),
                           verbatimTextOutput("pval")),
                  tabPanel("Volcano Plots",
                           selectInput("plottype", "Choose a plot:",
                                       choices = c("1", "2", "3")),
                           plotOutput('volcanoPlot', click='plot_click', brush='plot_brush'),
                           tableOutput('clickedPoints')),
                  tabPanel("Heatmap",
                           actionButton("show_heatmap", "Generate Heatmap"),
                           htmlOutput("heatmap_output"))
                  )
  )
)

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
# Function for deconvoluting signals, used for Random data
deconvolutef <- function(normRatios, firstFilter){
  Ratios <- normRatios[firstFilter,]
  strucvec <- c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,1)
  strucmat <- matrix(strucvec, ncol=3, byrow=T)

  A = length(Ratios[,1]) # of gene
  B = matrix(0, nrow=A,ncol=4);

  X <- strucmat

  colnames(X) <- c("A", "B", "AB")


  rsquare <- vector()
  ngen <- length(Ratios[,1])
  covB_l <- NULL
  covB <- data.frame(cov_int=double(ngen), cov_x1=double(ngen), cov_x2=double(ngen), cov_x3=double(ngen))

  #add progress bar
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
  modelStats$R2 <- rsquare
  return(modelStats)
}
#compute P value with threshold
Pval <- function(Betas, covB, Threshold){
  ttest_stat <- (abs(Betas) - log2(Threshold)) / sqrt(covB)
  ttest_stat <- data.frame(ttest_stat)
  Pvalue = 1 -apply(ttest_stat, 2, pt, df=4)
  return(Pvalue)
}

server = function(input,output, session){

  options(shiny.maxRequestSize=20*1024^2) #20MB max file size
  inputfile <- reactive({
    req(input$df)
    df <- read.csv(input$df$datapath, header=T,row.names=1)
    return(df)
  })


  normRatios <- eventReactive(input$normalise, {
    return(normalize(inputfile(), cntThres=input$cntThres, pseudo=input$pseudo, lowess=input$lowess))
    })
  # Function that filters based on q-value and Fold Change
  firstFilter <- eventReactive(input$deconvolute, {
    cntMat <- as.matrix(inputfile()[,1:8])

    DataPseudo <- MedianNorm(cntMat, CountThres=input$cntThres, pseudo=input$pseudo)
    F_AvsC = log2((DataPseudo[,3] + DataPseudo[,4])/(DataPseudo[,1]+ DataPseudo[,2]))
    F_ABvsB = log2((DataPseudo[,7] + DataPseudo[,8])/(DataPseudo[,5]+ DataPseudo[,6]))
    F_ABvsA = log2((DataPseudo[,7] + DataPseudo[,8])/(DataPseudo[,3]+ DataPseudo[,4]))
    F_BvsC = log2((DataPseudo[,5] + DataPseudo[,6])/(DataPseudo[,1]+ DataPseudo[,2]))

    Fold <- data.frame(F_AvsC, F_ABvsB, F_ABvsA, F_BvsC)

    qvals <- inputfile()[,9:12]
    countthresholdFilter <- rownames(DataPseudo)
    qvalDF <- qvals[rownames(qvals) %in% countthresholdFilter,]

    SigIDX = rowSums(abs(Fold[,1:4]) >= input$lfcThres & qvalDF <= 0.05)>0

    return(SigIDX)
  })
  #deconvolute for real data
  deconvolute <- eventReactive(input$deconvolute, {
    Ratios <- normRatios()[firstFilter(),]
    strucvec <- c(0,1,0,0,1,1,1,0,0,1,0,1,1,1,1,0,1,0,0,1,1,1,0,0,1,0,1,1,1,1)
    strucmat <- matrix(strucvec, ncol=3, byrow=T)

    A = length(Ratios[,1]) # of gene
    B = matrix(0, nrow=A,ncol=4);

    X <- strucmat
    cc <- strsplit(colnames(inputfile())[1:8], "_")
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

    #CompThreshold <- 1.5
    #Pvalue <- Pval(B, covB, CompThreshold)

    modelStats <- as.data.frame(cbind(B, covB))

    colnames(modelStats) <- c("int", cols[2], cols[1], cols[3], "p_int",
      paste0("cov_", cols[2]), paste0("cov_", cols[1]), paste0("cov_", cols[3]))
    modelStats$R2 <- rsquare
    modelStats <- round(modelStats, digits=7)
    print(dim(modelStats))
    return(modelStats)
  })



  FDR <- eventReactive(input$compFDR, {
    #If we generate random gene counts, how many of those will pass through our method and thresholds
    # FDR = fake genes passed through / total fake genes generated
    req(input$rdf)
    rdf <- read.csv(input$rdf$datapath, row.names=1, header=T)

    normRatios <- normalize(rdf, cntThres=input$cntThres, pseudo=input$pseudo, lowess=input$lowess)

    cntMat <- as.matrix(rdf[,1:8])
    DataPseudo <- MedianNorm(cntMat, CountThres=input$cntThres, pseudo=input$pseudo)
    F_AvsC = log2((DataPseudo[,3] + DataPseudo[,4])/(DataPseudo[,1]+ DataPseudo[,2]))
    F_ABvsB = log2((DataPseudo[,7] + DataPseudo[,8])/(DataPseudo[,5]+ DataPseudo[,6]))
    F_ABvsA = log2((DataPseudo[,7] + DataPseudo[,8])/(DataPseudo[,3]+ DataPseudo[,4]))
    F_BvsC = log2((DataPseudo[,5] + DataPseudo[,6])/(DataPseudo[,1]+ DataPseudo[,2]))

    Fold <- data.frame(F_AvsC, F_ABvsB, F_ABvsA, F_BvsC)

    lfct <- input$lfcThres
    SigIDX = rowSums(abs(Fold[,1:4]) >= lfct)>0
    print(sum(SigIDX))

    r_stats <- deconvolutef(normRatios, SigIDX)
    cc <- strsplit(colnames(rdf)[1:8], "_") #Extract column names from df
    cols <- unique(unlist(cc)[2*(1:length(cc))-1])[-1] #Extract col names
    colnames(r_stats)[1:4] <-  c("int", cols[2], cols[1], cols[3])
    print(colnames(r_stats))
    print(paste(input$H0Thres, "Pvalue input"))
    pvalues <- Pval(r_stats[,1:4], r_stats[,5:8], input$H0Thres)
    colnames(pvalues) <- c("p_int", "p_A", "p_B", "p_A+B")
    print(sum(rowSums(pvalues<=0.05)>0))

    return(cbind(r_stats, pvalues))
  })

  Pvalues_real <- reactive({
    real_stats <- deconvolute()
    p <- Pval(real_stats[,1:4], real_stats[,5:8], input$H0Thres)
    print(paste("real", sum(rowSums(p<=0.05)>0)))
    print(sum(rowSums(p>0.5)>0))
    colnames(p) <- paste0("p_", colnames(p))
    print(paste("Pvalues colnames:",colnames(p)))
    print(head(p))
    return(p)

  })

  Pvalues <- reactive({
    rand_stats <- FDR()
    p_r <- Pval(rand_stats[,1:4], rand_stats[,5:8], input$H0Thres)
    print(paste("rand", sum(rowSums(p_r<0.05)>0)))
    ps<- p_r

    return(p_r)
  })

  sigRand <- reactive({
    Pvalues <- Pvalues()
    rand_stats <- FDR()
    rr <- rand_stats$R2

    FCFDR <- log2(input$FCFDR)

    A_r<-(abs(rand_stats[,3])>FCFDR & rr>input$R2Thres & Pvalues[,3]<0.05)
    B_r<-(abs(rand_stats[,2])>FCFDR & rr>input$R2Thres & Pvalues[,2]<0.05)
    AB_r<-(abs(rand_stats[,4])>FCFDR & rr>input$R2Thres & Pvalues[,4]<0.05)

    return(list(A_r, B_r, AB_r))
  })

  sigReal <- reactive({
    Pvalues <- Pvalues_real()
    real_stats <- deconvolute()
    r2 <- real_stats$R2
    FCFDR <- log2(input$FCFDR)
    print(colnames(real_stats)[2:4])
    A <- (abs(real_stats[,2])>FCFDR & r2>input$R2Thres & Pvalues[,2]<0.05)
    B <- (abs(real_stats[,3])>FCFDR & r2>input$R2Thres & Pvalues[,3]<0.05)
    AB <- (abs(real_stats[,4])>FCFDR & r2>input$R2Thres & Pvalues[,4]<0.05)
    print(paste("real A|B|AB:", sum(A|B|AB)))
    print(paste("A B AB", sum(A), sum(B), sum(AB)))

    return(list(A, B, AB))
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

  output$betaR2_FDR <- renderPlot({
    rstats <- FDR()
    A_r <- ggplot(rstats, aes(A, R2)) +
      geom_point(size=0.5) +
      geom_point(data=rstats[sigRand()[[1]],], color="red", size=0.7) +
      coord_cartesian(xlim =c(-4, 4))
    B_r <- ggplot(rstats, aes(B, R2)) +
      geom_point(size=0.5) +
      geom_point(data=rstats[sigRand()[[2]],], color="red", size=0.7) +
      coord_cartesian(xlim =c(-4, 4))
    AB_r <- ggplot(rstats, aes(AB, R2)) +
      geom_point(size=0.5) +
      geom_point(data=rstats[sigRand()[[3]],], color="red", size=0.7) +
      coord_cartesian(xlim =c(-4, 4))


    plotr <- ggarrange(A_r, B_r, AB_r, ncol=3)
    annotate_figure(plotr, top=text_grob("Betas after LFC>X and qval<0.05 Selection"))
    plotr
  })

  output$files <- renderText({
    paste("Number of genes after Count Threshold:", length(normRatios()[,1]))
  })

  output$pval <- renderPrint({
    paste("N Genes after initial LFC threshold", nrow(FDR()))
  })

  output$cFDR <- renderPrint({
    cat(paste( paste("FDR: ", "A:", sum(sigRand()[[1]])/length(sigRand()[[1]])),
    paste("FDR: ", "B:", sum(sigRand()[[2]])/length(sigRand()[[2]])),
    paste("FDR: ","AB:", sum(sigRand()[[3]])/length(sigRand()[[3]])), sep="\n"))
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

      heatmap_obj <- Heatmap(as.matrix(hmdf), split=clus.split, col=col_fun, cluster_row_slices = F, cluster_columns = F)


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

      sigReal <- sigReal()
      boolfilt <- (sigReal[[1]] | sigReal[[2]] | sigReal[[3]])

      df <- deconvolute()
      df[,c(6,7,8)] <- Pvalues_real()[,c(2,3,4)]
      colnames(df)[6:8] <- colnames(Pvalues_real())[2:4]
      outcsv <- df[boolfilt,c(2,3,4,6,7,8,9)]

      write.csv(outcsv, file, row.names = T)
    }
  )

  output$coln <- renderText({outputDF()})

  observeEvent(input$norm, {print("apply norm")})
  observeEvent(input$deconvolute, {print("Deconvolute")})
  #observeEvent(input$lfcThres)
  observe(firstFilter())
  observeEvent(input$H0Thres, {print(input$H0Thres)})
  observe(FDR())
  observe(sigRand())
  observe(outputDF())
  observeEvent(input$R2Thres, {print(input$R2Thres)})
  observeEvent(input$plottype, {print(input$plottype)})

  observeEvent(input$show_heatmap, {
    ht1 <- HDF()
    InteractiveComplexHeatmapModal(input,output, session, ht1)
  })
}

shinyApp(ui, server)

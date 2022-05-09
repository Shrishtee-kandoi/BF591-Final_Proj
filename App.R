#BF591_Final_Project
#Shrishtee kandoi

#import libraries
library(tidyverse)
library(shiny)
library(ggplot2)
library(DT)
library(bslib)
library(RColorBrewer)
library(colourpicker)

#front-end
ui <- fluidPage(theme = shinythemes::shinytheme("cerulean"),
                titlePanel("BF591 Final project"),
                p("This website uses the Huntington's Disease dataset from (Labadorf et al., 2016) titled 'mRNA-Seq Expression profiling of human post-mortem BA9 brain tissue for Huntington's Disease and neurologically normal individuals'"),
                p("~ Shrishtee kandoi"),
                tabsetPanel(
                  #---------------------------------------------------------Samples------------------------------------------------------------
                  tabPanel("Samples", 
                           h3("Sample data: Summary, Table and Plots"),
                           p("Please use this tab to explore the 'sample_metadata.csv' file"),
                           p("1) Summary of Samples: a summary table that includes columns including data type, their mean, and Std. Dev "),
                           p("2) Table: a table of the sample information), "),
                           p("3) Plots of Samples: a bar_plot with count vs phenotype: Diagnosis, Age of Death and RIN (RNA integrity number)"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput(inputId = 'meta', label = 'Load sample metadata in csv format', accept = ".csv", placeholder = "sample_metadata.csv")
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel('Summary', tableOutput(outputId = 'summ')),
                                 tabPanel('Table', dataTableOutput(outputId = "metadata")),
                                 tabPanel('Plots', plotOutput(outputId = 'Diag'), plotOutput(outputId = 'aod'), plotOutput(outputId = 'rin'))
                               ))
                           )),
                  #---------------------------------------------------------Counts------------------------------------------------------------
                  tabPanel("Counts",
                           h3("Counts Matrix Exploration"),
                           p("Please use this tab to look at the count matrix to understand the count structure and aid in gene filtering. ",
                             "It loads a matrix of count information, as well as user input for percent variance and number of allowed non-zero samples,",
                             " and gives the following information in separate tabs:"),
                           p("1) Summary: a summary table that shows the effect of the user's filtering input, "),
                           p("2) Diagnostic plots: diagnostic scatter plots of median count vs variance and median count vs non-zeros (where genes that pass filters are shown in black), "),
                           p("3) Heatmap: a clustered heatmap of counts that remain after filtering, "),
                           p("4) PCA: PCA projections based on user input for which principle components to examine."),
                           br(),
                           p("Note: Please allow a few seconds for the 'norm_counts.csv' to load"),
                           
                           sidebarLayout(
                             sidebarPanel(
                               fileInput(inputId = 'countsID', label = 'Load normalized counts matrix CSV', accept = ".csv", placeholder = "norm_counts.csv"),
                               sliderInput(inputId = 'varID', label = 'Select the minimum percentile variance of genes', min = 0, max = 100, value = 80),
                               sliderInput(inputId = 'nonzero', label = 'Select the minimum number of non-zero samples', min = 0, max = 69, value = 5),
                               submitButton(text='Submit', width = '100%')
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel('Summary', tableOutput(outputId = 'normID')),
                                 tabPanel('Filter Plots', p('This can take a few seconds to load....Please be patient'), plotOutput(outputId = 'medvar'), plotOutput(outputId = 'medzero')),
                                 tabPanel('Heatmap', p('Loading.......'), plotOutput(outputId = "hmap")),
                                 tabPanel('PCA', selectInput(inputId = "comp1", label="Select X-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10")),
                                          selectInput(inputId = "comp2", label="Select Y-axis", choices = c("PC1", "PC2", "PC3", "PC4", "PC5", "PC5", "PC6", "PC7", "PC8", "PC9", "PC10"), selected = "PC2"),
                                          plotOutput(outputId = "PCAplot"))
                               ))
                           )),
                  #---------------------------------------------------------DE------------------------------------------------------------
                  tabPanel("DE",
                           h3("Differential Expression Results"),
                           p("Upload and explore the results of differential expression analysis"),
                           p("It loads a matrix of deseq results and user input for a desired p-value",
                             " and gives the following information in separate tabs: "),
                           p("1) Summary: a table with the differential expression results, "),
                           p("2) Volcano Plot: a volcano plot of log-fold change values"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput(inputId = 'DE_ID', label = "Load differential expression results", accept = ".csv", placeholder = "deseq_diff_exp_res.csv"),
                               HTML(paste(rep("<p> <strong> Note: </strong> This server allows you to view the table format of the data.csv files alongwith their Volcano plots </p>"), collapse = "")),
                               br(),
                               radioButtons(inputId = 'xaxis', label = 'Choose the column for the x-axis', choices = c('baseMean',
                                                                                                                       'log2FoldChange', 
                                                                                                                       'lfcSE',
                                                                                                                       'stat',
                                                                                                                       'pvalue',
                                                                                                                       'padj'), 
                                            selected = 'log2FoldChange'),
                               radioButtons(inputId = 'yaxis', label = 'Choose the column for the x-axis', choices = c('baseMean',
                                                                                                                       'log2FoldChange', 
                                                                                                                       'lfcSE',
                                                                                                                       'stat',
                                                                                                                       'pvalue',
                                                                                                                       'padj'), 
                                            selected = 'padj'),
                               colourpicker::colourInput(inputId = 'basecol', label='Base point color', value='#73C6B6'),
                               colourpicker::colourInput(inputId = 'highcol', label='Highlight point color', value='#D0ECE7'),
                               sliderInput(inputId = 'padjmag', label = 'Adjust magnitude of p adjusted subset coloring:', min = -50, max = 0, value = -10,),
                               submitButton("Run", icon("refresh"), width = '100%')
                             ),
                             mainPanel(
                               tabsetPanel(
                                 tabPanel('Table', dataTableOutput("difftable")),
                                 tabPanel('Volcano Plot', plotOutput("volcano"))
                               )
                             )
                           )),
                  
                  #---------------------------------------------------------IGE------------------------------------------------------------
                  tabPanel("IGE",
                           h3("Individual Gene Expresion Visualization"),
                           p("On this page you can visualize data about each gene individually"),
                           sidebarLayout(
                             sidebarPanel(
                               fileInput(inputId = 'IGE_countsID', label = 'Load normalized counts matrix CSV',accept = ".csv", placeholder = "norm_counts.csv"),
                               fileInput(inputId = 'IGE_metaID', label = 'Load sample information matrix CSV',accept = ".csv", placeholder = "sample_metadata.csv"),
                               selectInput("metachoice", choices = c("Diagnosis"="Diagnosis", "Diag"="Diag","Age of death"="Age_of_death", "RIN"="RIN", "Sequence Reads"="Seq_reads"),
                                           label = "Select metadata category", selected = "Age_of_death"),
                               #insert gene search box here
                               textInput("gene", label = "Enter Ensembl ID of gene (For ex: ENSG00000270773.1)", placeholder = "ENSG00000270773.1"),
                               #selectInput("plotType", label = "Choose what type of plot to make", choices = c("Bar", "Box", "Violin", "Beeswarm")),
                               submitButton(text='Plot', icon = icon("bar-chart-o"), width = '100%')
                             ),
                             mainPanel(
                               plotOutput("IGE_plot")
                             )
                           ))
                )
)

server <- function(input, output, session){
  options(shiny.maxRequestSize=30*1024^2)  #increase file upload limit
  #function to take metadata input file for first tab
  meta_load <- reactive({
    if (!is.null(input$meta)){
      meta <- read_csv(input$meta$datapath)
      return(meta)}
    else{stop('File is not a csv. Please upload the correct file format.')}
  })
  #function to take norm counts input file for counts data exploration
  counts_load <- reactive({
    if (!is.null(input$countsID)){
      counts <- read_csv(input$countsID$datapath)
      return(counts)}
    else{stop('File is not a csv. Please upload the correct file format.')}
  })
  #function to load Differential Expression .csv file
  de_load <- reactive({
    if (!is.null(input$DE_ID)){
      DE_ID <- read_csv(input$DE_ID$datapath)
      return(DE_ID)}
    else{return(NULL)}
  })
  #function to take inputs for fourth tab-
  IGE_meta <- reactive({
    if (!is.null(input$IGE_metaID)){
      meta <- read_csv(input$IGE_metaID$datapath)
      return(meta)}
    else{return(NULL)}
  })
  load_IGE_counts <- reactive({
    if (!is.null(input$IGE_countsID)){
      counts <- read_csv(input$IGE_countsID$datapath)
      return(counts)}
    else{return(NULL)}
  })
  #function to produce table
  summ_table <- function(meta_tib){
    if (!is.null(input$meta)){
      summ_tib <- tibble("Columns" = colnames(meta_tib), "Type" = sapply(meta_tib, class),
                         "Mean" = sapply(meta_tib, mean, na.rm = TRUE), "SD(+/-)" = sapply(meta_tib, sd, na.rm = TRUE))
      return(summ_tib)}
    else{return(NULL)}
  }
  
  #function to produce Diagnosis plot
  plot_Diag <- function(meta_tib){
    if (!is.null(input$meta)){
      plot <- ggplot(meta_tib, aes(Diagnosis))+
        geom_bar(bins = 10, color = "black", fill = "lavender")+
        labs(title = 'Barplot of phenotype: Diagnosis')+
        xlab('Diagnosis')+
        ylab('Count')+
        theme_bw()
      return(plot)}
    else{return(NULL)}
  }
  
  #function to produce AOD plot
  plot_aod <- function(meta_tib){
    if (!is.null(input$meta)){
      plot <- ggplot(meta_tib, aes(Age_of_death))+
        geom_bar(bins = 10, color = "black", fill = "azure")+
        labs(title = 'Barplot of phenotype: Age of Death')+
        xlab('Age of Death')+
        ylab('Count')+
        theme_bw()
      return(plot)}
    else{return(NULL)}
  }
  #function to produce RIN plot
  plot_rin <- function(meta_tib){
    if (!is.null(input$meta)){
      plot <- ggplot(meta_tib, aes(RIN))+
        geom_bar(bins = 10, color = "black", fill = "plum4")+
        labs(title = 'Barplot of phenotype: RNA Integrity Number')+
        xlab('RIN')+
        ylab('Count')+
        theme_bw()
      return(plot)}
    else{return(NULL)}
  }
  
  #function to produce summary table
  counts_sum <- function(counts_tib, perc_var, nz_genes){
    if (!is.null(input$countsID)){
      tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
      tot_genes <- nrow(counts_tib)
      counts_tib <- counts_tib %>% mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var))  #calculate variance and then percentile
      perc_val <- quantile(counts_tib$variance, probs = perc_var/100)   #calculate percentile
      counts_tib <- filter(counts_tib, variance >= perc_val)  #filter by percentile
      counts_tib <- na_if(counts_tib, 0)    #make zeroes NA's
      counts_tib$non_zero <- tot_samp-rowSums(is.na(counts_tib))  #calculate no. of genes with enough non zero samples
      counts_tib <- filter(counts_tib, non_zero >= nz_genes)  #filter by non-zero samples
      filt_genes <- nrow(counts_tib)    #calculate the number and % of genes passing the filters
      perc_pass_genes <- filt_genes/tot_genes*100
      fail_genes <- tot_genes-filt_genes
      perc_fail_genes <- fail_genes/tot_genes*100
      #produce the summary tibble
      summ_tib <- tibble('Measure' = c('Number of Samples', 'Number of Genes', 'No. of genes passing filters', "Percentage (%) passing filter", 'No. of genes not passing filters', 'Percentafe (%) not passing filter'),
                         'Value' = c(tot_samp, tot_genes, filt_genes, perc_pass_genes, fail_genes, perc_fail_genes))
      return(summ_tib)}
    else{return(NULL)}
  }
  #function to produce median vs variance plot
  med_vs_var <- function(counts_tib, perc_var){
    if (!is.null(input$countsID)){
      #make a plot tibble
      plot_tib <- counts_tib%>%
        mutate(Median = apply(counts_tib[-1], MARGIN = 1, FUN = median), 
               Variance = apply(counts_tib[-1], MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$Variance, probs = perc_var/100)   #calculate percentile
      plot_tib <- plot_tib %>% mutate(thresh = case_when(Variance >= perc_val ~ "TRUE", TRUE ~ "FALSE"))
      #scatter plot
      cols <- c("FALSE" = "blue", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, Variance))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        labs(title = 'Plot of Median vs Variance.', subtitle = "Filtered genes are marked with blue color")+
        scale_y_log10()+
        scale_x_log10()+
        theme_bw()+
        theme(legend.position = 'bottom')
      return(scatter)}
    else{return(NULL)}
  }
  #function to produce median vs non-zero samples plot
  med_vs_nz <- function(counts_tib, nz_genes){
    if (!is.null(input$countsID)){
      tot_samp <- ncol(counts_tib)-1  #store original number of samples and genes
      #make a plot tibble
      plot_tib <- counts_tib %>%   
        mutate(Median = apply(counts_tib[-1], MARGIN = 1, FUN = median)) %>% na_if(0)  #calc median, convert 0 to NA
      plot_tib$no_zeros <- rowSums(is.na(plot_tib))  #make new col, with counts.
      plot_tib <- plot_tib %>% mutate(thresh = case_when(no_zeros <= nz_genes ~ "TRUE", TRUE ~ "FALSE"))
      #plot scatter plot
      cols <- c("FALSE" = "blue", "TRUE" = "black")
      scatter <- ggplot(plot_tib, aes(Median, no_zeros))+
        geom_point(aes(color=thresh), alpha=0.75)+
        scale_color_manual(values = cols)+
        scale_x_log10()+
        labs(title = 'Plot of Median vs Number of Non-Zero genes', subtitle = "Filtered genes are marked with blue color")+
        theme_bw()+
        ylab('Number of samples with zero count')+
        theme(legend.position = 'bottom')
      return(scatter)}
    else{return(NULL)}
  }
  #function to produce heatmap
  plot_heatmap <- function(counts_tib, perc_var){
    if (!is.null(input$countsID)){
      counts_tib <- log10(counts_tib[-1])
      #produce plot_tib
      plot_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib, MARGIN = 1, FUN = var))
      perc_val <- quantile(plot_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
      plot_tib <- filter(plot_tib, variance >= perc_val) #filter the tibble
      hmap <- heatmap(as.matrix(plot_tib[-ncol(plot_tib)]), scale = "row", col= colorRampPalette(brewer.pal(9, "RdBu"))(25))
      legend <- legend(x = 'bottomright', legend=c('low','high'),
                       fill=colorRampPalette(brewer.pal(11, "RdBu"))(2))
      return(hmap)}
    else{return(NULL)}
  }
  #function to produce PCA plot
  plot_pca <- function(counts_tib, perc_var, comp1, comp2){
    if (!is.null(input$countsID)){
      #make plot tib-
      filt_tib <- counts_tib %>% 
        mutate(variance = apply(counts_tib[-1], MARGIN = 1, FUN = var), .after = gene)
      perc_val <- quantile(filt_tib$variance, probs = perc_var/100, na.rm = TRUE)   #calculate percentile
      filt_tib <- filter(filt_tib, variance >= perc_val) #filter the tibble
      pca_res <- prcomp(t(filt_tib[-c(1,2)]), scale = FALSE) #transpose the data and perform PCA
      #extract variance
      variance <- summary(pca_res)$importance[2,]
      x <- round(variance[comp1]*100, 2)
      y <- round(variance[comp2]*100, 2)
      #produce PCA plot
      plot_tib <- tibble(PC1 = pca_res$x[,comp1], PC2=pca_res$x[,comp2])
      pca <- ggplot(plot_tib, aes(PC1, PC2))+
        geom_point()+
        labs(title="Principle Component Analysis Plot")+
        xlab(str_c(comp1, x, "% variance", sep=" "))+
        ylab(str_c(comp2, y, "% variance", sep=" "))+
        theme_bw()
      return(pca)}
    else{return(NULL)}
  }
  #function to produce volcano plot
  volcano_plot <- function(dataf, x_name, y_name, slider, color1, color2) {
    if (!is.null(input$DE_ID)){
      vol_plot<-ggplot(data = dataf,aes(x = !!sym(x_name), y=-log10(!!sym(y_name)))) +
        geom_point(aes(color = padj< 1*10^(slider)))+
        theme(legend.position = "bottom") +
        labs( color = str_glue('{y_name} 1 x 10^ {slider}'))
      return(vol_plot) }
    else{return(NULL)}
  }
  #function to produce plot table
  draw_table <- function(dataf, slider) {
    if (!is.null(input$DE_ID)){
      df <- dataf %>% 
        dplyr::filter(padj <= 10^slider)
      df[] <- lapply(df, formatC, format="f", digits = 4)
      return(df)}
    else{return(NULL)}
  }
  #function to make distribution plots
  plot_distro <- function(counts_tib, meta_tib, meta_cat, selectgene){
    if (!is.null(input$IGE_metaID) & !is.null(input$IGE_countsID)){
      counts_tib <- column_to_rownames(counts_tib, var = "gene")
      gene_counts <- as.numeric(as.vector(counts_tib[selectgene,]))
      plot_tib <- tibble(Gene_Counts = gene_counts, meta_value = pull(meta_tib, meta_cat))
      if (meta_cat == "Diagnosis"){
        plot <- ggplot(plot_tib, aes(meta_value))+
          geom_bar()+
          theme_bw()+
          labs(title = "Plot of gene counts vs Diagnosis")
        return(plot)
      }
      else {
        plot <- ggplot(plot_tib, aes(meta_value,Gene_Counts))+
          geom_point()+
          theme_bw()+
          labs(title = str_c("Plot of gene counts vs ", meta_cat))
        return(plot)
      }}
    else{return(NULL)}
  }
  #methods to display output to front-end
  output$summ <- renderTable({
    summ_table(meta_load())
  })
  output$metadata <- DT::renderDataTable({
    meta_load()
  })
  output$aod <- renderPlot({
    plot_aod(meta_load())
  })
  output$rin <- renderPlot({
    plot_rin(meta_load())
  })
  output$Diag <- renderPlot({
    plot_Diag(meta_load())
  })
  output$normID <- renderTable({
    counts_sum(counts_load(), input$varID, input$nonzero)
  })
  output$medvar <- renderPlot({
    med_vs_var(counts_load(), input$varID)
  })
  output$medzero <- renderPlot({
    med_vs_nz(counts_load(), input$nonzero)
  })
  output$hmap <- renderPlot({
    plot_heatmap(counts_load(), input$varID)
  })
  output$PCAplot <- renderPlot({
    plot_pca(counts_load(), input$varID, input$comp1, input$comp2)
  })
  output$difftable <- DT::renderDataTable({
    de_load()
  })
  output$volcano <- renderPlot({
    volcano_plot(de_load(), input$xaxis, input$yaxis, input$padjmag, input$basecol, input$highcol)
  })
  output$IGE_plot <- renderPlot({
    plot_distro(load_IGE_counts(), IGE_meta(), input$metachoice, input$gene)
  })
}

shinyApp(ui=ui, server = server)
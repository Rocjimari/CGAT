#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.
#
# Find out more about building applications with Shiny here:
#
#    https://shiny.posit.co/
#

library(shiny)
library(shinydashboard)
library(DECIPHER)
#remotes::install_github("jonathonthill/sangerseqR", force = TRUE)
library(sangerseqR)
#remotes::install_github("mammerlin/sangeranalyseR")
library(sangeranalyseR)
library(DT)
library(tidyverse)
library(glue)

title <- tags$a(href ="https://www.google.com",
       tags$img(src="https://i.postimg.cc/qvZwS8PB/Logo-Ultimo.png", style = "width: auto; height: 100%;"))

ui <- dashboardPage(
  title = "App",
  skin = "green",
  dashboardHeader(title = tags$a(href = "https://postimg.cc/DSL1bbfR",
                                 tags$img(src="https://i.postimg.cc/qvZwS8PB/Logo-Ultimo.png",
                                          style = "width: auto; height: 100%;"))),
  dashboardSidebar(
    sidebarSearchForm(
      "searchText",
      "Introduce búsqueda",
      "Buscar",
      icon = shiny::icon("magnifying-glass-arrow-right")
    ),
    sidebarMenu(
      id = "panelID",
      menuItem("Alineamiento", 
               tabName = "plot_page",
               icon = icon("dna"),
               fileInput("files","Upload files",multiple = TRUE) )

    )
  ),
  dashboardBody(uiOutput("Alignment")
  ))
 

server <- function(input, output) {
  output$Alignment <- renderUI({
    req(input$files)
    list(
    fluidPage(
      tabItem(
        titlePanel(h3(id = "Alignment_read_1", "Alignment read 1", align = "center")),
        tags$style(
          HTML("#Alignment_read_1 { 
              color: #007f5f;
              font-family: 'Poppins';
              background: white;
              border: 1.5px solid black; 
              padding: 10px; 
            }")),
        tabName = "Alineamiento read 1",
        title = "Alineamiento", 
        fluidRow(
                 #column(width = 12,htmlOutput("alineamiento1")),
                 column(width = ifelse(!exists("error_ocurrido"), 12, 0), htmlOutput("alineamiento1")),
                 column(width = 6,fileInput("filesForward1", "Seleccione los archivos forward", multiple = TRUE)),
                 column(width = 6,fileInput("filesReverse1", "Seleccione los archivos reverse", multiple = TRUE))
                 #column(width = 12),htmlOutput("alineamiento1_manual")
                 )
      ),
      tabItem(titlePanel(
              h3(id="Alignment_read_2","Alignment read 2", align = "center")),
              tags$style(
                HTML(
                  "#Alignment_read_2 { 
                    color: #2b9348;
                    font-family: 'Poppins';
                    background: white;
                    border: 1.5px solid black; 
                    padding: 10px; 
                  }"
                )
              ),              tabName = "Alineamiento read 2",
              fluidRow(
                column(width = 12,htmlOutput("alineamiento2"))
              )),
      tabItem(titlePanel(
              h3(id="Alignment_read_3", "Alignment read 3", align = "center")),
              tags$style(
                HTML(
                  "#Alignment_read_3 { 
                    color: #d4d700;
                    font-family: 'Poppins';
                    background: white;
                    border: 1.3px solid black; 
                    padding: 10px; 
                  }"
                )
              ),
              tabName = "Alineamiento read 3",
              fluidRow(
                column(width = 12,htmlOutput("alineamiento3"))
              )),
      tabItem(titlePanel(
              h3(id="Consenso","Consenso final", align = "center")),
              tags$style(
                HTML(
                  "#Consenso { 
                    color: black;
                    font-family: 'Poppins';
                    background: #FFFFC1;
                    border: 1.5px solid black; 
                    padding: 10px; 
                  }"
                )
              ),
              tabName = "Alineamiento",
              fluidRow(
                column(width = 12, htmlOutput("alineamiento")))),
      tabItem(tabName = "numeroSeqs",
              fluidRow(
                column(width = 12,
                       tabBox(
                         width = "100%",
                         title ="Lectura archivos",
                         height = "250px",
                         id = "files_tab",
                         tabPanel(
                           "Read 1",
                           tabsetPanel(
                             tabPanel("Fowards", dataTableOutput("fw1")),
                             tabPanel("Reverse", dataTableOutput("rv1"))
                           )
                         ),
                         tabPanel(
                           "Read 2",
                           tabsetPanel(
                             tabPanel("Fowards", dataTableOutput("fw2")),
                             tabPanel("Reverse", dataTableOutput("rv2"))
                           )
                         ),
                         tabPanel(
                           "Read 3",
                           tabsetPanel(
                             tabPanel("Fowards", dataTableOutput("fw3")),
                             tabPanel("Reverse", dataTableOutput("rv3"))
                           )
                         )
                       )
              ) ))),
    fluidPage(downloadLink("downloadConsensus", "Descargar consenso")))
    })
  
  output$downloadConsensus <- downloadHandler(
    filename = function() {
      "consensus.fasta"  # Nombre del archivo que se descargará
    },
    content = function(file) {
      # Escribe el archivo consensus.fasta
      write.dna(merged.reads, file = file, format = 'fasta', nbcol = -1, colsep = "", colw = 10000000)
      
      # Lee el contenido del archivo y devuélvelo para la descarga
      content <- readLines(file)
      return(content)
    }
  )
  
  output$alineamiento1 <- renderUI({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- process_fragments1(file_info)
    
    if (is.null(result)) {
      # Si la función process_fragments1 no devuelve nada, cambia el file_info al input 
      file_info_forward <- data.frame(
        name = basename(input$filesForward1$name),
        path = input$filesForward1$datapath)
      
      file_info_reverse <- data.frame(
        name = basename(input$filesReverse1$name),
        path = input$filesReverse1$datapath)
      result <- process_fragments1_manual(file_info_forward, file_info_reverse)
      result <- result$alignment
    }
    
    htmlFile1 <- BrowseSeqs(result, openURL = FALSE)
    includeHTML(htmlFile1)
  })
  
  
  process_fragments1_manual <- function(file_info_forward, file_info_reverse) {
    if (nrow(file_info_forward) == 0 && nrow(file_info_reverse) == 0) {
      return(NULL)  # Devolver NULL si no hay datos para procesar
    }
    # Crear listas para almacenar las secuencias
    fwd_sequences <- list()
    rev_sequences <- list()
    read_names <- list()
    
    # Iterar sobre los archivos forward
    for (i in seq_len(nrow(file_info_forward))) {
      f <- file_info_forward$path[i]
      seq <- read.abif(f)
      # Quitamos bases de mala calidad de la secuencia --> recorte
      trims.seq <- trim.mott(seq, cutoff = 0.0001)
      seq.untrimmed <- seq@data$PBAS.2 
      name_1 <- file_info_forward$name[i]
      read_names[[name_1]] <- name_1
      fwd.untrimmed <- seq.untrimmed
      trims.fwd <- trims.seq
      fwd_sequences[[name_1]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
      
    }
    #return(fwd_sequences)
    # Iterar sobre los archivos reverse
    for (i in seq_len(nrow(file_info_reverse))) {
      f <- file_info_reverse$path[i]
      seq <- read.abif(f)
      # Quitamos bases de mala calidad de la secuencia --> recorte
      trims.seq <- trim.mott(seq, cutoff = 0.0001)
      seq.untrimmed <- seq@data$PBAS.2 
      name_1 <- file_info_reverse$name[i]
      read_names[[name_1]] <- name_1
      rev.untrimmed <- seq.untrimmed
      trims.rev <- trims.seq
      rev.trimmed <- substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
      complemento <- chartr("ACGT", "TGCA", rev.trimmed)
      rev_sequences[[name_1]] <- stringi::stri_reverse(complemento)
    }
    
    # Unimos secuencias fwd y rev
    fwd_reads <- DNAStringSet(unlist(fwd_sequences))
    rev_reads <- DNAStringSet(unlist(rev_sequences))
    
    reads <- DNAStringSet(c(unlist(fwd_sequences), unlist(rev_sequences)))
    
    for (f in names(read_names)) {
      if (f %in% names(reads)) {
        new_name <- read_names[[f]]
        names(reads)[names(reads) == f] <- new_name
      }
    }
    
    merged.reads <- merge.reads(reads)
    #return(merged.reads$alignment)
    print(merged.reads)
    return(merged.reads)
  }
  

  process_fragments1 <- function(file_info) {
    
    #Fragmento 1
    fragment_1_files <- file_info$path[grep(paste0("1-1F|-1F-|1f_|-1f-|_1F-|_1F-|-1-f|-1-F|1-1R|-1R-|1r_|-1r-|_1R-|_1R-|-1-R"), file_info$name)]
  
    if(length(fragment_1_files)>1){
      
      # Crear listas para almacenar las secuencias
      fwd_sequences <- list()
      rev_sequences <- list()
      read_names <- list()
            
      # Iterar sobre los archivos 1F
      for (f in fragment_1_files){
        seq <- read.abif(f)
        # Quitamos bases de mala calidad de la secuencia --> recorte
        trims.seq <- trim.mott(seq, cutoff = 0.0001)
        seq.untrimmed <- seq@data$PBAS.2 
        name_1  <- file_info$name[file_info$path == f]
        read_names[f]<- name_1
        if (str_detect(name_1,"1-1F|-1F-|1f_|-1f-|_1F|_1F-|-1-F")) {
          fwd.untrimmed <- seq.untrimmed
          trims.fwd <- trims.seq
          fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
        } else {
          rev.untrimmed <- seq.untrimmed
          trims.rev <- trims.seq
          rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
          complemento<- chartr("ACGT", "TGCA", rev.trimmed)
          rev_sequences[[f]] <- stringi::stri_reverse(complemento)
        }
      }
        
        
      # Unimos secuencias fwd y rev
      fwd_reads <- DNAStringSet(unlist(fwd_sequences))
      print(fwd_reads)
      rev_reads <- DNAStringSet(unlist(rev_sequences))
      print(rev_reads)
      
      
      reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
      print(reads)
      
      
      for (f in names(read_names)) {
        if (f %in% names(reads)) {
          new_name <- read_names[[f]]
          names(reads)[names(reads) == f] <- new_name
        }
      }
      merged.reads = merge.reads(reads)
      return(merged.reads$alignment )
      
    }
  }
  process_fragments2 <- function(file_info) {
    
    
    fragment_2_files <- file_info$path[grep(paste0("2-1F|-2F-|2f_|-2f-|_2F|_2F-|-2-F|-2F_|2-1R|-2R-|2r_|-2r-|_2R-|_2R-|-2-R|-2R_"), file_info$name)]
    if(length(fragment_2_files)>0){
      
      fwd_sequences <- list()
      rev_sequences <- list()
      read_names <- list()
      
      # Iterar sobre los archivos 1F
      for (f in fragment_2_files){
        seq <- read.abif(f)
        # Quitamos bases de mala calidad de la secuencia --> recorte
        trims.seq <- trim.mott(seq, cutoff = 0.0001)
        seq.untrimmed <- seq@data$PBAS.2 
        name_2  <- file_info$name[file_info$path == f]
        read_names[f]<- name_2
        if (str_detect(name_2,"2-1F|-2F-|2f_|-2f-|_2F|_2F-|-2-F|-2F_")) {
          fwd.untrimmed <- seq.untrimmed
          trims.fwd <- trims.seq
          fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
        } else {
          rev.untrimmed <- seq.untrimmed
          trims.rev <- trims.seq
          rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
          complemento<- chartr("ACGT", "TGCA", rev.trimmed)
          rev_sequences[[f]] <- stringi::stri_reverse(complemento)
        }
      }
      # Unimos secuencias fwd y rev
      reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
      
      for (f in names(read_names)) {
        if (f %in% names(reads)) {
          new_name <- read_names[[f]]
          names(reads)[names(reads) == f] <- new_name
        }
      }
      merged.reads = merge.reads(reads)
      return(merged.reads$alignment )
      
    }
    }
  process_fragments3 <- function(file_info) {
    
    fragment_3_files <- file_info$path[grep(paste0("3-1F|-3F-|3f_|-3f-|_3F-|_3F-|-3-F|-3F_|3-1R|-3R-|3r_|-3r-|_3R-|_3R-|-3-R|-3R_"), file_info$name)]
    print(fragment_3_files)
    if(length(fragment_3_files)>0){
      # Fragmento 3

      fwd_sequences <- list()
      rev_sequences <- list()
      read_names <- list()
      
      # Iterar sobre los archivos 1F
      for (f in fragment_3_files){
        seq <- read.abif(f)
        # Quitamos bases de mala calidad de la secuencia --> recorte
        trims.seq <- trim.mott(seq, cutoff = 0.0001)
        seq.untrimmed <- seq@data$PBAS.2 
        name_3  <- file_info$name[file_info$path == f]
        read_names[f]<- name_3
        if (str_detect(name_3,"3-1F|-3F-|3f_|-3f-|_3F|3F_|_3F-|-3F_")) {
          fwd.untrimmed <- seq.untrimmed
          trims.fwd <- trims.seq
          fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
        } else {
          rev.untrimmed <- seq.untrimmed
          trims.rev <- trims.seq
          rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
          complemento<- chartr("ACGT", "TGCA", rev.trimmed)
          rev_sequences[[f]] <- stringi::stri_reverse(complemento)
        }
      }
      # Unimos secuencias fwd y rev
      reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
      for (f in names(read_names)) {
        if (f %in% names(reads)) {
          new_name <- read_names[[f]]
          names(reads)[names(reads) == f] <- new_name
        }
      }
      
      merged.reads = merge.reads(reads)
      return(merged.reads$alignment )     
    }
  }
  

 
  
  
  output$alineamiento2 <- renderUI({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath
    )
    result <- process_fragments2(file_info)
    htmlFile2 <- BrowseSeqs(result, openURL = FALSE)
    includeHTML(htmlFile2)
  })
  output$alineamiento3 <- renderUI({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath
    )
    result <- process_fragments3(file_info)
    htmlFile3 <- BrowseSeqs(result, openURL = FALSE)
    includeHTML(htmlFile3)
  })
  
  
  
  process_fragments <- function(file_info) {
    
    #Fragmento 1
    fragment_1_files <- file_info$path[grep(paste0("1-1F|-1F-|1f_|-1f-|_1F-|1F_|-1-F|1-1R|-1R-|1r_|-1r-|_1R-|1R_|-1-R"), file_info$name)]
    
    # Crear listas para almacenar las secuencias
      fwd_sequences <- list()
      rev_sequences <- list()
      read_names <- list()
    
    # Iterar sobre los archivos 1F
    for (f in fragment_1_files){
      seq <- read.abif(f)
      # Quitamos bases de mala calidad de la secuencia --> recorte
      trims.seq <- trim.mott(seq, cutoff = 0.0001)
      seq.untrimmed <- seq@data$PBAS.2 
      name_1  <- file_info$name[file_info$path == f]
      read_names[f]<- name_1
      if (str_detect(name_1,"1-1F|-1F-|1f_|-1f-|_1F|1F_|-1-F")) {
        fwd.untrimmed <- seq.untrimmed
        trims.fwd <- trims.seq
        fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
      } else {
        rev.untrimmed <- seq.untrimmed
        trims.rev <- trims.seq
        rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
        complemento<- chartr("ACGT", "TGCA", rev.trimmed)
        rev_sequences[[f]] <- stringi::stri_reverse(complemento)
      }
    }
    # Unimos secuencias fwd y rev
    reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
    
    for (f in names(read_names)) {
      if (f %in% names(reads)) {
        new_name <- read_names[[f]]
        names(reads)[names(reads) == f] <- new_name
      }
    }
    merged.reads1 = merge.reads(reads)
    consenso1 <-merged.reads1$consensus
    if( is.null(consenso1)){
      if (length(input$filesForward1$name) > 0) {
        file_info_forward <- data.frame(
          name = basename(input$filesForward1$name),
          path = input$filesForward1$datapath
        )
      } else {
        file_info_forward <- data.frame(name = character(), path = character())
      }
      
      if (length(input$filesReverse1$name) > 0) {
        file_info_reverse <- data.frame(
          name = basename(input$filesReverse1$name),
          path = input$filesReverse1$datapath
        )
      } else {
        file_info_reverse <- data.frame(name = character(), path = character())
      }
      consenso1 = process_fragments1_manual(file_info_forward,file_info_reverse)$consensus
      print(consenso1)
    }
    #BrowseSeqs(merged.reads$alignment)
    
    
    
    # Fragmento 2
    fragment_2_files <- file_info$path[grep(paste0("2-1F|-2F-|2f_|-2f-|_2F-|2F_|-2-F|-2F_|2-1R|-2R-|2r_|-2r-|_2R-|2R_|-2-R|-2R_"), file_info$name)]
    
    read_names <- list()
    fwd_sequences <- list()
    rev_sequences <- list()
    
    # Iterar sobre los archivos 1F
    for (f in fragment_2_files){
      seq <- read.abif(f)
      # Quitamos bases de mala calidad de la secuencia --> recorte
      trims.seq <- trim.mott(seq, cutoff = 0.0001)
      seq.untrimmed <- seq@data$PBAS.2 
      name_2  <- file_info$name[file_info$path == f]
      read_names[f]<- name_2
      
      if (str_detect(name_2,"2-1F|-2F-|2f_|-2f-|_2F|2F_|-2-F|-2F_|-2F_")) {
        fwd.untrimmed <- seq.untrimmed
        trims.fwd <- trims.seq
        fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
      } else {
        rev.untrimmed <- seq.untrimmed
        trims.rev <- trims.seq
        rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
        complemento<- chartr("ACGT", "TGCA", rev.trimmed)
        rev_sequences[[f]] <- stringi::stri_reverse(complemento)
      }
    }
    # Unimos secuencias fwd y rev
    reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
    
    for (f in names(read_names)) {
      if (f %in% names(reads)) {
        new_name <- read_names[[f]]
        names(reads)[names(reads) == f] <- new_name
      }
    }
    merged.reads2 = merge.reads(reads)
    consenso2 <-merged.reads2$consensus
    #BrowseSeqs(merged.reads$alignment)
    print(consenso2)
    
    
    # Fragmento 3
    fragment_3_files <- file_info$path[grep(paste0("3-1F|-3F-|3f_|-3f-|_3F-|3F_|-3-F|-3F_|3-1R|-3R-|3r_|-3r-|_3R-|3R_|-3-R|-3R_"), file_info$name)]
    
    fwd_sequences <- list()
    rev_sequences <- list()
    read_names <- list()
    
    
    # Iterar sobre los archivos 1F
    for (f in fragment_3_files){
      seq <- read.abif(f)
      # Quitamos bases de mala calidad de la secuencia --> recorte
      trims.seq <- trim.mott(seq, cutoff = 0.0001)
      seq.untrimmed <- seq@data$PBAS.2 
      name_3  <- file_info$name[file_info$path == f]
      read_names[f]<- name_3
      
      if (str_detect(name_3,"3-1F|-3F-|3f_|-3f-|_3F|3F_|-3-F|-3F_")) {
        fwd.untrimmed <- seq.untrimmed
        trims.fwd <- trims.seq
        fwd_sequences[[f]] <- substring(fwd.untrimmed, trims.fwd$start, trims.fwd$finish)
      } else {
        rev.untrimmed <- seq.untrimmed
        trims.rev <- trims.seq
        rev.trimmed = substring(rev.untrimmed, trims.rev$start, trims.rev$finish)
        complemento<- chartr("ACGT", "TGCA", rev.trimmed)
        rev_sequences[[f]] <- stringi::stri_reverse(complemento)
      }
    }
    # Unimos secuencias fwd y rev
    reads <- DNAStringSet(c(unlist(fwd_sequences),unlist(rev_sequences)))
    
    for (f in names(read_names)) {
      if (f %in% names(reads)) {
        new_name <- read_names[[f]]
        names(reads)[names(reads) == f] <- new_name
      }
    }
    
    merged.reads3 = merge.reads(reads)
    consenso3 <-merged.reads3$consensus
    #BrowseSeqs(merged.reads$alignment)  
    print(consenso3)
    
    if (is.null(consenso1) && is.null(consenso2)) {
      return(merged.reads3$alignment)
      
    } else if (is.null(consenso1) && is.null(consenso3)) {
      return(merged.reads2$alignment)
      
    } else if (is.null(consenso2) && is.null(consenso3)) {
      return(merged.reads1$alignment)
      
    } else{
      # Unión consenso  
      reads = DNAStringSet(c(as.character(consenso1),as.character(consenso2),as.character(consenso3)))
      merged.reads = merge.reads(reads)
      #BrowseSeqs(merged.reads$alignment)
      print(merged.reads)
      return(merged.reads$alignment)
      
    }
    
  }
  
  
  count_files1 <- function(file_info) {
    fragment_1_files <- file_info$path[grep(paste0("1-1F|-1F-|1f_|-1f-|_1F-|-1-F|1-1R|-1R-|1r_|-1r-|_1R-|-1-R"), file_info$name)]
    
    fwd1 <- list()
    rev1 <- list()
    
    for (f in fragment_1_files) {
      name_1 <- file_info$name[file_info$path == f]
      
      if (str_detect(name_1, "1-1F|-1F-|1f_|-1f-|_1F|-1-F")) {
        fwd1[[f]] <- name_1
      } else {
        # Store in the second list if condition is not met
        rev1[[f]] <- name_1
      }
    }
    
    return(list(fwd1 = fwd1, rev1 = rev1))
  }  
  count_files2 <- function(file_info) {
    fragment_2_files <- file_info$path[grep(paste0("2-1F|-2F-|2f_|-2f-|_2F-|-2-F|-2F_|2-1R|-2R-|2r_|-2r-|_2R-|-2-R|-2F_"), file_info$name)]
    
    fwd2 <- list()
    rev2 <- list()
    
    for (f in fragment_2_files) {
      name_2 <- file_info$name[file_info$path == f]
      
      if (str_detect(name_2, "2-1F|-2F-|2f_|-2f-|_2F|-2-F|-2F_")) {
        fwd2[[f]] <- name_2
      } else {
        # Store in the second list if condition is not met
        rev2[[f]] <- name_2
      }
    }
    
    return(list(fwd2 = fwd2, rev2 = rev2))
  }  
  count_files3 <- function(file_info) {
    fragment_3_files <- file_info$path[grep(paste0("3-1F|-3F-|3f_|-3f-|_3F-|-3-F|-3F_|3-1R|-3R-|3r_|-3r-|_3R-|-3-R|-3R_"), file_info$name)]
    
    fwd3 <- list()
    rev3 <- list()
    
    for (f in fragment_3_files) {
      name_3 <- file_info$name[file_info$path == f]
      print("name_3")
      print(name_3)
      
      if (str_detect(name_3, "3-1F|-3F-|3f_|-3f-|_3F|-3-F|-3F_")) {
        fwd3[[f]] <- name_3
        print("fwd3")
        print(fwd3)
      } else {
        # Store in the second list if condition is not met
        rev3[[f]] <- name_3
        print("rev3")
        print(rev3)
      }
    }
    
    return(list(fwd3 = fwd3, rev3 = rev3))
  }
  
  output$alineamiento <- renderUI({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath
    )
    result <- process_fragments(file_info)
    htmlFile <- BrowseSeqs(result, openURL = FALSE)
    includeHTML(htmlFile)
  })
  
  
  output$fw1 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files1(file_info)
    #if( length(result) == 0){
     # result <- data.frame(
      #  name = basename(input$filesForward1$name),
       # path = input$filesForward1$datapath)
      #fwd1 <- unlist(fwd1_list)
    #  return(data.frame(FileName = fwd1))  
      
      
    #}else{
      fwd1_list <- result$fwd1
      fwd1 <- unlist(fwd1_list)
      return(data.frame(FileName = fwd1)) 
    #}
    
      
    },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
    rownames= FALSE)
  output$rv1 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files1(file_info)
    
    rev1_list <- result$rev1
    rev1 <- unlist(rev1_list)
    return(data.frame(FileName = rev1))   
    
  },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
  rownames= FALSE)
  
  
  
  output$fw2 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files2(file_info)
    
    fwd2_list <- result$fwd2
    fwd2 <- unlist(fwd2_list)
    return(data.frame(FileName = fwd2))   
  },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
  rownames= FALSE)
  output$rv2 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files2(file_info)
    
    rev2_list <- result$rev2
    rev2 <- unlist(rev2_list)
    return(data.frame(FileName = rev2))   
    
  },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
  rownames= FALSE)
  

  
  output$fw3 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files3(file_info)
    
    fwd3_list <- result$fwd3
    fwd3 <- unlist(fwd3_list)
    return(data.frame(FileName = fwd3))   
  },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
  rownames= FALSE)
  output$rv3 <- renderDT({
    file_info <- data.frame(
      name = basename(input$files$name),
      path = input$files$datapath)
    result <- count_files3(file_info)
    
    rev3_list <- result$rev3
    rev3 <- unlist(rev3_list)
    return(data.frame(FileName = rev3))   
    
  },options = list(pageLength = 5, autoWidth = TRUE,info = FALSE, paging = FALSE),
  rownames= FALSE)
  
  
  
}
    

shinyApp(ui = ui, server = server)

library(shiny)
library(shinydashboard)
library(DT)
library(shinyjs)
library(sodium)
library(data.table)
library(openxlsx)
library(shinyWidgets)

germline_df = as.data.table(read.xlsx('germlineNUMTs.xlsx'))
cancer_df = as.data.table(read.xlsx('cancerNUMTs.xlsx'))

function(input, output, session) {
  
  login = FALSE
  USER <- reactiveValues(login = login)
  
  observe({ 
    if (USER$login == FALSE) {
      if (!is.null(input$login)) {
        if (input$login > 0) {
          Username <- isolate(input$userName)
          Password <- isolate(input$passwd)
          if(length(which(credentials$username_id==Username))==1) { 
            pasmatch  <- credentials["passod"][which(credentials$username_id==Username),]
            pasverify <- password_verify(pasmatch, Password)
            if(pasverify) {
              USER$login <- TRUE
            } else {
              shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
              shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
            }
          } else {
            shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade")
            shinyjs::delay(3000, shinyjs::toggle(id = "nomatch", anim = TRUE, time = 1, animType = "fade"))
          }
        } 
      }
    }    
  })
  
  output$logoutbtn <- renderUI({
    req(USER$login)
    tags$li(a(icon("fa fa-sign-out"), "Logout", 
              href="javascript:window.location.reload(true)"),
            class = "dropdown", 
            style = "background-color: #eee !important; border: 0;
            font-weight: bold; margin:5px; padding: 10px;")
  })
  
  output$sidebarpanel <- renderUI({
    if (USER$login == TRUE ){ 
      sidebarMenu(
        menuItem("Sample Information", tabName = "dashboard", icon = icon("info")),
        menuItem("Germline NUMTs", tabName = "GermlineNUMTs", icon = icon("table")),
        menuItem("Cancer NUMTs", tabName = "CancerNUMTs", icon = icon("table"))
      )
    }
  })
  
  output$body <- renderUI({
    if (USER$login == TRUE ) {
      tabItems(
        # First tab
        tabItem(tabName ="dashboard", class = "active",
                fluidRow(
                  tabBox(
                    title = 'Sample information',br(),br(),
                    id = 'sampleInfor', height = '250px',
                    tabPanel(title = "germline", width = 4, status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,icon = icon("dna"),
                             img(src=b64_germline, height = 250)
                    ),
                    tabPanel(title = "cancer", width = 4, status = "primary", solidHeader = TRUE,
                             collapsible = TRUE,icon = icon("dna"),
                             img(src=b64_cancer, height = 350)
                    )
                  )
                )
        ),
        
        # Second tab
        tabItem(tabName = "GermlineNUMTs",
                fluidRow(
                  #box(width = 12, dataTableOutput('GermlineTable'))
                  sidebarLayout(
                    sidebarPanel(width = 3,
                                 pickerInput("populations", "Select population(s):", population_list,
                                             multiple = TRUE, options = list(`actions-box` = TRUE), selected = population_list
                                 ),
                                 
                                 checkboxGroupInput("frequency", "Select frequency (full dataset):", frequency_list,
                                                    selected = frequency_list
                                 ),br(),
                                 
                                 selectizeInput("chromosome", "Select chromosome(s):", chromosome_list, multiple = TRUE,
                                                selected = chromosome_list,
                                                options = list(
                                                  plugins = list("remove_button")
                                                )
                                 ),
                                 actionButton("Select_All", label = "Select_All"),
                                 actionButton("Deselect_All", label = "Deselect_All"), br(),br(),
                                 
                                 #selectInput("mtGene", "Select mito region(s):", mtGene_list, multiple = TRUE),
                                 
                                 pickerInput("putative_concatenated_NUMTs", "concatenatedNUMTs:", NUMTtype_list,
                                             multiple = TRUE, options = list(`actions-box` = TRUE), selected = NUMTtype_list
                                 ),br(),
                                 
                                 pickerInput("Longread_NUMTs", "Long-read validated NUMTs:", longRead_list,
                                             multiple = TRUE, options = list(`actions-box` = TRUE), selected = longRead_list
                                 ),br()
                                 
                                 #textInput("downloadGermline", "Enter download file name:"),
                                 #downloadButton("downloadGermline", "Download germline NUMTs")
                    ),
                    mainPanel(
                      fluidRow(
                        column(width=12,br(),
                               dataTableOutput('germlineOut'),br(),br(),br()
                               #h3("Overlap view plots")
                        )
                      )
                    )
                  )
                )
        ),
        
        # Third tab
        tabItem(tabName = "CancerNUMTs",
                fluidRow(
                  sidebarLayout(
                    sidebarPanel(width = 3,
                                 selectizeInput("cancerType", "Select cancer type(s):", cancerType_list, multiple = TRUE,
                                                selected = cancerType_list,
                                                options = list(
                                                  plugins = list("remove_button")
                                                )),
                                 actionButton("Select_Allcancers", label = "Select_All"),
                                 actionButton("Deselect_Allcancers", label = "Deselect_All"), br(),br(),
                                 
                                 #selectizeInput("chromosome", "Select chromosome(s):", chromosome_list, multiple = TRUE), 
                                 
                                 #selectizeInput("mtGene", "Select mito region(s):", mtGene_list, multiple = TRUE),
                                 #pickerInput("cancerGenes", "Select cancer gene(s):", cancerGene_list,
                                 #             multiple = TRUE, options = list(`actions-box` = TRUE), selected = cancerGene_list
                                 #),br(),
                                 
                                 pickerInput("complexNUMTs", "Select NUMTs complexity:", complexNUMTs_list,
                                             multiple = TRUE, options = list(`actions-box` = TRUE), selected = complexNUMTs_list
                                 ),br()
                                 
                                 #textInput("downloadCancer", "Enter download file name:"),
                                 #downloadButton("downloadCancer", "Download cancer-specific NUMTs")
                                 
                    ),             
                    mainPanel(
                      fluidRow(
                        column(width=12,br(),
                               dataTableOutput('cancerOut'),br(),br(),br()
                               #h3("Overlap view plots")
                        )
                      )
                    )
                  )
                )
        )
      )
      
    }
    else {
      loginpage
    }
  })
  
  observeEvent(input$Select_All,{
    updateSelectInput(session, inputId = "chromosome", choices = chromosome_list, selected = chromosome_list)
  })
  observeEvent(input$Deselect_All,{
    updateSelectInput(session, inputId = "chromosome", choices = chromosome_list, selected = "")
  })
  
  observeEvent(input$Select_Allcancers,{
    updateSelectInput(session, inputId = "cancerType", choices = cancerType_list, selected = cancerType_list)
  })
  observeEvent(input$Deselect_Allcancers,{
    updateSelectInput(session, inputId = "cancerType", choices = cancerType_list, selected = "")
  })
  
  
  output$germlineOut = renderDataTable({
    
    germlineTable <- germline_df
    populationPattern <- paste(input$populations, collapse="|")
    germlineTable_out <- germlineTable[Frequency_group %in% input$frequency]
    germlineTable_out <- subset(germlineTable_out, grepl(populationPattern, Ethnicity))
    germlineTable_out <- germlineTable_out[chromosome %in% input$chromosome]
    germlineTable_out <- germlineTable_out[concatenatedNUMTs %in% input$putative_concatenated_NUMTs]
    germlineTable_out <- germlineTable_out[longRead_validation %in% input$Longread_NUMTs]
    
    datatable(
      germlineTable_out,
      filter = 'top', 
      escape = FALSE,
      options = list(
        options = list(searchHighlight = F), 
        scrollX=TRUE,
        scrollY=600,
        pageLength = 25)
    )
  })
  
  output$cancerOut = renderDataTable({
    
    cancerTable <- cancer_df
    cancerTable_out <- cancerTable[tumourType %in% input$cancerType]
    #cancerTable_out <- cancerTable_out[complexNUMTs %in% c('UN','complexNUMTs')]
    cancerTable_out <- cancerTable_out[complexNUMTs %in% input$complexNUMTs]
    #germlineTable_out <- germlineTable_out[longRead_validation %in% input$Longread_NUMTs]
    
    datatable(
      cancerTable_out,
      filter = 'top', 
      escape = FALSE,
      options = list(
        options = list(searchHighlight = F), 
        scrollX=TRUE,
        scrollY=600,
        pageLength = 25)
    )
  })
}  
  #output$CancerTable <-  DT::renderDataTable({
  #  datatable(cancer_df, options = list(autoWidth = TRUE,
  #                                   searching = FALSE))
  #})
  
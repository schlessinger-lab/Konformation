## app.R ##
library(shiny)
library(bio3d)
library(shinyWidgets)
library(DT)
library(shinycssloaders)
source("predict.conformation.R")
source("dfgmod.R")
library(plotly)
library(randomForest)
library(clusterSim)
library(tidyverse)
library(dplyr)
library(data.table)
library(shinythemes)
library(bsplus)
library(shinyBS)
#library(rcdk)
addline_format <- function(x, ...) {
  gsub('\\s', '\n', x)    
}  
frag.annot.merged = fread("./fragments.merged.csv") 
fragments = as.character(unique(as.character(frag.annot.merged$Fragment)))
frag.annot.merged$`Odds Ratio` = as.character(frag.annot.merged$`Odds Ratio`)
frag.annot.merged$`Log Fold Change` = as.character(frag.annot.merged$`Log Fold Change`)
frag.annot.merged = frag.annot.merged[,c(1:3,8,6, 4:5, 10,11)]
row.names(frag.annot.merged) =  c(1:nrow(frag.annot.merged))
frag.annot = read.table("./fragments.annotated.csv.2",
                        header = T,
                        sep = ",")
ligand.pdb = read.table("./pdb.ligand.smi.csv.complete.cases",
                        header = F ,
                        sep = ",")
ligand.lib = read.table("./ligand.library.csv", sep = "," , header = T, stringsAsFactors = F)
colorst =  c("#7C2600", "#CE5A28" , "#00D2F1" , "#AA96DA" , "#7F7F7F")
data = read.table("./full.data.px.vec.merged.csv", sep = "," , header = T, stringsAsFactors = F)
data$key = data$PDBID
proteinsel = as.character(data[,39]) 
coln = gsub(".y$" , "", gsub(".x$", ".normalized", names(data)))
coln[1] = "PDBid"
colnames(data) = coln
twitter <- 
  "https://twitter.com/intent/tweet?text=Check%20out%20KinaMetrix%20by%20the%20@schlessingerlab%20to%20classify%20kinases!&url=http://kinametrix.com"
featlist  = coln[c(2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30)]

ui <-    
  navbarPage(  
    theme = shinytheme("flatly"), 
    title = HTML("<a href=\"http://kinametrix.com\">Kina<b>Metrix</b></a>"),
    windowTitle = "KinaMetrix.com",  
     id = "nav",
    position = "fixed-top",  
    fluid = T ,
    collapsible = T, 
    #########################################################
    ##########################################################
    tabPanel(
      "Home",
      
      h2(
        "   Kinases are dynamic proteins that can adopt several distinct conformational states.",
        align = "center"
      ),

      includeHTML("./home.html"),  
      h3(HTML("Thats why we built Kina<b>Metrix</b>."), align = "center"),
      h4( 
        "A webserver to investigate kinase conformations and inhibitor space.",
        align = "center"
      ),
      tags$hr(),
      h4(HTML("On Kina<b>Metrix</b>.com you can:"), align = "center"),
      tags$br(),
      fluidRow(
        column(2),
        column(
          4,
          actionBttn(
            inputId = "kon",
            label =
              h4(
                "Investigate our database of kinase structures annotated with their conformation.",
                align = "center",
                style = "cursor:pointer;"
              ),
            size = "lg",
            icon = icon("search", lib = "font-awesome") ,
            color = "primary"
          )
        ),
        column(
          4,
          actionBttn(
            inputId = "kin",
            label =
              h4(
                HTML(
                  "Upload a custom kinase structure and have it's conformation classified by <i>Kinformation</i>."
                ),
                align = "center"
              ),
            size = "lg",
            icon = icon("check-circle", lib = "font-awesome"),
            color = "warning"
          )
          
        )
      ),
      tags$br(),
      fluidRow(
        column(2),
        column(
          4,
          actionBttn(
            inputId = "frag",
            label =
              h4(
                HTML(
                  "Explore our library of over 10,000 chemical fragments associated with conformation."
                ),
                align = "center"
              ),
            size = "lg",
            icon = icon("database", lib = "font-awesome"),
            color = "royal"
          )
        ),
        column(
          4,
          actionBttn(
            inputId = "models",
            label =
              h4(
                HTML(
                  "\nObtain homology models of kinases in inactive conformation states using DFGModel."
                ),
                align = "center"
              ),
            size = "lg",
            icon = icon("compass", lib = "font-awesome"),
            color = "success"
          )
        ),
        column(2)
      ),
      tags$hr(),
      includeHTML("./home2.html") ,
      tags$hr(),
      h4("Other sites you should visit to study kinases!", align = "center"),
      h4(
        HTML("<a href=\"http://klifs.vu-compmedchem.nl\">KLIFS</a>"),
        HTML(
          "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
        ),
        HTML(
          "<a href=\"https://www.ebi.ac.uk/chembl/sarfari/kinasesarfari\">KinaseSARfari</a>"
        ),
        HTML(
          "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
        ),
        HTML("<a href=\"http://kinhub.org/kinmap/\">KinMap</a>"),
        align = "center"
      ),
      tags$br(),
      h4(
        HTML("For any questions please head over to our <b>FAQ</b>"),
        align = "center"
      ),
      h4(
        HTML(
          "<a href=\"http://schlessingerlab.org\">Schlessinger Lab</b></a>"
        ),
        HTML(
          "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
        ),
        tags$a(href = twitter, "Tweet us out!", class =
                 "twitter-share-button"),
        #includeScript("http://platform.twitter.com/widgets.js"),
        HTML(
          "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
        ),
        HTML(
          "Made with ❤ using <a href=\"https://shiny.rstudio.com\">Shiny</a>"
        ),
        HTML(
          "&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;"
        ),
        tags$a(
          img(src = "http://icahn.mssm.edu/sites/ismms/assets/images/IcahnLogo.png"),
          href = "https://icahn.mssm.edu"
        ),
        align = "center"
      )
    ),


    ##############################
    tabPanel(
      "Data explorer",
      id = "datapanel",
      tabsetPanel(type = "tabs",
      id = "datapanel",
        tabPanel(
          "Search Classified Kinase Structures", 
          id = "kinclass",
          tags$style(type = "text/css", "body {padding-top: 70px;}"),
          h3(
            HTML(
              "Database of Kinase Conformation Predictions From <i>Kinformation</i>"
            ),
            align = "center"
          ),
          withSpinner(dataTableOutput("datat"), type = 8),
          bsModal("nglview", "", "go", htmlOutput("frame"), size = "medium")
        ),
        tabPanel(
          "Search Geometric Descriptors", 
          tags$style(type = "text/css", "body {padding-top: 70px;}"),
          h3(
            HTML(
              "Database of Geometric Descriptors Describing Kinase Conformations </i>"
            ),
            align = "center"
          ),
          dataTableOutput("geom")
        ), 
        tabPanel(
            "Search Kinase Bound Ligands", 
            h3("Database of Kinase Bound Ligands", align = "center"),
            fluidRow(
              tags$br(),
              column( 8,
                dataTableOutput("liglib")
              ),
              column( 4,
                plotOutput( outputId = "ligandsmirender", 
                  height = "400px",
                  width = "375px"
                )
              )
            )
        ),
        tabPanel(
          "Visualize the Kinase Conformational Space",
          h3(
            "3D scatterplot of kinase structure predictions using Kinformation",
            align = "center"
          ),
          sidebarLayout(
                  sidebarPanel(
                    selectInput(
                      "featureInputgenename", 
                      label = "Select a Protein:" , 
                      choices = c("all", proteinsel ), 
                      selected = "all"
                    ), 
                    selectInput(
                      inputId = "featureInput4",
                      label = "Please Choose data type:",
                      choices = c("all", "training", "test"),
                      selected = "all"
                    ),
                    selectInput(
                      inputId = "featureInput1",
                      label = "Select first feature",
                      choices = featlist,
                      selected = featlist[11]
                    ),
                    selectInput(
                      inputId = "featureInput2",
                      label = "Select second feature",
                      choices = featlist ,
                      selected = featlist[12]
                    ),
                    
                    selectInput(
                      inputId = "featureInput3",
                      label = "Select third feature",
                      choices = featlist ,
                      selected = featlist[13]
                    ),
                    
                    uiOutput("click")
                  ),
                  mainPanel(withSpinner(
                    plotlyOutput(
                      "scaplott",
                      height = "600px",
                      width = "750px",
                      inline = T
                    ),
                    color = "#0dc5c1",
                    type = 8
                 ))
          )
        )
      )
    ),
    ###############################
    tabPanel(
      HTML("<i>Kinformation</i>"),
      tags$head(tags$style(
        ".progress-bar{background-color:green;}"
      )),
      h3(
        HTML("Annotate Kinase Conformation using <i>Kinformation</i>"),
        align = "center"
      ),
      sidebarLayout( 
        sidebarPanel(
          fileInput(
            inputId = "file1",
            label = h4("Upload Kinase Structure",
                       tags$style(type = "text/css", "#q1 {vertical-align: top;}"),
                       bsButton("q1", label = "", icon = icon("question-circle", lib = "font-awesome"), style = "info", size = "extra-small")
            ),
            width = "400px",
            buttonLabel = "Browse:",
            multiple = FALSE,
            accept = c(".pdb")
          ),
          checkboxInput("pairwise",
                        p("Compute pairwise geometric similarity?", 
                          bsButton("q2", label = "", icon = icon("question-circle", lib = "font-awesome"), style = "info", size = "extra-small"))
                        , FALSE),
          bsTooltip(id = "q2", title = HTML("This will compute the euclidean distance between the uploaded structure and all other kinase structures available. <br><i> <b>This will increase the computational time. </b></i></br> "),
                    placement = "right" , 
                    trigger = "hover", 
                    options = list(container = "body")
          ),
          bsTooltip(id = "q1", title = HTML("Please upload a <b>typical</b> kinase structure in PDB format"),
                    placement = "right" , 
                    trigger = "hover", 
                    options = list(container = "body")
          ),
          useSweetAlert(),
          tags$hr(),
          span(textOutput("jobid"), style = "color:black"),
          dataTableOutput("contents"), 
          bsTooltip(id = "contents", title = HTML("Each row displays the probablity of each conformation for the uploaded kinase structure"),
                    placement = "right" , 
                    trigger = "hover", 
                    options = list(container = "body")
          ), 
          uiOutput('ui.action'),
          tags$br(), 
          uiOutput('ui.actiontwo'),
          uiOutput('ui.actionthree'),
          uiOutput('ui.actionfive')
        ),
        mainPanel(
          imageOutput("kclass"),
          plotlyOutput(
            "scaplot" ,
            height = "500px",
            width = "750px",
            inline = T
          ), 
          dataTableOutput("similar"),
          bsTooltip(id = "similar", 
                    title = HTML(" This table compares all 
                                 of the available kinase structures to the uploaded kinase structure based on the geometric descriptors
                                 in terms of their euclidean distance. The smaller the distance, the more similar the Kinase structure
                                 is to the uploaded structure"),
                    placement =  "bottom" , 
                    trigger = "hover", 
                    options = list(container = "body")
                    )
          )
        )
    ),
    ##################
    tabPanel( 
      "Fragment Libraries",
      tabsetPanel(type = "tabs",
        tabPanel(
          "Interactive figure",
          h3("Enrichment of Substructures per Kinase Conformation" , align = "center"),
          tags$blockquote(
            "The red line indicates enriched substructures below a p-value threshold of .05 "
          ),
          tags$blockquote(
            "Click on a point to visualize a substructure of interest as well as get links to the parent ligand and interaction maps seen in the structure (scroll down)"
          ), 
          sidebarLayout(
            sidebarPanel(withSpinner(
              plotlyOutput(
                outputId = "timeseries",
                height = "600px",
                width = "425px"
              ),
              type = 8
            )
            , width = 5),
            mainPanel(
              withSpinner(
                plotOutput(
                  outputId = "correlation",
                  height = "200px",
                  width = "200px"
                ),
                type = 8
              ),
              withSpinner(dataTableOutput("cor"), type = 8),
              width = 5
            )
          )
        ), 
        ###############
        tabPanel("Search by SMILES", 
          sidebarLayout(
            sidebarPanel(
              h5("Search for fragments by SMILES string"),
              textInput("insmi", "Paste a valid SMILES string to search", ""),
              uiOutput('ui.smi.action.one')   , 
              bsTooltip(id = "searchtype", title = HTML("Please use <b> Maximum Common Substructure</b> when searching with a drug-like molecule <br>Use <b> Tanimoto Distance </b> when searching with a fragment-like molecule</br>"),
                  placement = "right" , 
                  trigger = "hover", 
                  options = list(container = "body")
              ),
              actionButton("smisubmit", label = "Submit"),
              tags$hr(),
              uiOutput('ui.smi.action.two')   , 
              uiOutput('ui.smi.action.draw.in.smi') ,
              uiOutput('ui.smi.action.three')   , 
              uiOutput('showselsmi')
            ),
          mainPanel(
            dataTableOutput("smirend")
          )
        )
      )
    )  
    ),
    #############
    tabPanel("Models using DFGModel",
             dfgui() 
             ) ,
    tabPanel(
      "FAQ",
      h3("How to cite us:"),
      h5(
        "There are currently three articles that KinaMetrix is a part of, please cite all three if possible. " ),
      h5("If you utilized this website by accessing the database,
        using Kinformation, or viewing our fragment database please refer to our following papers:"
      ),
      h5(
        HTML(
          "<ul> <li> <a href=\"https://www.cell.com/cell-chemical-biology/fulltext/S2451-9456(18)30149-1\">Ung, Rahman and Schlessinger; Redefining the Protein Kinase Conformational Space with Machine Learning</a> (2018) </li> 
          <li><a href=\"https://academic.oup.com/nar/advance-article/doi/10.1093/nar/gky916/5128920\">Rahman, Ung and Schlessinger; KinaMetrix: A Webserver to Investigate Kinase Conformations and Inhibitor Space</a> (2018) </li></ul>"
        )
        ),
      h5(
        "If you used any of the homology models of kinases in missing inactive conformations please cite our KinaMetrix article as well as the following: "
      ),
      h5(
        HTML(
          "<ul> <li> <a href=\"https://pubs.acs.org/doi/abs/10.1021/cb500696t\">Ung and Schlessinger; DFGmodel: Predicting Protein Kinase Structures in Inactive States for Structure-Based Discovery of Type-II Inhibitors</a> (2015) </li>"
        )
      ),
      h3("How to contact us:"),
      h5(
        "If you have any issues with the webserver in general, running Kinformation, accesing the fragment library, find any bugs, or have an idea for a helpful feature, please contact:"
      ),
      h5(
        HTML(
          "<ul><li><b>rayees(dot)rahman@icahn(dot)mssm(dot)edu</b>
          <br>Rayees is the principal server developer and maintainer. He also trained Kinformation and built the fragment libaries and network.</br></li><ul>"
        )
        ),
      h5(
        "For any questions about the kinase homology models, the geometric descriptors used to train Kinformation or about DFGModel, please contact:"
      ),
      h5(
        HTML(
          "<ul><li><b>peter(dot)ung@mssm(dot)edu</b>
          <br>Dr. Ung develped the algorithm to derive the geometric descriptors describing the DFG and αC-Helix motifs in kinases. <br> He also created DFGModel and generated all the homology models available on this webserver.</br></br></li><ul>"
        )
        ),
      h5(
        "For general inquires please contact:"
      ),
      h5(
        HTML(
          "<ul><li><b>avner(dot)schlessinger@mssm(dot)edu</b>
          <br>Dr. Schlessinger is the principal investigator of this project.</br></li></ul>"
        )
        ),
      h3("Our funding sources:"),
      h5("This server is funded by:"),
      h5(
        HTML(
          "<ul><li>National Institute of General Medical Sciences Integrated Pharmacological Sciences Training Program (T32GM062754) [Rayees Rahman]</li></ul>"
        )
        ),
      h3("Software, databases and packages we use:"),
      HTML(
        "
        <a href=\"http://www.rcsb.org/\">
        <img src=\"https://cdn.rcsb.org/rcsb-pdb/v2/common/images/rcsb_logo.png\" height=100 width=250> </a>
        <a href=\"https://blast.ncbi.nlm.nih.gov/Blast.cgi\">
        <img src=\"https://blast.ncbi.nlm.nih.gov/images/protein-blast-cover.png\" height=100 width=250> </a>
        <a href=\"https://www.r-project.org/\">
        <img src=\"https://www.r-project.org/Rlogo.png\"> </a>
        <a href=\"https://shiny.rstudio.com/\">
        <img src=\"https://www.rstudio.com/wp-content/uploads/2014/04/shiny.png\"height=200 width=200> </a>
        <a href=\"https://pymol.org/2/\">
        <img src=\"https://upload.wikimedia.org/wikipedia/commons/thumb/8/87/PyMOL_logo.svg/2000px-PyMOL_logo.svg.png\" height=200 width=200> </a>
        <br></br>
        <a href=\"https://github.com/arose/nglview\">nglview</a>
        <a href=\"https://www.python.org/\">
        <img src=\"https://www.python.org/static/img/python-logo@2x.png\" height=100 width=300>> </a>
        <a href=\"http://www.rdkit.org\">
        <img src=\"http://www.rdkit.org/Images/logo.png\"> </a>
        
        "
      )
      
      
      
      )
    
      )
drawsmi = function(smiles)
{
    #smi = as.character(smiles)
    #smi = gsub("\\[n", "\\[N", smi)
    #mol = parse.smiles(smi)
    mol = smiles
    par(mar = c(0, 0, 0, 0))
    d = get.depictor(
      height = 1000,
      width = 1000,
      zoom = 6,
      style = "cow"
    )
    print(mol[[1]])
    i = view.image.2d(mol[[1]], depictor = d)
    plot(
      NA,
      NA,
      xlim = c(0, 1000),
      ylim = c(0, 1000),
      xaxt = 'n',
      yaxt = 'n',
      xlab = '',
      ylab = ''
    )
    rasterImage(i, 1, 1, 1000, 1000)
}
da = data[, c(35,36, 37,38, 39, 41, 18, 20,  40, 42, 43 , 44 ,45, 19, 2:14)]
da[,7] = as.numeric(as.character(da[,7]))
da[,8] =  as.factor(da[,8])
frag.merge.sub = frag.annot.merged

##@@@@@@@@@@@@@@@@@@@@@@@#########################################################
##@@@@@@@@@@@@@@@@@@@@@@@###########################################################
##@@@@@@@@@@@@@@@@@@@@@@@#############################################################
server <- function(input, output, session) {
  output$downloadmodels <- downloadHandler(
        filename = function() {      
          name = models[input$models_row_last_clicked,2]
          fn = paste0(name, ".models.tar.bz2")
        } ,
        content = function(file) {
          name = models[input$models_row_last_clicked,2]
          upid = models[input$models_row_last_clicked,1]
          path = paste0("./1_models/cido/", name, "/cido.", upid, ".", name,".top_pdb.tar.bz2" )
          file.copy(path, file)
         } , 
         contentType = "application/tar.bz2"
      )
  output$models = renderDT({
    names(models) = c("Uniprot ID", "Gene Name", "Modeled Conformation", "Download")
    models$Download = NULL
    datatable(models, plugins = 'searchHighlight', 
              rownames= FALSE, escape = FALSE,  class = 'cell-border stripe',  
              style = 'bootstrap', 
              selection = 'single' ,
              options = list(pageLength = 15, 
                #columnDefs = list(list(targets = c(9:10,13), searchable = FALSE)),
                autoWidth = T,      
                #scrollX =T,  
                lengthMenu = c(5, 10, 15 , 30, 50 , nrow(models)),
                searchHighlight = TRUE
                )
      )
    })

  
  
 #       filename = function() {"descriptors.csv"} ,
  #      content = function(file) { write.csv(desc, file, row.names = TRUE, sep = ",") } 
  #    )
  output$smirend = renderDataTable({
       datatable(frag.annot.merged ,  caption = HTML("For more information about the table columns see the <b>FAQ</b>" ),
              rownames= FALSE, escape = FALSE, selection = 'single' , extensions = c( 'Buttons','Responsive'), 
              options = list(pageLength = 10, autoWidth = F,      
              dom = 'Bfrtip',
              buttons = c('csv',  I('colvis')))) })  
  observeEvent( input$smisubmit , { 
    output$smirend = renderDataTable({
      inputsmiles = parse.smiles(input$insmi, kekulise=TRUE) 
      if ( is.null(inputsmiles[[1]]))
      {
        output$ui.smi.action.one <- renderUI({
          h5(HTML("<font color = \"red\"> Invalid SMILES!</font>"), align = "center")
        })
        output$ui.smi.action.two <- renderUI({ return()  })
        output$ui.smi.action.draw.in.smi <- renderUI({ return() })
        datatable(frag.annot.merged ,  caption = HTML("For more information about the table columns see the <b>FAQ</b>" ),
              rownames= FALSE, escape = FALSE, selection = 'single' , extensions = c( 'Buttons','Responsive'), 
              options = list(pageLength = 10, autoWidth = F,      
              dom = 'Bfrtip',
              buttons = c( 'csv', I('colvis'))))
      }
      else 
      {
        output$insmiimg = renderPlot(drawsmi(inputsmiles))
        output$ui.smi.action.two <- renderUI({
          h5("Input Structure", align = "center")
        })
        output$ui.smi.action.draw.in.smi <- renderUI({
          plotOutput(
              outputId = "insmiimg",
              height = "300px",
              width = "300px"
          )
        })
        s.m = parse.smiles(input$insmi)
        t = c()
        withProgress(
          message = 'Searching Fragment Database',
          value = 0,
          {
            i = 1
            t = sapply(fragments, FUN = function(x){  
              prog <<- i/ (length(fragments) *4000) ;
              incProgress(prog, detail = paste0("Working on ", i , " of ", length(fragments)))  ;  
              i <<- i + 1 ; 
              return( rcdk::matches(x, s.m ) )} )
            setProgress(1)
          }
        )
        zz = fragments[as.numeric(as.character(which(t)))]
        frag.merge.sub <<- data.table(frag.annot.merged[frag.annot.merged$Fragment %in% zz, ])
        datatable(frag.merge.sub , caption = HTML("For more information about the table columns see the <b>FAQ</b>" ),
              rownames= FALSE, escape = FALSE, selection = 'single' , extensions = c('Responsive'), 
              options = list(pageLength = 15, order = list(4, 'desc'), autoWidth = F))
      } 
    })  
  })
  output$showselsmi = renderUI({
      validate(need(
        !is.null(input$smirend_row_last_clicked),
        "Click on rows to visualize chemical structure"
      ))   
      rwsel = input$smirend_row_last_clicked
      inputsmiles = parse.smiles(unique(as.character(frag.merge.sub[rwsel,1])), kekulise=TRUE) 
      output$insmiimgsel = renderPlot({
        par(mar = c(0, 0, 0, 0))
        d = get.depictor(
          height = 1000, 
          width = 1000,
          zoom = 6,
          style = "cow"
        )
        iii = view.image.2d(inputsmiles[[1]], depictor = d)
        plot(
          NA,
          NA,
          xlim = c(0, 1000),
          ylim = c(0, 1000),
          xaxt = 'n',
          yaxt = 'n',
          xlab = '',
          ylab = ''
        )
        rasterImage(iii, 1, 1, 1000, 1000)
      })
      output$ui.smi.action.three <- renderUI({
          h5("Selected Structure", align = "center")
      })  
      plotOutput(
          outputId = "insmiimgsel",
          height = "300px",
          width = "300px"
      )
  })
  output$datat = renderDT({
    da$PDBID = as.character(paste0("<a href=\"#\" onclick = \"$('#nglview').modal('show')\",>" , da$PDBID , "</a>"))
    names(da) = c("PDB ID", "Chain", "UniProt ID", "Gene Name", "Protein Name",  "Species",  "Conformation Probablity",  "Conformation Classification" , 
     "Mutation Status", "Structure Release Date", "Latest Version", "Resolution", "PMID", "Source", 
     "D1",	"D2", 	"h_cgvc",	"ang_NHs",	"ang_CHs",	"dist_NC",	"dist_NH", 	"dist_CH", 	"h_scvc", 	"h_norm",	"h_sc_x",	"n_psi",	"r_curv")
    da = da[,c(1:14)] 
    da = da[,c(1,5,2:4,14,6:8,10:13,9)]
    datatable(da, plugins = 'searchHighlight', 
              rownames= FALSE, escape = FALSE, filter = "top" , class = 'cell-border stripe',  
              style = 'bootstrap', 
              selection = 'single' , extensions = c( 'Buttons','FixedColumns'),
              options = list(pageLength = 15, 
                columnDefs = list(list(targets = c(9:10,13), searchable = FALSE)),
                autoWidth = T,      
                dom = 'lBftipr',
                scrollX =T, 
                fixedColumns = list(leftColumns = 2),
                lengthMenu = c(5, 10, 15 , 30, 50 , nrow(da)),
                buttons = c('csv', 'excel',  I('colvis')),
                searchHighlight = TRUE
                )
              )
  }) 
  output$liglib = renderDT({
    ligand.lib$Ligand.ID = as.character(paste0(
      #'<a href="http://www.rcsb.org/pdb/ligand/ligandsummary.do;?hetId=' ,
      '<a href=', '"http://www.rcsb.org/pdb/ligand/ligandsummary.do;?hetId=', ligand.lib$Ligand.ID , '" target="_blank">',ligand.lib$Ligand.ID ,  '</a>' ))#, 
      #'">' , 
      #"</a>"))
      names(ligand.lib) = c("PDB Ligand ID", "PDB Name", "SMILES", "CIDI Structures", "CIDO Structures", "CODI Structures", "CODO Structures", "ω CD Structures")

    datatable(ligand.lib, plugins = 'searchHighlight', 
              rownames= FALSE, escape = FALSE,  class = 'cell-border stripe',  
              style = 'bootstrap', 
              selection = 'single' , extensions = c( 'Buttons','Responsive'),
              options = list(pageLength = 15, 
                #columnDefs = list(list(targets = c(9:10,13), searchable = FALSE)),
                autoWidth = T,      
                dom = 'lBftipr',
                scrollX =T, 
                lengthMenu = c(5, 10, 15 , 30, 50 , nrow(ligand.lib)),
                buttons = c('csv', 'excel',  I('colvis')),
                searchHighlight = TRUE
                )
              )
  })    

  

  output$ligandsmirender  = renderPlot({
    validate(need(
      !is.null(input$liglib_row_last_clicked),
      "Click on rows to visualize chemical structure"
    ))
    rwsel = ligand.lib[input$liglib_row_last_clicked,] 
    ligandsmi = rwsel$SMILES
    ligandsmi = as.character(ligandsmi)
    ligandsmi  = gsub("\\[n", "\\[N", ligandsmi )
    ligandmol = parse.smiles(ligandsmi)
    par(mar = c(0, 0, 0, 0))
    d = get.depictor(
      height = 1000,
      width = 1000,
      zoom = 6,
      style = "cow"
    )
    ii = view.image.2d(ligandmol[[1]], depictor = d)
    plot(
      NA,
      NA,
      xlim = c(0, 1000),
      ylim = c(0, 1000),
      xaxt = 'n',
      yaxt = 'n',
      xlab = '',
      ylab = ''
    )
    rasterImage(ii, 1, 1, 1000, 1000)
  })
  
  output$geom = renderDT({
    #da$PDBID = as.character(paste0("<a href=\"#\" onclick = \"$('#nglview').modal('show')\",>" , da$PDBID , "</a>"))
    names(da) = c("PDB ID", "Chain", "UniProt ID", "Gene Name", "Protein Name",  "Species",  "Conformation Probablity",  "Conformation Classification" , 
     "Mutation Status", "Structure Release Date", "Latest Version", "Resolution", "PMID", "Source", 
     "D1",	"D2", 	"h_cgvc",	"ang_NHs",	"ang_CHs",	"dist_NC",	"dist_NH", 	"dist_CH", 	"h_scvc", 	"h_norm",	"h_sc_x",	"n_psi",	"r_curv")
    da = da[,c(1,5,2,7,8, 15:ncol(da))]
    datatable(da, plugins = 'searchHighlight', 
              rownames= FALSE, escape = FALSE, filter = "top" , class = 'cell-border stripe',  
              style = 'bootstrap', 
              selection = 'single' , extensions = c( 'Buttons','FixedColumns'),
              options = list(pageLength = 15, 
              autoWidth = T,      
              dom = 'lBftipr',
              scrollX =T, 
              fixedColumns = list(leftColumns = 2),
              lengthMenu = c(5, 10, 15 , 30, 50 , nrow(da)),
              buttons = c('csv', 'excel',  I('colvis')),
              searchHighlight = TRUE
              ))
  })

  output$frame <- renderUI({
      rwsel = input$datat_row_last_clicked
      pdbidsel = unique(as.character(da[rwsel,1]))
      x = list( 
        h2(pdbidsel, align = "center"),
        includeScript("www/ngl.embedded.min.js"),
        HTML("<script> if( !Detector.webgl ) Detector.addGetWebGLMessage(); </script>") ,
        HTML("<script> NGL.mainScriptFilePath = \"ngl.embedded.min.js\"; </script>") ,
        HTML("<script> var stage = new NGL.Stage( \"viewport\" ); </script>") ,
        HTML(paste0("<script> stage.loadFile( \"https://files.rcsb.org/view/" , pdbidsel , ".pdb\", { defaultRepresentation: true } ); </script>")) ,
        HTML("<script> NGL.init( ); </script>") ,
        tags$div(  id="viewport" , style="width:555px; height:320px;")
      )
    x 
  })
  observeEvent(input$kon, {
    updateTabsetPanel(session, "nav",
                      selected = "Data explorer")
  })
  observeEvent(input$kin, {
    updateNavlistPanel(session, "nav",
                       selected = HTML("<i>Kinformation</i>"))
  })
  observeEvent(input$frag, {
    updateTabsetPanel(session, "nav",
                      selected = "Fragment Libraries" )
  })
  observeEvent(input$dl, {
    updateTabsetPanel(session, "nav",
                      selected = "Data explorer")
  })
  observeEvent(input$models, {
    updateTabsetPanel(session, "nav",
                      selected = "Models using DFGModel")
  })
    
  
  enrich = read.table("./tc.9.all.fragments.enriched.csv",
                      sep = ",",
                      header = T)
  enrich$key = as.character(paste(as.character(enrich$conf), as.character(enrich$frag), sep =
                                    ">>"))
  
  #enrich$key = as.character(enrich$frag )
  enrich$lp = as.numeric(as.character(enrich$lp))
  enrich.f = enrich[enrich$pval < .05, ]
  enrich.f$lp2 = enrich.f$lp + runif(n = nrow(enrich.f),
                                     min = -.01 ,
                                     max = .01)
  enrich.f$lp2 = round(enrich.f$lp2, digits = 4)
  enrich.f$key = as.character(enrich.f$key)
  head(enrich.f$lp2)
  thres1 = .05
  thres2 = .05 / 10535
  
  # Set some colors
  plotcolor <- "#F5F1DA"
  papercolor <- "#E3DFC8"
  
  colors = c( "#f43605","#7C2600", "#CE5A28" , "#00D2F1" , "#AA96DA" , "Green")  
  
  colorst = c("#f43605","#7C2600", "#CE5A28" , "#00D2F1" , "#AA96DA" , "#7F7F7F")
  # Plot time series chart
  output$timeseries <- renderPlotly({
    p = ggplot(enrich.f,
               aes(
                 x = factor(conf),
                 y = lp2,
                 key = key,
                 color = factor(conf),
                 text = frag
               ))
    p = p +
      geom_jitter(width = 0.2,
                  height = 0,
                  size = 3) +
      theme_bw() +
      scale_color_manual(values = colorst) +
      scale_y_continuous("-log(p-value)") +
      labs(y = "-log(p-value)") +
      scale_x_discrete(
        "",
        labels = c(
          "cidi" = "cHelix-in\nDFG-in\n(type 1)" ,
          "cido" = "cHelix-in\nDFG-out\n(type 2)" ,
          "codi" = "cHelix-out\nDFG-in\n(type 1/2)" ,
          "codo" = "cHelix-out\nDFG-out",
          "omega" = "omega"
        )
      ) +
      geom_hline(yintercept = -log10(.05),
                 color = "red",
                 size = 1.5)  +
      geom_hline(
        yintercept = -log10(thres2),
        color = "purple",
        size = 1.5
      )  +
      theme(
        panel.grid.major.x = element_blank() ,
        panel.grid.major.y = element_line(color = "grey", size = 1.25),
        panel.grid.minor = element_blank(),
        axis.text.y = element_text(size = 10),
        axis.text.x = element_text(
          size = 13,
          angle = 65,
          hjust = -1,
          face = "bold",
          color = "black"
        ) ,
        axis.title.y = element_text(size = 10),
        legend.position = "none"
      ) + scale_fill_discrete(guide = FALSE)
    ggplotly(p, key = enrich.f$key)#, source = "source" )
    
  })
  
  # Coupled hover event
  i = NULL
  output$correlation <- renderPlot({
    # Read in hover data
    eventdata <- event_data("plotly_click")#, source = "source")
    validate(need(
      !is.null(eventdata),
      "Click on the points to vizualize the enriched substructures"
    ))
    sel = unlist(strsplit(eventdata$key, ">>"))
    lp = sel[2]
    rw = enrich.f[enrich.f$frag == lp, ]
    #rw
    smi = as.character(rw$frag)
    smi = gsub("\\[n", "\\[N", smi)
    mol = parse.smiles(smi)
    par(mar = c(0, 0, 0, 0))
    d = get.depictor(
      height = 1000,
      width = 1000,
      zoom = 6,
      style = "cow"
    )
    i = view.image.2d(mol[[1]], depictor = d)
    plot(
      NA,
      NA,
      xlim = c(0, 1000),
      ylim = c(0, 1000),
      xaxt = 'n',
      yaxt = 'n',
      xlab = '',
      ylab = ''
    )
    rasterImage(i, 1, 1, 1000, 1000)
    
  })
  output$cor <- DT::renderDataTable(DT::datatable({
    event =   event_data("plotly_click")
    if (is.null(event))
    {
      #print("Click on a point to return the SMILEs string of the substructure")
    }
    else
    {
      sel = unlist(strsplit(event$key, ">>"))
      all = frag.annot[which(frag.annot$smi == sel[2]), ]
      all = all[all$conf == sel[1],]
      all.pdb = unlist(strsplit(as.character(all$pdb), '\\|'))
      lig.pdb = ligand.pdb[which(ligand.pdb$V1 %in% all.pdb),]
      ligs = unlist(strsplit(as.character(all$lig), "-"))
      lig.pdb = lig.pdb[which(lig.pdb$V2 %in% ligs), ]
      lig.pdb$V3 = NULL
      lig.pdb$pdb = gsub("_[A-Z]", "", lig.pdb$V1)
      lig.pdb$V1 = NULL
      lig.pdb$linklig = paste0(
        '<a href="http://www.rcsb.org/pdb/ligand/ligandsummary.do;?hetId=' ,
        lig.pdb$V2,
        '">',
        as.character(lig.pdb$V2),
        '</a>'
      )
      lig.pdb$link = paste0(
        '<a href="https://www.rcsb.org/3d-view/',
        as.character(lig.pdb$pdb),
        '?preset=ligandInteraction&sele=',
        as.character(lig.pdb$V2),
        '">',
        as.character(lig.pdb$pdb),
        '</a>'
      )
      lig.pdb$pdb = NULL
      lig.pdb$V2 = NULL
      lig.pdb$Conformation = rep(sel[1], nrow(lig.pdb))
      names(lig.pdb) = c(
        "smiles",
        "Link to parent ligand for substructure",
        "Link to structure interaction",
        "Structure conformation"
      )
      lig.pdb$smiles = NULL
      DT = setDT(lig.pdb)
      DT
      
    }
  },
  escape = FALSE))
  observeEvent(input$featureInput4, {
    if (input$featureInput4 == "all")
    {
      data3 = data
      row.names(data3) = data3$PDBid
      data3$opa = rep(1,nrow(data3))
    }
    else
    {
      data3 = data[which(data$source == input$featureInput4) ,]
      data3$opa = rep(1,nrow(data3))
    }
    observeEvent(input$featureInputgenename, {
      if ( input$featureInputgenename == "all")
      {
        data3$opa = rep(1,nrow(data3))
      }
      else 
      {
        data3$opa = rep(.1, nrow(data3))
        data3$opa = as.numeric(data3$opa)
        sel = data3[which(data3[,39] == input$featureInputgenename ),] 
        sel$opa = as.numeric(rep(.9, nrow(sel))) 
        sel$V19 = rep("Selected", nrow(sel))
        #data3non = data3[which(data3$Protein.Name != input$featureInputgenename ),] 
        #data3non$opa = rep(.1,nrow(data3non))
        data3 = rbind(data3, sel)

        data3$V19 = as.factor(as.character(data3$V19))
      }
      output$scaplott <- renderPlotly({
      key = data3$PDBid
      plot_ly(
        data3,
        x = data3[, input$featureInput1] ,
        y = data3[, input$featureInput2] ,
        z = data3[, input$featureInput3],
        color = ~ V19,
        key = ~ key ,
        colors = colors,
        opacity = data3$opa , 
        mode = 'markers',
        hoverinfo = 'text',
        text = ~ paste('Id:', PDBid, '<b>', V19, '<br>', Protein.Name)
      ) %>%
        layout(
          title = paste(
            input$featureInput1,
            "vs",
            input$featureInput2,
            "vs",
            input$featureInput3,
            sep = " "
          ),
          scene = list(
            xaxis = list(title = input$featureInput1),
            yaxis = list(title = input$featureInput2),
            zaxis = list(title = input$featureInput3)
          )
        )
      })
    })

    output$click <- renderText({
      d <- event_data("plotly_click")
      if (is.null(d))
      {
        "Click points on the plot to get a link to the PDB"
      }
      else
      {
        pd = gsub("_[A-Z]$", "",  d$key)
        lin = paste("https://www.rcsb.org/3d-view/", pd, "/1" , sep = "")
        mod = paste0("<a href='",  lin, "' target='_blank'>View PDB</a>")
        as.character(paste(d$key, mod, sep = ": "))
        #browseURL(lin)
      }
    })
  })
  observeEvent(input$file1, {
    z = try(read.pdb(input$file1$datapath) , silent = T, TRUE)
    if (!inherits(z, "try-error"))
    {
      confirmSweetAlert(
        session = session,
        title = "Looks like a valid PDB file",
        text = "Do you wish to continue running Kinformation?
        It may take some time to run.",
        type = "success",
        danger_mode = T,
        inputId = "confirm" ,
        closeOnClickOutside = T ,
        btn_labels = c("No", "Yes!")
      )
    }
    else
    {
      sendSweetAlert(
        session = session,
        title = "Error!",
        text = "This file does not look like a valid PDB file",
        type = "error" ,
        btn_labels = c("Okay") ,
        closeOnClickOutside = T
      )
    }
  })
  observeEvent(req(input$confirm),  {
    if (req(input$confirm))
    {
      withProgress(
        message = 'Running Kinformation',
        detail = "Creating Files",
        value = 0.1,
        {
          pdb = read.pdb(input$file1$datapath)
          Sys.sleep(1)
          rand = sample(1:100000, 1)
          output$jobid <- renderText({
            paste0("JobID: ", rand)
            
          })
          fn = paste0("/srv/shiny-server/4_Konformation/", rand, ".pdb")
          incProgress(0.1, detail = "Preparing Data")
          file.copy(input$file1$datapath, fn)
          fn = paste0(rand, ".pdb")
          cmd = paste(
            "/srv/shiny-server/4_Konformation/prep_and_run_konformation.sh",
            fn,
            sep = " "
          )
          Sys.sleep(1)
          incProgress(0.01, detail = "Generating Descriptors")
          o = system(cmd, intern = T)
          Sys.sleep(1)
          incProgress(0.01, detail = "Running Kinformation")
          desc = read.table(
            paste0(
              "/srv/shiny-server/4_Konformation/1_result/",
              rand,
              ".csv"
            ),
            sep = ",",
            header = T,
            stringsAsFactors = F
          )
          pred = predict_conformation(desc)
          data2 = data[, c(1, 20:29)]
          protnam = as.character(data[,39])
          data2$dist_NC = NULL
          row.names(data2) = data2$PDBid
          #data2$PDBid = NULL
          row.names(desc) = desc$pdb_id
          data2$conformation = data2$V19
          data2$V19 = NULL
          desc = desc[2, c(3, 4, 5, 7, 8, 9, 10, 6)]
          desc$PDBid = "Input Structure"
          desc$conformation = "Added Input Structure"
          desc$opa = 1 
          data2$opa = rep(.3,  nrow(data2))
          total = rbind(desc, data2)
          similarity = data.frame()
          numprog = .1 / nrow(data2)
          if ( input$pairwise )
          {
            incProgress(0.1, detail = "Comparing Descriptors")
            for (inum in 1:nrow(data2))
            {
              rw = data2[inum,]
              prot = protnam[inum]
              combos=as.data.frame(rbind(rw,desc))
              combos[,c("PDBid", "conformation")] = NULL 
              dis = as.numeric(dist(data.matrix(combos)))
              rw2 = as.data.frame(rbind(c(as.character(unique(rw$PDBid)), round(dis, digits = 3), prot ,  as.character(unique(rw$conformation)))))
              names(rw2)= c("PDBid", "EuclideanDistance", "Protein Name", "Conformation")
              similarity = as.data.frame(rbind(similarity, rw2))
              similarity$EuclideanDistance = as.numeric(as.character(similarity$EuclideanDistance)) 
              incProgress(numprog , detail = "Comparing Descriptors")
            }
          }
          
          setProgress(1)
        }
      )
      output$kclass = renderImage({
        kclass = names(which.max(apply(pred,2, max)))
        if ( kclass == "cidi" )
        {
          list(src = 'www/cidi.png',
               contentType = 'image/png',
               width = 610,
               height = 320,
               align = "center",
               alt = "chelix-in dfg-in")
        }
        else if (kclass == "cido" ) {
          list(src = 'www/cido.png',
               contentType = 'image/png',
               width = 610,
               height = 320,
               align = "center",
               alt = "chelix-in dfg-out")
        }
        else if (kclass == "codo" ) {
          list(src = 'www/codo.png',
               contentType = 'image/png',
               width = 610,
               height = 320,
               align = "center",
               alt = "chelix-out dfg-out")
        }
        else if (kclass == "codi" ) {
          list(src = 'www/codi.png',
               contentType = 'image/png',
               width = 610,
               height = 320,
               align = "center",
               alt = "chelix-out dfg-in")
        }
        else {
          list(src = 'www/omega.png',
               contentType = 'image/png',
               width = 610,
               height = 320,
               
               align = "center",
               alt = "omegaCD")
        }
      }, deleteFile = FALSE)
      output$contents = renderDT({
        names(pred) = c(
          "cHelix in/DFG in",
          "cHelix in/DFG out",
          "cHelix out/DFG in",
          "cHelix out/DFG out",
          "ωCD"
        )
        datatable(
          t(pred),
          caption = 'Table of Probablities for each conformation',
          options = list(
            autoWidth = T,
            paging = FALSE,
            searching = FALSE
          )
        )
      })
      output$downloadData <- downloadHandler(
        filename = function() {"probablities.csv"} ,
        content = function(file) { write.csv(t(pred), file, row.names = TRUE, sep = ",") } 
      )
      output$downloadDesc <- downloadHandler(
        filename = function() {"descriptors.csv"} ,
        content = function(file) { write.csv(desc, file, row.names = TRUE, sep = ",") } 
      )
      if ( input$pairwise )
      {
        output$downloadComp <- downloadHandler(
          filename = function() {"pairwisecomparisons.csv"} ,
          content = function(file) { write.csv(similarity, file, row.names = TRUE, sep = ",") } 
        )
        output$ui.actionfive <- renderUI({
          if (is.null(input$pairwise)) return()
          downloadButton('downloadComp', HTML('Download euclidean distances between <br> uploaded structure and available </br>kinase structures'))
        })
      }
      
      output$ui.action <- renderUI({
        if (is.null(pred)) return()
        downloadButton('downloadData', 'Download probablities of conformations')
      })
      output$ui.actiontwo <- renderUI({
        if (is.null(pred)) return()
        downloadButton('downloadDesc', HTML('Download geometric descriptors used by<br>Kinformation</br>'))
      })
      output$ui.actionthree <- renderUI({
        if (is.null(pred)) return()
        h5(HTML("For more information about the descriptors used, please refer to the <b>FAQ</b> or our Cell Chemical Biology article"))
      })
         
      output$scaplot <- renderPlotly({
        plot_ly(
          total,
          
          x = total$p1p1x,
          y = total$p2p2x,
          z = total$h_cgvc,
          opacity = ~ opa ,
          group = ~ factor(conformation) ,
          color = ~ conformation,
          colors = colors,
          text = ~ paste('Id:', PDBid)
        )  %>%
          layout(title = "Comparison of uploaded structure to other kinase structures", scene = list(
            xaxis = list(title = "DFG Vector 1"),
            yaxis = list(title = "DFG Vector 2"),
            zaxis = list(title = "cHelix Vector"),
            camera = list(eye = list(x = desc$p1p1x , y = desc$p2p2x, z = desc$h_cgvc ))
          ))
      })
      if ( input$pairwise )
      {
        output$similar <- renderDT({
          names(similarity) = c("PDB ID",  "Euclidean Distance" , "Protein Name",  "Conformation")
          DTab = datatable(similarity, options = list(order = list(list(2, 'asc'))), rownames= FALSE)
          DTab
        })
      }
      
    }
  })
}


shinyApp(ui, server)

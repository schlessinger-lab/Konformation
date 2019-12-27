models = read.table("./available.models.csv", sep = ",", header = T)

dfgui = function()
{
    dfguicode = list(
        h3("Models of Kinases in Unseen Conformations ", align = "center"),
        tags$br(),
        h4(HTML("These Models are available for <b>academic</b> and <b>educational</b> use under a Academic Public License and is intended for non-commercial use."), align = "center"),
        h4("For a commercial license please contact the Schlessinger laboratory.", align = "center"),
        h4("To download models, select a row from the table below and click the 'Download selected models' button below", align="center"),
        downloadButton("downloadmodels", "Download selected models"),
        tags$br(), 
        dataTableOutput("models")
        # fluidRow(
        #     column(
        #         8,
        #         dataTableOutput("models")
        #     ),
        #     column(
        #         4,
        #         h3("hello")
        #     )
        # )
    )
    return(dfguicode) 
}

dfgserver = function()
{
    
    names(models) = c("Uniprot ID", "Gene Name", "Modeled Conformation", "Download")
    #models$Download = rep(downloadLink("test",label = "download"), nrow(models))
    print(panda)
    print(length( rep(downloadLink("test",label = "download"), nrow(models))))
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
}
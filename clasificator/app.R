library(shiny)
library(caret)
library(poppr)
library(nnet)

newPredict <- function(newdata) {
    #Load the population to be classified
    newpop <- read.genalex(newdata,ploidy =3)
    
    #Load the trained model 
    trainedModel <- readRDS("Trained_model.rds")
    
    #Load the allele table 
    newtable <- newpop@tab
    
    #Remove uninformative alleles
    uninf <- grep("\\.0$",colnames(newtable))
    newtable <- newtable[,-uninf]
    
    #Equate the predictors in the model to the variables in the new dataset
    missingAl <- which(!trainedModel$finalModel$xNames %in% colnames(newtable))
    tempMat <- matrix(data = 0, nrow = length(newtable[,1]), ncol = length(missingAl), dimnames = list(rownames(newtable), trainedModel$finalModel$xNames[missingAl]))
    newtable <- cbind(newtable,tempMat)
    
    #Generate the prediction and arrange in a results table
    resTable <- data.frame(Id = rownames(newtable), Lineage = predict(trainedModel,newtable), Probability = unlist(apply(predict(trainedModel,newtable,type = "prob"),MARGIN = 1,FUN = max)))
    
    #Return the resulting table
    return(resTable)
    
}

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Phytophthora infestans lineage clasificator"),
    sidebarLayout(
        sidebarPanel(
            fileInput("file", "Genalex File", buttonLabel = "Upload", accept = ".csv"),
            textOutput("accuracy"),
            downloadButton("report", "Download results"),
            width = 3
        ),
        mainPanel(
            tableOutput("prediction")
        )
    )
)

# Define server logic required to draw a histogram
server <- function(input, output) {
    trainedModel <- readRDS("Trained_model.rds")
    accMl <- round(mean(trainedModel$resample$Accuracy) * 100, digits = 4)
    output$accuracy <- renderText(paste("Training accuracy is ",accMl, "%",sep = ""))
    
    popdata <- reactive({
        req(input$file)
        
        input$file$datapath
    })

    output$prediction <- renderTable({
        newPredict(popdata())
    })
    
    output$report <- downloadHandler(
        filename = "Predicted_lineages.csv",
        content = function(file) {
            write.csv(newPredict(popdata())[,-1],file)
        }
    )
    
}

# Run the application 
shinyApp(ui = ui, server = server)

# this shinyapp may work in tensorflow 2.4.0
library(shiny)
library(keras)
library(tfhub)

# Load the model,define where your models are 
model <- load_model_tf("/home/zyh/Applications/rstudio/mymodel/")

# Define the UI
ui <- fluidPage(
  # App title ----
  titlePanel("Hello TensorFlow!"),
  # Sidebar layout with input and output definitions ----
  sidebarLayout(
    # Sidebar panel for inputs ----
    sidebarPanel(
      # Input: File upload
      fileInput("image_path", label = "Input a JPEG image")
    ),
    # Main panel for displaying outputs ----
    mainPanel(
      # Output: Histogram ----
      textOutput(outputId = "prediction"),
      plotOutput(outputId = "image")
    )
  )
)

# Define server logic required to draw a histogram ----
server <- function(input, output) {
  
  image <- reactive({
    req(input$image_path)
    jpeg::readJPEG(input$image_path$datapath)
  })
  
  output$prediction <- renderText({
    
    img <- image() %>% 
      array_reshape(., dim = c(1, dim(.), 1)) 
      
    #predict_class <- classifier(tf$constant(img, tf$float32))
    #result <- np.squeeze(model.predict(img))
    #result <- tf.keras.layers.Softmax()(result)
    #predict_class <- tf$argmax(img, axis = 1L) %>% as.integer()
    #paste0("The predicted class number is ", predict(model, img))
    paste0("The predicted class number is ", predict_classes(model, img))
  })
  
  output$image <- renderPlot({
    plot(as.raster(image()))
  })
  
}

shinyApp(ui, server)
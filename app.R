library(shiny)
library(shinythemes)
library(readxl)
library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(ggtext)
library(RColorBrewer)
library(ggsci)
library(emmeans)


options(shiny.maxRequestSize = 200 * 1024^2)  # 200 MB

#### FUNCTIONS ####


# Defining the exponential model function based on the Nature 2015 paper methods section
model_function <- function(t, dt, A0, A1, tauA, tauB) {
  tau <- t - dt
  H_tau <- ifelse(tau > 0, 1, 0)
  
  # Avoid division by zero or negative time constants
  tauA <- max(tauA, .Machine$double.eps)
  tauB <- max(tauB, .Machine$double.eps)
  
  exp_neg_tauA <- exp(-max(tau, 0) / tauA)
  exp_neg_tauB <- exp(-max(tau, 0) / tauB)
  
  (1 - exp_neg_tauA) * (A0 + A1 * exp_neg_tauB) * H_tau
}

# Function to calculate an R-squared value of the exponential fit.
calculate_R_squared <- function(params, data) {
  names(params) <- c("dt", "A0", "A1", "tauA", "tauB")
  
  predicted <- sapply(data$time, function(t) {
    model_function(t, params["dt"], params["A0"], params["A1"], params["tauA"], params["tauB"])
  })
  
  # Check for NaN or Inf in predicted
  if(any(is.na(predicted) | is.infinite(predicted))) {
    cat("Non-finite values detected in predicted speeds.\n")
    return(Inf)  # Return Inf to indicate a bad result
  }
  
  # Checking how close the predicted speed (based on the model function) is to the actual calculated relative speed. 
  residuals <- data$relative_speed - predicted
  SSE <- sum(residuals^2)
  SST <- sum((data$relative_speed - mean(data$relative_speed))^2)
  
  R_squared <- 1 - SSE / SST
  
  if(is.finite(R_squared)){
    return(-R_squared)
  } else {
    cat("Non-finite R_squared value detected.\n")
    return(0)
  }
}


#setwd("/Volumes/tuxworri-ddr-neuropathology/Isabelle Clowes/tracking")
#experiments <- excel_sheets("./Habituation_data_collated.xlsx")

#### UI ####


ui <- fluidPage(
  theme = shinytheme("simplex"),
  titlePanel("DART data analysis and curve fitting"),
  
  sidebarLayout(
    sidebarPanel(
      width = 3, 
      tags$h4("Calculation options"),
      numericInput("dtflexibility", "dt flexibility:", min = 0, max = 10, step = 0.5, value = 2),
      tags$hr(),
      
      tags$h4("Plotting options"),
      numericInput("ylimit", "y-axis upper limit:", min = 1, max = 20, step = 1, value = 5),
      numericInput("ylowlimit", "y-axis lower limit:", min = -10, max = 5, step = 1, value = -1), 
      numericInput("axislinewidth", "Axis line width:", min = 0.1, max = 5, step = 0.1, value = 1),
      numericInput("axistitlefontsize", "Axis title font size:", min = 0, max = 60, step = 1, value = 22),
      numericInput("axistextfontsize", "Axis text font size:", min = 0, max = 60, step = 1, value = 20),
      numericInput("plottitlefontsize", "Plot title font size:", min = 0, max = 60, step = 1, value = 0),
      numericInput("legendtextfontsize", "Legend text font size:", min = 0, max = 60, step = 1, value = 18),
      numericInput("legendtitlefontsize", "Legend title font size:", min = 0, max = 60, step = 1, value = 18),
      numericInput("keysize", "Legend key size:", min = 0, max = 20, step = 0.5, value = 1.5),
      numericInput("keywidth", "Legend key width:", min = 0, max = 20, step = 0.5, value = 4),
      numericInput("rawlinewidth", "Raw speed line thickness:", min = 0, max = 10, step = 0.1, value = 0.5),
      numericInput("plotlinewidth", "Relative speed line thickness:", min = 0, max = 10, step = 0.1, value = 0.5),
      numericInput("fitlinewidth", "Fit line thickness:", min = 0, max = 10, step = 0.1, value = 1.5),
      tags$hr(),
      tags$h4("Linear regression (peaks)"),
      checkboxInput("lm_enable", "Overlay linear regression across peak responses", value = FALSE),
      selectInput(
        "lm_y",
        "Peak series to regress",
        choices = c(
          "Fitted (mean predicted_speeds)" = "mean_speed",
          "Relative (mean speed - pre-stim baseline)" = "mean_rel_speed",
          "Raw (mean speed)" = "mean_raw_speed"
        ),
        selected = "mean_speed"
      ),
      checkboxInput("lm_show_ci", "Show 95% CI ribbon", value = TRUE),
      numericInput("lm_ci_mult", "CI multiplier (1.96 ≈ 95%)", min = 0, max = 10, step = 0.01, value = 1.96),
      
      
      tags$hr(),
      downloadButton("download_data", "Download full data")
    ),
    
    mainPanel(
      width = 9,
      tabsetPanel(
        
        tabPanel(
          "Import data",
          
          fluidRow(
            h2("File input"),
            column(
              6,
              fileInput(
                "file",
                "Choose Excel file (.xlsx / .xls)",
                multiple = FALSE,
                accept = c(".xlsx", ".xls")
              ),
              verbatimTextOutput("inputfile")
            )
          ),
          
          fluidRow(
            h2("Selections"),
            column(3, selectInput("data", "Data", choices = NULL)),
            column(3, selectInput("group", "Select Group (only 1 for calculations)", choices = c("test"), multiple = TRUE)),
            column(3, numericInput("numberstim", "How many stimulations in total?", value = 6, min = 1, step = 1)),
            column(3, selectInput("stimulations", "Select stimulation number(s)", choices = NULL, multiple = TRUE))
          ),
          
          fluidRow(
            h2("Set defaults"),
            column(6, actionButton("resethabit", "Set habituation default")),
            column(6, actionButton("resetramp", "Set ramp default"))
          ),
          
          fluidRow(
            h2("Parameter settings"),
            column(
              3,
              numericInput(
                "pretime",
                "Time window before stimulation (minutes)",
                min = 0, max = 30, step = 1, value = 2
              )
            ),
            column(
              3,
              numericInput(
                "prestim",
                "Pre-stimulation time window for calculations (seconds)",
                min = 5, max = 300, step = 10, value = 60
              )
            ),
            column(
              3,
              numericInput(
                "stimgap",
                "Time gap between stimulations (seconds)",
                min = 10, max = 3600, step = 10, value = 150
              )
            ),
            column(
              3,
              numericInput(
                "poststim",
                "Post-stimulation time window for calculations (seconds)",
                min = 10, max = 300, step = 10, value = 120
              )
            )
          ),
          
          fluidRow(
            column(6, tableOutput("output")),
            column(4, tableOutput("gofoutput")),
            column(
              2,
              actionButton(
                "cleartable",
                HTML("Clear tables<br>(do this between groups<br>or when changing stims)")
              )
            )
          ),
          
          fluidRow(
            column(12, actionButton("calcparam", "Calculate/update parameters"))
          ),
          
          fluidRow(
            h2("Per-fly curve calculation"),
            column(3, actionButton("sequentialcalc", "Calculate per-fly curves"))
          ),
          
          fluidRow(
            h2("Random fly selection"),
            column(3, numericInput("randomno", "How many flies to randomly sample?", value = 10, min = 1, step = 1, max = 20)),
            column(3, numericInput("times", "How many times to repeat calculations?", value = 10, min = 1, step = 1, max = 20)),
            column(3, actionButton("randomcalc", "Select x flies randomly and calculate parameters"))
          )
        ),
        
        tabPanel(
          "Data table",
          fluidRow(tableOutput("datatable"))
        ),
        
        tabPanel(
          "Raw speed plot",
          fluidRow(column(4, actionButton("plotfull", "Update full data graph"))),
          fluidRow(plotOutput("fulldata", width = "90%"))
        ),
        
        tabPanel(
          "Full fitted plot",
          fluidRow(column(4, actionButton("plotfullfitted", "Update full fitted data graph"))),
          fluidRow(plotOutput("fulldatafitted", width = "1200px", height = "600px"))
        ),
        
        tabPanel(
          "Regression (peaks)",
          fluidRow(
            column(12, tableOutput("peak_table"))
          ),
          fluidRow(
            column(12, verbatimTextOutput("lm_summary"))
          )
        )
        
        
      )
    )
  )
)

##### SERVER #####

server <- function(input, output, session) {
  
  # Print info about the uploaded file
  output$inputfile <- renderPrint({
    req(input$file)
    list(
      name = input$file$name,
      size_bytes = input$file$size,
      type = input$file$type
    )
  })
  
  # Path to the uploaded temp file on the server
  selected_file <- reactive({
    req(input$file)
    input$file$datapath
  })
  
  # Update sheet choices when a new file is uploaded
  observeEvent(selected_file(), {
    file_path <- selected_file()
    
    # Optional: basic extension check
    ext <- tolower(tools::file_ext(input$file$name))
    validate(need(ext %in% c("xlsx", "xls"), "Please upload an .xlsx or .xls file."))
    
    experiments <- excel_sheets(file_path)
    updateSelectInput(session, "data", choices = experiments, selected = experiments[1])
  }, ignoreNULL = TRUE)
  
  # Read the selected sheet
  sheet <- reactive({
    req(input$data)
    file_path <- selected_file()
    req(file_path)
    
    read_excel(file_path, sheet = input$data)
  })
  # Reactive expression to track the number of stimulations
  numberstim <- reactive({
    req(input$numberstim)
    input$numberstim
  })
  
  # Observe the number of stimulations and update the selectInput choices
  observe({
    stim_choices <- 1:numberstim()
    updateSelectInput(session, "stimulations", choices = stim_choices)
  })
  
  # Update the "control" and "groups" selectInput choices based on the unique 
  # "Group" names in the data
  observe({
    data <- sheet()
    data <- data[,-1]
    unique_groups <- as.factor(colnames(data))
    updateSelectInput(session, "group", choices = unique_groups)
  })
  
  # Reset choices to defaults for the 2 experiment types (habituation and ramp)
  
  observeEvent(input$resethabit, {
    updateNumericInput(inputId = "numberstim", value = "6")
    updateNumericInput(inputId = "pretime", value = "2")
    updateNumericInput(inputId = "prestim", value = "60")
    updateNumericInput(inputId = "stimgap", value = "150")
    updateNumericInput(inputId = "poststim", value = "120")
  })
  
  observeEvent(input$resetramp, {
    updateNumericInput(inputId = "numberstim", value = "5")
    updateNumericInput(inputId = "pretime", value = "2")
    updateNumericInput(inputId = "prestim", value = "120")
    updateNumericInput(inputId = "stimgap", value = "300")
    updateNumericInput(inputId = "poststim", value = "300")
  })
  
  # Function to get data based on stimulation number
  get_data <- function(data, stim_number, pretime, prestim, stimgap, poststim) {
    stim_start <- (stim_number - 1) * stimgap + (pretime * 60)
    stim_end <- stim_start + poststim
    data[(data$time >= stim_start-prestim) & (data$time <= stim_end), ]
  }
  
  # Initialize reactive values
  values <- reactiveValues(
    data = NULL,
    merged_data = NULL, ## this is to make a data frame with all of the predicted values
    optimized_params = NULL,
    pre_stim_speed = NULL,
    pre_stim_speed_list = NULL,
    stim_number = NULL,
    stim_number_list = NULL,
    optim_results = NULL,
    max_time = NULL,
    max_speed = NULL,
    max_time_rel = NULL,
    max_speed_rel = NULL,
    max_speed_list = NULL,
    max_rel_speed_list = NULL,
    GOFs = NULL
  )
  
  # ---------------------------
  # Linear regression on peaks
  # ---------------------------
  
  stimulus_locations <- reactive({
    req(input$numberstim, input$stimgap, input$pretime)
    vapply(
      X = seq_len(as.numeric(input$numberstim)),
      FUN = function(i) (i - 1) * as.numeric(input$stimgap) + (as.numeric(input$pretime) * 60),
      FUN.VALUE = numeric(1)
    )
  })
  
  summarized_data <- reactive({
    req(values$merged_data)
    values$merged_data %>%
      group_by(time, StimNo) %>%
      summarize(
        mean_speed     = mean(predicted_speeds, na.rm = TRUE),
        mean_raw_speed = mean(speed, na.rm = TRUE),
        se_speed       = sd(predicted_speeds, na.rm = TRUE) / sqrt(n()),
        mean_rel_speed = mean(relative_speed, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  lm_peak_fit <- reactive({
    if (!isTRUE(input$lm_enable)) return(NULL)
    
    df <- summarized_data()
    y_col <- req(input$lm_y)
    
    peaks <- df %>%
      group_by(StimNo) %>%
      slice_max(.data[[y_col]], n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        stim_no   = as.numeric(StimNo),
        stim_no_c = as.numeric(StimNo) - 1,
        stim_time = (as.numeric(StimNo) - 1) * as.numeric(input$stimgap) + (as.numeric(input$pretime) * 60),
        time_peak = time,
        peak      = .data[[y_col]]
      ) %>%
      arrange(stim_no)
    
    if (nrow(peaks) < 2 || all(is.na(peaks$peak))) return(NULL)
    
    fit <- lm(peak ~ stim_no_c, data = peaks)
    
    pred <- peaks %>%
      mutate(
        pred = predict(fit, newdata = ., se.fit = TRUE)$fit,
        se   = predict(fit, newdata = ., se.fit = TRUE)$se.fit,
        lo   = pred - as.numeric(input$lm_ci_mult) * se,
        hi   = pred + as.numeric(input$lm_ci_mult) * se
      )
    
    list(peaks = peaks, pred = pred, fit = fit)
  })
  
  output$peak_table <- renderTable({
    x <- lm_peak_fit()
    validate(need(!is.null(x), "Enable the regression overlay and generate curves/plots to see results here."))
    x$peaks
  })
  
  output$lm_summary <- renderPrint({
    x <- lm_peak_fit()
    validate(need(!is.null(x), "Enable the regression overlay and generate curves/plots to see results here."))
    summary(x$fit)
  })
  
  
  # Code to calculate the parameters for the fitted line upon pressing the
  # calculate parameters button
  observeEvent(input$calcparam, {
    
    
    for(i in as.numeric(input$stimulations)) {
      values$stim_number <- i
      showNotification("Retrieving data for stimulation number:", i, type = "message")
      cat("Retrieving data for stimulation number:", i, "\n")
      
      data <- sheet()
      data$time <- as.numeric(data$time)
      stim_data <- get_data(data, i, as.numeric(input$pretime), as.numeric(input$prestim), as.numeric(input$stimgap), as.numeric(input$poststim))
      req(stim_data)  # Ensure data is not NULL before continuing
      
      showNotification("Data retrieved successfully.", type = "message")
      cat("Data retrieved successfully.\n")
      
      # Data needs to be tidied into long format
      stim_data <- stim_data %>% 
        pivot_longer(!time, names_to = "Group", values_to = "speed") %>% 
        filter(Group %in% input$group)
      # 60 second time window prior to stimulation taken to calculate the pre-stimulation speed
      pre_stimuli_avg_speed <- mean(stim_data$speed[1:as.numeric(input$prestim)])
      # Mean of pre-stimulation speed subtracted from raw speed to give relative speed. Necessary for
      # curve fitting.
      stim_data$relative_speed <- stim_data$speed - pre_stimuli_avg_speed
      
      showNotification("Optimising parameters...", type = "message")
      cat("Optimising parameters...\n")
      # Optimizing parameters of the fitted exponential curve.
      # Initial parameter estimates
      initial_params <- c(dt = input$prestim, A0 = 0.1, A1 = 2, tauA = 0.1, tauB = 3)
      
      optim_results <- optim(
        par = initial_params,
        fn = calculate_R_squared,
        data = stim_data,
        method = "L-BFGS-B",
        lower = c(dt = stim_data$time[isolate(input$prestim)-isolate(input$dtflexibility)], A0 = 0, A1 = 0, tauA = -10, tauB = 0),
        upper = c(dt = stim_data$time[isolate(input$prestim)+isolate(input$dtflexibility)], A0 = 10, A1 = 20, tauA = 20, tauB = 50),
        control = list(maxit = 9999999)
      )
      
      optimized_params <- optim_results$par
      names(optimized_params) <- c("dt", "A0", "A1", "tauA", "tauB")
      
      # Making a new column in data which contains the predicted speeds - i.e, the speeds calculated
      # through the exponential model.
      stim_data$predicted_speeds <- sapply(stim_data$time, function(t) model_function(t, optimized_params["dt"], optimized_params["A0"], optimized_params["A1"], optimized_params["tauA"], optimized_params["tauB"]))
      stim_data$GOF <- -optim_results$value
      stim_data$StimNo <- values$stim_number
      # Find the index of the maximum speeds
      max_index <- which.max(stim_data$predicted_speeds)
      max_index_rel <- which.max(stim_data$relative_speed)
      
      # Find the time and value of the maximum speeds
      max_time <- stim_data$time[max_index]
      max_speed <- stim_data$predicted_speeds[max_index]
      max_time_rel <- stim_data$time[max_index_rel]
      max_speed_rel <- stim_data$relative_speed[max_index_rel]
      
      # Storing results in reactive values
      values$data <- stim_data
      values$pre_stim_speed <- pre_stimuli_avg_speed
      values$pre_stim_speed_list <- c(values$pre_stim_speed_list, pre_stimuli_avg_speed)
      values$merged_data <- rbind(values$merged_data, values$data)
      values$optimized_params <- optim_results$par
      values$optim_results <- optim_results
      values$max_time <- max_time
      values$max_speed <- max_speed
      values$max_speed_list <- c(values$max_speed_list, max_speed)
      values$max_rel_speed_list <- c(values$max_rel_speed_list, max_speed_rel)
      values$stim_number_list <- c(values$stim_number_list, values$stim_number)
      values$max_time_rel <- max_time_rel
      values$max_speed_rel <- max_speed_rel
      values$GOFs <- c(values$GOFs, stim_data$GOF[1])
      
      # Outputting the top of the data table with the newly calculated 
      output$output <- renderTable({
        head(values$data)
      })
      
      output$gofoutput <- renderTable({
        collated_data <- cbind(values$stim_number_list, values$pre_stim_speed_list, values$max_rel_speed_list, values$max_speed_list, values$GOFs)
        colnames(collated_data) <- c("Stim number", "Mean pre-stim speed", "Max relative amplitude", "Max fitted amplitude", "GOF")
        collated_data
      })
      showNotification("Parameters should now be optimised.", type = "message")
      cat("Parameters should now be optimised.\n")
      
      # clear stimdata 
      
      stim_data <- c()
      #i <- i+1
    }
  })
  
  observeEvent(input$cleartable, {
    collated_data <- NULL
    values$GOFs <- NULL
    values$pre_stim_speed <- NULL
    values$pre_stim_speed_list <- NULL
    values$max_speed_list <- NULL
    values$max_rel_speed_list <- NULL
    values$stim_number_list <- NULL
    values$merged_data <- NULL
    output$gofoutput <- renderTable({
      collated_data
    })
  })
  
  observe({
    # Ensure the optim_results value is available and not NULL
    if (!is.null(values$optim_results) && -values$optim_results$value < 0.2) {
      showNotification("GOF is less than 0.2", type = "message")
      cat("GOF is less than 0.2\n")
    }
  })
  
  # Code to randomly select x flies from the data, take a mean and then calculate
  # the fit. 
  observeEvent(input$randomcalc, {
    
    for(i in as.numeric(input$stimulations)) {
      for(j in 1:input$times){
        values$stim_number <- i
        showNotification("Retrieving data for stimulation number:", i, type = "message")
        cat("Retrieving data for stimulation number:", i, "\n")
        
        data <- sheet()
        data$time <- as.numeric(data$time)
        stim_data <- get_data(data, i, as.numeric(input$pretime), as.numeric(input$prestim), as.numeric(input$stimgap), as.numeric(input$poststim))
        req(stim_data)  # Ensure data is not NULL before continuing
        
        showNotification("Data retrieved successfully.", type = "message")
        cat("Data retrieved successfully.\n")
        
        
        
        # Vector of random columns for data selection
        if (ncol(stim_data) > 9) {
          randomcols <- sample(2:ncol(stim_data), isolate(input$randomno), replace = FALSE)
          
          # subset the data table including only these random columns
          stim_data <- stim_data[, c(1, randomcols)]
          #print(stim_data)
          
          # Calculate the row means (excluding the first row for each mean calculation)
          row_means <- rowMeans(stim_data[, -1], na.rm = TRUE)
          
          
          # Add the new column to the data frame
          stim_data$tenmeans <- row_means
          
          # Select only the time and tenmeans columns
          stim_data <- stim_data[, c(1, ncol(stim_data))]
        }
        else{
          stop("Not enough columns in data")
        }
        
        
        
        # Data needs to be tidied into long format
        stim_data <- stim_data %>% 
          pivot_longer(!time, names_to = "Group", values_to = "speed")
        # 60 second time window prior to stimulation taken to calculate the pre-stimulation speed
        pre_stimuli_avg_speed <- mean(stim_data$speed[1:as.numeric(input$prestim)])
        # Mean of pre-stimulation speed subtracted from raw speed to give relative speed. Necessary for
        # curve fitting.
        stim_data$relative_speed <- stim_data$speed - pre_stimuli_avg_speed
        
        showNotification("Optimising parameters...", type = "message")
        cat("Optimising parameters...\n")
        # Optimizing parameters of the fitted exponential curve.
        # Initial parameter estimates
        initial_params <- c(dt = input$prestim, A0 = 0.1, A1 = 2, tauA = 0.1, tauB = 3)
        
        optim_results <- optim(
          par = initial_params,
          fn = calculate_R_squared,
          data = stim_data,
          method = "L-BFGS-B",
          lower = c(dt = stim_data$time[isolate(input$prestim)-isolate(input$dtflexibility)], A0 = 0, A1 = 0, tauA = -10, tauB = 0),
          upper = c(dt = stim_data$time[isolate(input$prestim)+isolate(input$dtflexibility)], A0 = 10, A1 = 20, tauA = 20, tauB = 50),
          control = list(maxit = 9999999)
        )
        
        optimized_params <- optim_results$par
        names(optimized_params) <- c("dt", "A0", "A1", "tauA", "tauB")
        
        # Making a new column in data which contains the predicted speeds - i.e, the speeds calculated
        # through the exponential model.
        stim_data$predicted_speeds <- sapply(stim_data$time, function(t) model_function(t, optimized_params["dt"], optimized_params["A0"], optimized_params["A1"], optimized_params["tauA"], optimized_params["tauB"]))
        stim_data$GOF <- -optim_results$value
        stim_data$StimNo <- values$stim_number
        # Find the index of the maximum speeds
        max_index <- which.max(stim_data$predicted_speeds)
        max_index_rel <- which.max(stim_data$relative_speed)
        
        # Find the time and value of the maximum speeds
        max_time <- stim_data$time[max_index]
        max_speed <- stim_data$predicted_speeds[max_index]
        max_time_rel <- stim_data$time[max_index_rel]
        max_speed_rel <- stim_data$relative_speed[max_index_rel]
        
        # Updating all of the reactive values
        values$data <- stim_data
        values$pre_stim_speed <- pre_stimuli_avg_speed
        values$pre_stim_speed_list <- c(values$pre_stim_speed_list, pre_stimuli_avg_speed)
        values$optimized_params <- optim_results$par
        values$optim_results <- optim_results
        values$max_time <- max_time
        values$max_speed <- max_speed
        values$max_speed_list <- c(values$max_speed_list, max_speed)
        values$max_rel_speed_list <- c(values$max_rel_speed_list, max_speed_rel)
        values$stim_number_list <- c(values$stim_number_list, values$stim_number)
        values$max_time_rel <- max_time_rel
        values$max_speed_rel <- max_speed_rel
        values$GOFs <- c(values$GOFs, stim_data$GOF[1])
        values$merged_data <- rbind(values$merged_data, values$data)
        
        # Outputting the top of the data table with the newly calculated 
        output$output <- renderTable({
          head(values$data)
        })
        
        output$gofoutput <- renderTable({
          collated_data <- cbind(values$stim_number_list, values$pre_stim_speed_list, values$max_rel_speed_list, values$max_speed_list, values$GOFs)
          colnames(collated_data) <- c("Stim number", "Mean pre-stim speed", "Max relative amplitude", "Max fitted amplitude", "GOF")
          collated_data
        })
        showNotification("Parameters should now be optimised.", type = "message")
        cat("Parameters should now be optimised.\n")
        
        # clear stimdata 
        
        stim_data <- c()
        j <- j+1
      }
      #i <- i+1
    }
  })
  
  # Code to sequentially select individual flies and calculate curves for them 
  observeEvent(input$sequentialcalc, {
    
    for(i in as.numeric(input$stimulations)) {
      values$stim_number <- i
      showNotification("Retrieving data for stimulation number:", i, type = "message")
      cat("Retrieving data for stimulation number:", i, "\n")
      
      data <- sheet()
      data$time <- as.numeric(data$time)
      stim_data <- get_data(data, i, as.numeric(input$pretime), as.numeric(input$prestim), as.numeric(input$stimgap), as.numeric(input$poststim))
      req(stim_data)  # Ensure data is not NULL before continuing
      
      showNotification("Data retrieved successfully.", type = "message")
      cat("Data retrieved successfully.\n")
      
      
      # Selecting individual flies
      for(fly in 2:ncol(stim_data)) { 
        if (ncol(stim_data) < 10) {
          stop("Not enough columns in data")
        }
        
        else{
          fly_data <- stim_data[, c(1, fly)]
          #print(stim_data)
          
          
          
          
          
          # Data needs to be tidied into long format
          fly_data <- fly_data %>% 
            pivot_longer(!time, names_to = "Group", values_to = "speed")
          # 60 second time window prior to stimulation taken to calculate the pre-stimulation speed
          pre_stimuli_avg_speed <- mean(fly_data$speed[1:as.numeric(input$prestim)])
          # Mean of pre-stimulation speed subtracted from raw speed to give relative speed. Necessary for
          # curve fitting.
          fly_data$relative_speed <- fly_data$speed - pre_stimuli_avg_speed
          
          showNotification(paste("Optimising parameters for", unique(fly_data$Group), sep = " "), type = "message")
          cat("Optimising parameters...\n")
          # Optimizing parameters of the fitted exponential curve.
          # Initial parameter estimates
          initial_params <- c(dt = input$prestim, A0 = 0.1, A1 = 2, tauA = 0.1, tauB = 3)
          
          optim_results <- optim(
            par = initial_params,
            fn = calculate_R_squared,
            data = fly_data,
            method = "L-BFGS-B",
            lower = c(dt = fly_data$time[isolate(input$prestim)-isolate(input$dtflexibility)], A0 = 0, A1 = 0, tauA = -10, tauB = 0),
            upper = c(dt = fly_data$time[isolate(input$prestim)+isolate(input$dtflexibility)], A0 = 10, A1 = 20, tauA = 20, tauB = 50),
            control = list(maxit = 9999999)
          )
          
          optimized_params <- optim_results$par
          names(optimized_params) <- c("dt", "A0", "A1", "tauA", "tauB")
          
          # Making a new column in data which contains the predicted speeds - i.e, the speeds calculated
          # through the exponential model.
          fly_data$predicted_speeds <- sapply(stim_data$time, function(t) model_function(t, optimized_params["dt"], optimized_params["A0"], optimized_params["A1"], optimized_params["tauA"], optimized_params["tauB"]))
          fly_data$GOF <- -optim_results$value
          fly_data$StimNo <- values$stim_number
          # Find the index of the maximum speeds
          max_index <- which.max(fly_data$predicted_speeds)
          max_index_rel <- which.max(fly_data$relative_speed)
          
          # Find the time and value of the maximum speeds
          max_time <- fly_data$time[max_index]
          max_speed <- fly_data$predicted_speeds[max_index]
          max_time_rel <- fly_data$time[max_index_rel]
          max_speed_rel <- fly_data$relative_speed[max_index_rel]
          
          # Updating all of the reactive values
          values$data <- fly_data
          values$pre_stim_speed <- pre_stimuli_avg_speed
          values$pre_stim_speed_list <- c(values$pre_stim_speed_list, pre_stimuli_avg_speed)
          values$optimized_params <- optim_results$par
          values$optim_results <- optim_results
          values$max_time <- max_time
          values$max_speed <- max_speed
          values$max_speed_list <- c(values$max_speed_list, max_speed)
          values$max_rel_speed_list <- c(values$max_rel_speed_list, max_speed_rel)
          values$stim_number_list <- c(values$stim_number_list, values$stim_number)
          values$max_time_rel <- max_time_rel
          values$max_speed_rel <- max_speed_rel
          values$GOFs <- c(values$GOFs, fly_data$GOF[1])
          values$merged_data <- rbind(values$merged_data, values$data)
          
          # Outputting the top of the data table with the newly calculated 
          output$output <- renderTable({
            head(values$data)
          })
          
          output$gofoutput <- renderTable({
            collated_data <- cbind(values$stim_number_list, values$pre_stim_speed_list, values$max_rel_speed_list, values$max_speed_list, values$GOFs)
            colnames(collated_data) <- c("Stim number", "Mean pre-stim speed", "Max relative amplitude", "Max fitted amplitude", "GOF")
            collated_data
          })
          showNotification("Parameters should now be optimised.", type = "message")
          cat("Parameters should now be optimised.\n")
          
          # clear stimdata 
          
          fly_data <- c()
        }
      }
    }
  })
  
  
  output$datatable <- renderTable({
    req(values$merged_data)
    values$merged_data
  })
  
  output$download_data <- downloadHandler(
    filename = function() {
      req(values$merged_data)
      paste0("summary_mean_moving_speed_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$merged_data)
      readr::write_csv(values$merged_data, file)
    }
  )
  
  
  
  # Plot full data graph when button is pressed:
  observeEvent(input$plotfull, {
    
    # Ensure the data is properly imported and the correct group is selected
    data <- sheet()
    data$time <- as.numeric(data$time)
    data <- data %>% pivot_longer(!time, names_to = "Group", values_to = "speed") %>%
      filter(Group %in% input$group)
    
    # Stimulus locations (seconds from t=0)
    stimulus_locations <- stimulus_locations()
    
    # Plot the graph
    output$fulldata <- renderPlot({
      
      # Step 1: Summarize the data dynamically (mean and standard error per time point)
      summarized_data <- summarized_data()
      
      p <- ggplot(data = summarized_data, aes(x = time, y = mean_raw_speed)) +
        geom_line(linewidth = isolate(input$rawlinewidth), colour = "grey80") +
        ylim(c(isolate(input$ylowlimit), isolate(input$ylimit)))+
        labs(title = "Mean Raw Movement Speed", x = "Time (seconds)", y = "Speed (mm/s)") +
        geom_vline(xintercept = stimulus_locations, linetype = "dashed", color = "black", lwd = 0.5) +
        #theme_minimal()+
        theme(
          #legend.position = "none",
          aspect.ratio = 0.625,
          panel.background = element_blank(),
          axis.text.x=element_markdown(angle=60,hjust=0.5, size = isolate(input$axistextfontsize), vjust = 0.5, face = "bold"),
          axis.text.y=element_text(size = isolate(input$axistextfontsize), face = "bold"),
          axis.ticks = element_line(linewidth = 1),
          axis.title.x = element_markdown(face="bold", size = isolate(input$axistitlefontsize)),
          axis.title.y = element_markdown(face="bold", size = isolate(input$axistitlefontsize)),
          plot.title = element_text(size=isolate(input$plottitlefontsize), face="bold.italic", hjust = 0.5),
          axis.line = element_line(colour = "black", linewidth = isolate(input$axislinewidth))
        )
      
      # Optional: overlay a linear regression across per-stimulation peaks (like the Quarto workflow)
      lm_obj <- lm_peak_fit()
      if (!is.null(lm_obj)) {
        if (isTRUE(input$lm_show_ci)) {
          p <- p + geom_ribbon(
            data = lm_obj$pred,
            aes(x = stim_time, ymin = lo, ymax = hi),
            inherit.aes = FALSE,
            alpha = 0.2
          )
        }
        p <- p +
          geom_line(
            data = lm_obj$pred,
            aes(x = stim_time, y = pred),
            inherit.aes = FALSE,
            linetype = "solid",
            linewidth = 1
          ) +
          geom_point(
            data = lm_obj$peaks,
            aes(x = stim_time, y = peak),
            inherit.aes = FALSE,
            size = 2
          )
      }
      
      p
      
    })
  })
  
  # Plot full data graph when button is pressed:
  observeEvent(input$plotfullfitted, {
    
    # Ensure the data is properly imported and the correct group is selected
    data <- sheet()
    data$time <- as.numeric(data$time)
    data <- data %>% pivot_longer(!time, names_to = "Group", values_to = "speed") %>%
      filter(Group %in% input$group)
    
    # Stimulus locations (seconds from t=0)
    stimulus_locations <- stimulus_locations()
    
    # Plot the graph
    output$fulldatafitted <- renderPlot({
      
      # Step 1: Summarize the data dynamically (mean and standard error per time point)
      summarized_data <- summarized_data()
      
      p <- ggplot(data = isolate(values$merged_data), aes(x = time, color = as.factor(StimNo))) +
        #theme_minimal()+
        geom_line(data = summarized_data, aes(y = mean_rel_speed), color = "black", linewidth = isolate(input$plotlinewidth)) +
        #geom_line(aes(y = predicted_speeds), linewidth = isolate(input$fitlinewidth)) +
        # Smoothed line with shaded error band for predicted speeds
        # geom_smooth(aes(y = predicted_speeds), method = "loess", fill = "grey", 
        #             alpha = 0.3, linewidth = isolate(input$fitlinewidth)) +
        
        # Add ribbon for error band and line for mean predicted speeds
        geom_ribbon(data = summarized_data, 
                    aes(ymin = mean_speed - se_speed, ymax = mean_speed + se_speed, 
                        fill = as.factor(StimNo)), alpha = 0.3, show.legend = FALSE) +
        geom_line(data = summarized_data, aes(y = mean_speed), linewidth = isolate(input$fitlinewidth)) +
        geom_line(data = summarized_data, aes(y = mean_raw_speed), linewidth = isolate(input$rawlinewidth),
                  colour = "grey80") +
        
        ylim(c(isolate(input$ylowlimit), isolate(input$ylimit)))+
        labs(title = "Mean Relative Movement Speed with Fit", x = "Time (seconds)", y = "Speed (mm/s)") +
        geom_vline(xintercept = stimulus_locations, linetype = "dashed", color = "black", lwd = 0.5) +
        labs(color = "Stimulation number")+
        #theme_minimal()+
        theme(
          legend.title = element_markdown(size = isolate(input$legendtitlefontsize), face = "bold"),
          legend.text = element_markdown(face = "bold", size = isolate(input$legendtextfontsize)),
          legend.key.size = unit(isolate(input$keysize), "cm"),
          legend.key.width = unit(isolate(input$keywidth),"cm"),
          aspect.ratio = 0.625,
          panel.background = element_blank(),
          axis.text.x=element_markdown(angle=60,hjust=0.5, size = isolate(input$axistextfontsize), vjust = 0.5, face = "bold"),
          axis.text.y=element_text(size = isolate(input$axistextfontsize), face = "bold"),
          axis.ticks = element_line(linewidth = 1),
          axis.title.x = element_markdown(face="bold", size = isolate(input$axistitlefontsize)),
          axis.title.y = element_markdown(face="bold", size = isolate(input$axistitlefontsize)),
          plot.title = element_text(size=isolate(input$plottitlefontsize), face="bold.italic", hjust = 0.5),
          axis.line = element_line(colour = "black", linewidth = isolate(input$axislinewidth))
        )
      
      # Optional: overlay a linear regression across per-stimulation peaks (like the Quarto workflow)
      lm_obj <- lm_peak_fit()
      if (!is.null(lm_obj)) {
        if (isTRUE(input$lm_show_ci)) {
          p <- p + geom_ribbon(
            data = lm_obj$pred,
            aes(x = stim_time, ymin = lo, ymax = hi),
            inherit.aes = FALSE,
            alpha = 0.2
          )
        }
        p <- p +
          geom_line(
            data = lm_obj$pred,
            aes(x = stim_time, y = pred),
            inherit.aes = FALSE,
            linetype = "solid",
            linewidth = 1
          ) +
          geom_point(
            data = lm_obj$peaks,
            aes(x = stim_time, y = peak),
            inherit.aes = FALSE,
            size = 2
          )
      }
      
      p
      
    })
  })
}

shinyApp(ui, server)

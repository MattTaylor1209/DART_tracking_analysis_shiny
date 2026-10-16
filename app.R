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
library(lme4)
library(lmerTest)


options(shiny.maxRequestSize = 200 * 1024^2)  # 200 MB

#### FUNCTIONS ####


# Single-inactivation model from Faville et al. (2015) Sci Rep, eqn 13:
#   V(tau) = (1 - exp(-tau/tauA)) * (A0 + A1 * exp(-tau/tauB)) * H(tau),  tau = t - dt
# Vectorised over t. No clamping of time constants here: positivity and
# tauA <= tauB are enforced by the parameterisation/bounds in fit_dart().
# pmax(t - dt, 0) makes the rise factor exactly 0 before dt, so H(tau) is implicit.
model_function <- function(t, dt, A0, A1, tauA, tauB) {
  tau <- pmax(t - dt, 0)
  (1 - exp(-tau / tauA)) * (A0 + A1 * exp(-tau / tauB))
}

# Objective: sum of squared errors. (Minimising SSE is equivalent to maximising
# R^2 for a given trace, since SST is fixed.)
# p = c(dt, A0, A1, tauA, dtau); tauB = tauA + dtau, so tauB >= tauA always.
sse_fn <- function(p, data) {
  pred <- model_function(data$time, p[["dt"]], p[["A0"]], p[["A1"]],
                         p[["tauA"]], p[["tauA"]] + p[["dtau"]])
  sse <- sum((data$relative_speed - pred)^2)
  if (!is.finite(sse)) return(1e12)  # keep L-BFGS-B away from non-finite regions
  sse
}

# Fit one trace. Returns a list of parameters, fit statistics and bound flags.
#   stim_time  : stimulus onset (s, same clock as data$time)
#   flex       : allowed +/- window for dt around stim_time (s)
#   tauA_min   : smallest tauA allowed (s); ~ sampling interval is sensible
fit_dart <- function(data, stim_time, flex, tauA_min = 0.1) {
  data <- data[is.finite(data$time) & is.finite(data$relative_speed), ]
  
  lower <- c(dt = stim_time - flex, A0 = 0,  A1 = 0,  tauA = tauA_min, dtau = 0)
  upper <- c(dt = stim_time + flex, A0 = 10, A1 = 20, tauA = 20,       dtau = 50)
  start <- c(dt = stim_time,        A0 = 0.1, A1 = 2, tauA = max(1, tauA_min), dtau = 3)
  
  fit <- optim(
    par = start, fn = sse_fn, data = data, method = "L-BFGS-B",
    lower = lower, upper = upper,
    control = list(maxit = 10000)
  )
  
  p   <- fit$par
  sst <- sum((data$relative_speed - mean(data$relative_speed))^2)
  
  # Flag parameters sitting on a bound (tolerance relative to the bound range)
  tol      <- 1e-4 * (upper - lower)
  at_lower <- abs(p - lower) <= tol
  at_upper <- abs(p - upper) <= tol
  
  list(
    dt        = p[["dt"]],
    latency   = p[["dt"]] - stim_time,
    A0        = p[["A0"]],
    A1        = p[["A1"]],
    tauA      = p[["tauA"]],
    tauB      = p[["tauA"]] + p[["dtau"]],
    R2        = 1 - fit$value / sst,
    converged = fit$convergence == 0,
    # tauA at its floor = activation faster than the sampling can resolve
    tauA_at_lower = unname(at_lower["tauA"]),
    # dtau at 0 = tauB collapsed onto tauA (inactivation not separable)
    tauB_eq_tauA  = unname(at_lower["dtau"]),
    tauB_at_upper = unname(at_upper["tauA"] | at_upper["dtau"]),
    dt_at_bound   = unname(at_lower["dt"] | at_upper["dt"]),
    any_at_bound  = any(at_lower | at_upper)
  )
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
      numericInput("dtflexibility", "dt flexibility (± s around stimulus):", min = 0, max = 10, step = 0.5, value = 2),
      numericInput("tauA_min", "Minimum tauA (s; ~ sampling interval):", min = 0.01, max = 5, step = 0.05, value = 0.1),
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
      checkboxInput("lm_random_slope", "Allow per-fly random slope (stimulus)", value = FALSE),
      checkboxInput("lm_compare_models", "Compare random-intercept vs random-slope models", value = TRUE),
      
      
      tags$hr(),
      downloadButton("download_data", "Download full data"),
      downloadButton("download_fits", "Download fit table")
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
            column(8, tableOutput("output"))
          ),
          fluidRow(
            column(
              2,
              actionButton(
                "cleartable",
                HTML("Clear tables<br>(do this between groups<br>or when changing stims)")
              )
            )),
          
          fluidRow(
            tableOutput("gofoutput")
          ),
          
          
          # fluidRow(
          #   column(12, actionButton("calcparam", "Calculate/update parameters"))
          # ),
          
          fluidRow(
            h2("Per-fly curve calculation"),
            column(3, actionButton("sequentialcalc", "Calculate per-fly curves"))
          ),
          
          # fluidRow(
          #   h2("Random fly selection"),
          #   column(3, numericInput("randomno", "How many flies to randomly sample?", value = 10, min = 1, step = 1, max = 20)),
          #   column(3, numericInput("times", "How many times to repeat calculations?", value = 10, min = 1, step = 1, max = 20)),
          #   column(3, actionButton("randomcalc", "Select x flies randomly and calculate parameters"))
          # )
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
          "Individual fly",
          fluidRow(
            column(4, selectInput("fly_select", "Fly", choices = NULL)),
            column(2, br(), actionButton("fly_prev", "◀ Previous")),
            column(2, br(), actionButton("fly_next", "Next ▶")),
            column(4, selectInput("fly_stims", "Stimulation(s)", choices = NULL, multiple = TRUE))
          ),
          fluidRow(
            column(4, checkboxInput("fly_show_raw", "Show raw speed", value = TRUE)),
            column(4, checkboxInput("fly_show_dt", "Mark fitted onset (dt)", value = TRUE)),
            column(4, checkboxInput("fly_free_y", "Auto-scale y-axis", value = FALSE))
          ),
          fluidRow(plotOutput("fly_plot", width = "1200px", height = "600px")),
          fluidRow(column(12, h4("Fit parameters for this fly"), tableOutput("fly_params")))
        ),
        
        tabPanel(
          "Regression (peaks)",
          fluidRow(
            column(12, verbatimTextOutput("lm_summary"))
          ),
          fluidRow(
            column(12, verbatimTextOutput("lm_compare"))
          ),
          fluidRow(
            column(12, tableOutput("peak_table"))
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
  available_sheets <- reactiveVal(character(0))
  
  observeEvent(input$file, {
    req(input$file$datapath)
    
    shs <- readxl::excel_sheets(input$file$datapath)
    available_sheets(shs)
    
    # Keep previous selection if still valid, otherwise fall back to first sheet
    current <- isolate(input$sheet_name)  # change to your input id
    new_sel <- if (!is.null(current) && current %in% shs) current else shs[[1]]
    
    updateSelectInput(
      session,
      inputId = "data",             # change to your input id
      choices = shs,
      selected = new_sel
    )
  }, ignoreInit = TRUE)
  
  
  sheet <- reactive({
    req(input$file$datapath)
    
    shs <- available_sheets()
    req(length(shs) > 0)
    
    # Validate selection
    if (is.null(input$data) || !(input$data %in% shs)) {
      validate(need(FALSE, "Selected sheet not found in the uploaded file. Please choose a valid sheet."))
    }
    
    readxl::read_excel(input$file$datapath, sheet = input$data)
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
    GOFs = NULL,
    R2 = NULL,
    fit_table = NULL   # one row per fly x stim: parameters, GOF and bound flags
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
  
  peak_data <- reactive({
    req(values$merged_data)
    values$merged_data %>%
      group_by(time, StimNo) %>%
      mutate(
        fly_id = row_number()
      ) %>%
      ungroup() %>%
      group_by(time, StimNo, fly_id) %>%
      summarize(
        #fly_id = row_number(),
        mean_speed     = mean(predicted_speeds, na.rm = TRUE),
        mean_raw_speed = mean(speed, na.rm = TRUE),
        se_speed       = sd(predicted_speeds, na.rm = TRUE) / sqrt(n()),
        mean_rel_speed = mean(relative_speed, na.rm = TRUE),
        .groups = "drop"
      )
  })
  
  peak_table <- reactive({
    req(values$merged_data, input$lm_y, input$stimgap, input$pretime)
    
    validate(need("Group" %in% names(values$merged_data),
                  "Couldn't find a fly identifier column. Expected a 'Group' column (from pivot_longer names_to='Group')."))
    
    # Map UI choices to a column that actually exists in values$merged_data
    y_col <- switch(
      input$lm_y,
      "mean_speed"     = "predicted_speeds",
      "mean_rel_speed" = "relative_speed",
      "mean_raw_speed" = "speed",
      input$lm_y  # fallback: allow direct column names if you pass them
    )
    
    validate(need(y_col %in% names(values$merged_data),
                  paste0("Regression variable '", input$lm_y,
                         "' maps to '", y_col,
                         "', but that column wasn't found in merged_data.")))
    
    values$merged_data %>%
      mutate(
        fly_id   = .data$Group,  # fly IDs are the original wide column names (except 'time')
        stim_no  = as.numeric(StimNo),
        stim_no_c = stim_no - 1,
        stim_time = (stim_no - 1) * as.numeric(input$stimgap) + (as.numeric(input$pretime) * 60),
        y = .data[[y_col]]
      ) %>%
      group_by(fly_id, StimNo) %>%
      slice_max(y, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(
        fly_id,
        StimNo,
        stim_no,
        stim_no_c,
        stim_time,
        time_peak = time,
        peak = y
      ) %>%
      arrange(stim_no, fly_id)
  })
  
  peak_summary <- reactive({
    peaks <- peak_table()
    peaks %>%
      group_by(stim_no, stim_no_c, stim_time) %>%
      summarise(
        mean_peak = mean(peak, na.rm = TRUE),
        se_peak   = sd(peak, na.rm = TRUE) / sqrt(sum(!is.na(peak))),
        n         = sum(!is.na(peak)),
        .groups = "drop"
      ) %>%
      arrange(stim_no)
  })
  
  lmer_fits <- reactive({
    if (!isTRUE(input$lm_enable)) return(NULL)
    peaks <- peak_table()
    validate(need(nrow(peaks) >= 2, "Not enough peaks to fit model."))
    
    # Random-intercept model (REML for estimates)
    fit_ri <- lmerTest::lmer(peak ~ stim_no_c + (1 | fly_id), data = peaks, REML = TRUE)
    
    # Random-slope + intercept model (may fail / be singular for some datasets)
    fit_rs <- tryCatch(
      lmerTest::lmer(peak ~ stim_no_c + (stim_no_c | fly_id), data = peaks, REML = TRUE),
      error = function(e) NULL
    )
    
    list(peaks = peaks, fit_ri = fit_ri, fit_rs = fit_rs)
  })
  
  lmer_fit_selected <- reactive({
    x <- lmer_fits()
    if (is.null(x)) return(NULL)
    
    if (isTRUE(input$lm_random_slope) && !is.null(x$fit_rs)) {
      x$fit_rs
    } else {
      x$fit_ri
    }
  })
  
  lmer_compare <- reactive({
    if (!isTRUE(input$lm_enable) || !isTRUE(input$lm_compare_models)) return(NULL)
    x <- lmer_fits()
    if (is.null(x) || is.null(x$fit_rs)) return(NULL)
    
    # Use ML fits for likelihood ratio test
    fit_ri_ml <- lme4::refitML(x$fit_ri)
    fit_rs_ml <- lme4::refitML(x$fit_rs)
    anova(fit_ri_ml, fit_rs_ml)
    
  })
  
  lmer_line_ci <- reactive({
    fit <- lmer_fit_selected()
    if (is.null(fit)) return(NULL)
    
    grid <- peak_summary() %>% select(stim_no, stim_no_c, stim_time)
    
    X  <- model.matrix(~ stim_no_c, data = grid)
    b  <- lme4::fixef(fit)
    V  <- as.matrix(stats::vcov(fit))     # fixed-effect var-cov
    se <- sqrt(diag(X %*% V %*% t(X)))
    
    grid$pred <- as.numeric(X %*% b)
    mult <- as.numeric(input$lm_ci_mult)
    grid$lo <- grid$pred - mult * se
    grid$hi <- grid$pred + mult * se
    
    grid
  })
  
  
  output$peak_table <- renderTable({
    peaks <- peak_table()
    validate(need(!is.null(peaks) && nrow(peaks) > 0, "Enable the regression overlay and calculate parameters to generate peak data."))
    peaks
  })
  
  output$lm_summary <- renderPrint({
    fit <- lmer_fit_selected()
    validate(need(!is.null(fit), "Enable the regression overlay and calculate parameters to fit the model."))
    
    cat("Model used:", if (isTRUE(input$lm_random_slope)) "Random slope + intercept (if available)" else "Random intercept only", "\n")
    if (inherits(fit, "merMod") && lme4::isSingular(fit, tol = 1e-4)) {
      cat("NOTE: Model is singular (one or more random-effect variances ~ 0, or correlations at boundary).\n\n")
    }
    print(summary(fit))
  })
  
  output$lm_compare <- renderPrint({
    cmp <- lmer_compare()
    validate(need(!is.null(cmp), "Enable model comparison (and ensure random-slope model fits) to see results."))
    print(cmp)
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
    values$fit_table <- NULL
    values$R2 <- NULL
    output$gofoutput <- renderTable({
      collated_data
    })
  })
  
  observe({
    # R^2 is now stored directly (the optimiser minimises SSE, not -R^2)
    if (!is.null(values$R2) && is.finite(values$R2) && values$R2 < 0.2) {
      showNotification("GOF is less than 0.2", type = "message")
      cat("GOF is less than 0.2\n")
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
          
          # Stimulus onset for this stim, on the same clock as data$time
          stim_time <- (i - 1) * as.numeric(input$stimgap) + as.numeric(input$pretime) * 60
          
          # Pre-stimulation baseline: everything in the window before stimulus onset
          # (selected by time, not row index, so it doesn't assume 1 Hz sampling)
          pre_stimuli_avg_speed <- mean(fly_data$speed[fly_data$time < stim_time], na.rm = TRUE)
          # Mean of pre-stimulation speed subtracted from raw speed to give relative speed. Necessary for
          # curve fitting.
          fly_data$relative_speed <- fly_data$speed - pre_stimuli_avg_speed
          
          showNotification(paste("Optimising parameters for", unique(fly_data$Group), sep = " "), type = "message")
          cat("Optimising parameters...\n")
          
          fit <- fit_dart(
            fly_data,
            stim_time = stim_time,
            flex      = as.numeric(isolate(input$dtflexibility)),
            tauA_min  = as.numeric(isolate(input$tauA_min))
          )
          
          # Predicted speeds from the fitted model (vectorised, on this fly's own time column)
          fly_data$predicted_speeds <- model_function(fly_data$time, fit$dt, fit$A0, fit$A1, fit$tauA, fit$tauB)
          fly_data$GOF <- fit$R2
          fly_data$StimNo <- values$stim_number
          fly_data$dt <- fit$dt
          fly_data$latency <- fit$latency
          fly_data$A0 <- fit$A0
          fly_data$A1 <- fit$A1
          fly_data$tauA <- fit$tauA
          fly_data$tauB <- fit$tauB
          fly_data$tauA_at_lower <- fit$tauA_at_lower
          fly_data$tauB_eq_tauA <- fit$tauB_eq_tauA
          fly_data$any_at_bound <- fit$any_at_bound
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
          values$optimized_params <- unlist(fit[c("dt", "A0", "A1", "tauA", "tauB")])
          values$R2 <- fit$R2
          values$max_time <- max_time
          values$max_speed <- max_speed
          values$max_speed_list <- c(values$max_speed_list, max_speed)
          values$max_rel_speed_list <- c(values$max_rel_speed_list, max_speed_rel)
          values$stim_number_list <- c(values$stim_number_list, values$stim_number)
          values$max_time_rel <- max_time_rel
          values$max_speed_rel <- max_speed_rel
          values$GOFs <- c(values$GOFs, fly_data$GOF[1])
          values$merged_data <- rbind(values$merged_data, values$data)
          
          # One row per fly x stim. tauA is reported as NA when it sits on its
          # lower bound: the activation is faster than the sampling resolves,
          # so the number is the bound, not an estimate.
          values$fit_table <- rbind(
            values$fit_table,
            data.frame(
              Fly                      = unique(fly_data$Group),
              `Stim number`            = values$stim_number,
              `Mean pre-stim speed`    = pre_stimuli_avg_speed,
              `Max relative amplitude` = max_speed_rel,
              `Max fitted amplitude`   = max_speed,
              GOF                      = fit$R2,
              `Latency (s)`            = fit$latency,
              A0                       = fit$A0,
              A1                       = fit$A1,
              tauA                     = if (fit$tauA_at_lower) NA_real_ else fit$tauA,
              tauB                     = fit$tauB,
              `tauA at floor`          = fit$tauA_at_lower,
              `tauB = tauA`            = fit$tauB_eq_tauA,
              `Any param at bound`     = fit$any_at_bound,
              Converged                = fit$converged,
              check.names = FALSE
            )
          )
          
          # Outputting the top of the data table with the newly calculated
          output$output <- renderTable({
            head(values$data)
          })
          
          output$gofoutput <- renderTable({
            values$fit_table
          }, na = "< floor")
          
          # clear stimdata 
          
          fly_data <- c()
        }
      }
      showNotification("Parameters should now be optimised.", type = "message")
      cat("Parameters should now be optimised.\n")
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
  
  # ---------------------------
  # Individual fly inspection
  # ---------------------------
  
  fly_ids <- reactive({
    req(values$merged_data)
    unique(values$merged_data$Group)
  })
  
  # Keep the fly and stim choices in sync with whatever has been fitted so far
  observe({
    ids <- fly_ids()
    current <- isolate(input$fly_select)
    updateSelectInput(session, "fly_select", choices = ids,
                      selected = if (!is.null(current) && current %in% ids) current else ids[1])
    
    stims <- sort(unique(values$merged_data$StimNo))
    current_s <- isolate(input$fly_stims)
    keep <- intersect(current_s, as.character(stims))
    updateSelectInput(session, "fly_stims", choices = stims,
                      selected = if (length(keep)) keep else stims)
  })
  
  step_fly <- function(direction) {
    ids <- fly_ids()
    pos <- match(input$fly_select, ids)
    if (is.na(pos)) pos <- 1
    new <- ((pos - 1 + direction) %% length(ids)) + 1  # wraps around
    updateSelectInput(session, "fly_select", selected = ids[new])
  }
  observeEvent(input$fly_prev, step_fly(-1))
  observeEvent(input$fly_next, step_fly(+1))
  
  fly_subset <- reactive({
    req(values$merged_data, input$fly_select, input$fly_stims)
    values$merged_data %>%
      filter(Group == input$fly_select, StimNo %in% as.numeric(input$fly_stims))
  })
  
  output$fly_plot <- renderPlot({
    d <- fly_subset()
    validate(need(nrow(d) > 0, "No fitted data for this fly / stimulation."))
    
    stim_times <- stimulus_locations()[sort(unique(d$StimNo))]
    onsets <- d %>% distinct(StimNo, dt)
    
    p <- ggplot(d, aes(x = time, colour = as.factor(StimNo), group = StimNo))
    
    if (isTRUE(input$fly_show_raw)) {
      p <- p + geom_line(aes(y = speed), colour = "grey80", linewidth = input$rawlinewidth)
    }
    
    p <- p +
      geom_hline(yintercept = 0, colour = "grey60", linewidth = 0.3) +
      geom_line(aes(y = relative_speed), colour = "black", linewidth = input$plotlinewidth) +
      geom_line(aes(y = predicted_speeds), linewidth = input$fitlinewidth) +
      geom_vline(xintercept = stim_times, linetype = "dashed", colour = "black", lwd = 0.5)
    
    if (isTRUE(input$fly_show_dt)) {
      p <- p + geom_vline(data = onsets, aes(xintercept = dt, colour = as.factor(StimNo)),
                          linetype = "dotted", lwd = 0.7, show.legend = FALSE)
    }
    
    if (!isTRUE(input$fly_free_y)) {
      p <- p + coord_cartesian(ylim = c(input$ylowlimit, input$ylimit))
    }
    
    p +
      labs(title = input$fly_select, x = "Time (seconds)", y = "Speed (mm/s)",
           colour = "Stimulation number") +
      theme(
        legend.title = element_markdown(size = input$legendtitlefontsize, face = "bold"),
        legend.text = element_markdown(face = "bold", size = input$legendtextfontsize),
        legend.key.size = unit(input$keysize, "cm"),
        legend.key.width = unit(input$keywidth, "cm"),
        aspect.ratio = 0.625,
        panel.background = element_blank(),
        axis.text.x = element_text(size = input$axistextfontsize, face = "bold"),
        axis.text.y = element_text(size = input$axistextfontsize, face = "bold"),
        axis.ticks = element_line(linewidth = 1),
        axis.title = element_text(face = "bold", size = input$axistitlefontsize),
        # always show the fly ID here, even if the pooled plot title size is 0
        plot.title = element_text(size = max(input$plottitlefontsize, 16), face = "bold", hjust = 0.5),
        axis.line = element_line(colour = "black", linewidth = input$axislinewidth)
      )
  })
  
  output$fly_params <- renderTable({
    req(values$fit_table, input$fly_select)
    values$fit_table %>%
      filter(Fly == input$fly_select, `Stim number` %in% as.numeric(input$fly_stims))
  }, na = "< floor", digits = 3)
  
  output$download_fits <- downloadHandler(
    filename = function() {
      paste0("per_fly_fit_parameters_", Sys.Date(), ".csv")
    },
    content = function(file) {
      req(values$fit_table)
      readr::write_csv(values$fit_table, file)
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
      # Optional: overlay lmer line fitted on per-fly peaks,
      # but PLOT only mean±SE peaks
      if (isTRUE(input$lm_enable)) {
        ln <- lmer_line_ci()
        ps <- peak_summary()
        
        if (!is.null(ln)) {
          
          # mean peak ± SE points (what you want to see)
          p <- p +
            geom_errorbar(
              data = ps,
              aes(x = stim_time, ymin = mean_peak - se_peak, ymax = mean_peak + se_peak),
              inherit.aes = FALSE,
              width = 0
            ) +
            geom_point(
              data = ps,
              aes(x = stim_time, y = mean_peak),
              inherit.aes = FALSE,
              size = 3
            )
          
          # line + CI for fixed effect
          if (isTRUE(input$lm_show_ci)) {
            p <- p + geom_ribbon(
              data = ln,
              aes(x = stim_time, ymin = lo, ymax = hi),
              inherit.aes = FALSE,
              alpha = 0.2
            )
          }
          
          p <- p + geom_line(
            data = ln,
            aes(x = stim_time, y = pred),
            inherit.aes = FALSE,
            linewidth = 1
          )
        }
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
      # Optional: overlay lmer line fitted on per-fly peaks,
      # but PLOT only mean±SE peaks
      if (isTRUE(input$lm_enable)) {
        ln <- lmer_line_ci()
        ps <- peak_summary()
        
        if (!is.null(ln)) {
          
          # mean peak ± SE points (what you want to see)
          p <- p +
            geom_errorbar(
              data = ps,
              aes(x = stim_time, ymin = mean_peak - se_peak, ymax = mean_peak + se_peak),
              inherit.aes = FALSE,
              width = 0
            ) +
            geom_point(
              data = ps,
              aes(x = stim_time, y = mean_peak),
              inherit.aes = FALSE,
              size = 3
            )
          
          # line + CI for fixed effect
          if (isTRUE(input$lm_show_ci)) {
            p <- p + geom_ribbon(
              data = ln,
              aes(x = stim_time, ymin = lo, ymax = hi),
              inherit.aes = FALSE,
              alpha = 0.2
            )
          }
          
          p <- p + geom_line(
            data = ln,
            aes(x = stim_time, y = pred),
            inherit.aes = FALSE,
            linewidth = 1
          )
        }
      }
      
      p
      
    })
  })
}

shinyApp(ui, server)
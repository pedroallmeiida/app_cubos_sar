# App Shiny: clique em um pixel do cubo e veja histograma, ACF, PACF,
# serie com ajustes e tabela de metricas.
# Execucao: shiny::runApp("app_cubos.R")

library(shiny)
library(terra)
library(ggplot2)
library(forecast)
library(dplyr)
library(purrr)
library(tibble)
library(maxLik)
library(invgamma)
library(pso)

repo_root <- "C:/Users/pedro/Documents/github/artigo_kb_arma/codigo_cubos"; repo_cubos <- repo_root; 
#repo_cubos <- "E:/Snap/dados_area_preservacao_brasilia_spk"


invisible(lapply(list.files(file.path(repo_root, "modulos"), pattern = "\\.R$",
                            full.names = TRUE), source))

cubos <- list(
  VH  = file.path(repo_cubos, "cubos/cube_vh.tif"),
  VV  = file.path(repo_cubos, "cubos/cube_vv.tif"),
  RVI = file.path(repo_cubos, "cubos/cube_rvi.tif")
)

cores <- c("Observado" = "black", "KB ARMA (sqrt)" = "red", "KB ARMA (W)" = "blue",
           "KB ARMA (log)" = "green", "Gama ARMA (sqrt)" = "deepskyblue",
           "Gama ARMA (W)" = "hotpink", "Gama ARMA (log)" = "orange",
           "G0 ARMA (sqrt)" = "purple", "G0 ARMA (W)" = "brown",
           "G0 ARMA (log)" = "darkcyan", "ARIMA" = "gray40", "ETS" = "gold3")

ajustar <- function(z, L, criterio, incluir) {
  f <- list(
    "KB ARMA (sqrt)"   = function() auto.K_sqrt(z, L = L, criterio = criterio)$mu.hat,
    "KB ARMA (W)"      = function() auto.K_W(z, L = L, criterio = criterio)$mu.hat,
    "KB ARMA (log)"    = function() auto.K_log(z, L = L, criterio = criterio)$mu.hat,
    "Gama ARMA (sqrt)" = function() auto.gamma_sqrt(z, L = L, criterio = criterio)$mu.hat,
    "Gama ARMA (W)"    = function() auto.gamma_W(z, L = L, criterio = criterio)$mu.hat,
    "Gama ARMA (log)"  = function() auto.gamma_log(z, L = L, criterio = criterio)$mu.hat,
    "G0 ARMA (sqrt)"   = function() auto.g0_sqrt(z, L = L, criterio = criterio)$mu.hat,
    "G0 ARMA (W)"      = function() auto.g0_W(z, L = L, criterio = criterio)$mu.hat,
    "G0 ARMA (log)"    = function() auto.g0_log(z, L = L, criterio = criterio)$mu.hat,
    "ARIMA"            = function() as.numeric(fitted(auto.arima(z))),
    "ETS"              = function() as.numeric(ets(z)$fitted)
  )[incluir]
  out <- list()
  for (nm in names(f)) {
    withProgress(message = paste("Ajustando", nm), value = 0, {
      out[[nm]] <- tryCatch({
        v <- as.numeric(f[[nm]]())
        if (length(v) < length(z)) c(rep(NA_real_, length(z) - length(v)), v) else v
      }, error = function(e) rep(NA_real_, length(z)))
    })
  }
  out
}

metricas <- function(z, ajustes) {
  imap_dfr(ajustes, function(previsto, modelo) {
    ok <- is.finite(z) & is.finite(previsto)
    if (!any(ok)) return(tibble(Modelo = modelo, RMSE = NA_real_, MAE = NA_real_,
                                MAPE = NA_real_, sMAPE = NA_real_, MASE = NA_real_))
    e <- z[ok] - previsto[ok]
    esc <- mean(abs(diff(z)), na.rm = TRUE)
    den <- abs(z[ok]) + abs(previsto[ok])
    tibble(Modelo = modelo,
           RMSE  = sqrt(mean(e^2)),
           MAE   = mean(abs(e)),
           MAPE  = mean(abs(e / z[ok])) * 100,
           sMAPE = mean(ifelse(den == 0, 0, 2 * abs(e) / den)),
           MASE  = if (is.finite(esc) && esc > 0) mean(abs(e)) / esc else NA_real_)
  })
}

ui <- fluidPage(
  titlePanel("Cubo de imagens: analise de series por pixel"),
  sidebarLayout(
    sidebarPanel(
      width = 3,
      selectInput("cubo", "Cubo", names(cubos), "VH"),
      numericInput("L", "Numero de looks (L)", 1, min = 1),
      selectInput("criterio", "Criterio", c("AIC", "BIC"), "AIC"),
      checkboxGroupInput("modelos", "Modelos", names(cores)[-1],
                         selected = c("KB ARMA (sqrt)", "KB ARMA (W)", "KB ARMA (log)",
                                      "Gama ARMA (sqrt)", "Gama ARMA (W)", "Gama ARMA (log)")),
      actionButton("ajustar", "Ajustar modelos", class = "btn-primary"),
      helpText("Clique na imagem para escolher o pixel; depois clique em 'Ajustar modelos'."),
      textOutput("coord")
    ),
    mainPanel(
      width = 9,
      fluidRow(
        column(3, plotOutput("img", click = "clique", height = 300)),
        column(3, plotOutput("hist", height = 300)),
        column(3, plotOutput("acf", height = 300)),
        column(3, plotOutput("pacf", height = 300))
      ),
      plotOutput("serie", height = 380),
      h4("Metricas"),
      tableOutput("tabela")
    )
  )
)

server <- function(input, output, session) {
  cubo <- reactive(rast(cubos[[input$cubo]]))
  media <- reactive(app(cubo(), mean, na.rm = TRUE))

  pt <- reactiveVal(NULL)
  observeEvent(input$clique, pt(c(x = input$clique$x, y = input$clique$y)))
  observeEvent(input$cubo, pt(NULL))

  serie <- reactive({
    req(pt())
    v <- as.numeric(terra::extract(cubo(), matrix(pt(), ncol = 2))[1, -1])
    validate(need(sum(is.finite(v)) > 10, "Pixel fora da imagem ou sem dados."))
    v
  })

  resultado <- eventReactive(input$ajustar, {
    z <- serie()
    validate(need(length(input$modelos) > 0, "Selecione ao menos um modelo."))
    aj <- ajustar(z, input$L, input$criterio, input$modelos)
    list(aj = aj, tab = metricas(z, aj), z = z)
  })

  output$coord <- renderText({
    if (is.null(pt())) "Nenhum pixel selecionado."
    else sprintf("Pixel: x = %.2f, y = %.2f", pt()[1], pt()[2])
  })

  output$img <- renderPlot({
    plot(media(), main = paste("Media -", input$cubo),
         col = colorRampPalette(c("black", "white"))(100))
    if (!is.null(pt())) points(pt()[1], pt()[2], col = "red", pch = 15, cex = 1.5)
  })

  output$hist <- renderPlot({
    hist(serie(), breaks = 15, main = "Histograma", xlab = input$cubo,
         ylab = "Frequencia", col = "lightblue", border = "black")
  })
  output$acf  <- renderPlot(Acf(serie(), main = "ACF"))
  output$pacf <- renderPlot(Pacf(serie(), main = "PACF"))

  output$serie <- renderPlot({
    r <- resultado()
    df <- data.frame(t = seq_along(r$z), Observado = r$z, r$aj, check.names = FALSE)
    long <- tidyr::pivot_longer(df, -t, names_to = "Modelo", values_to = "y")
    ggplot(long, aes(t, y, colour = Modelo)) +
      geom_line(linewidth = 0.7) +
      scale_colour_manual(values = cores) +
      labs(title = paste("Serie", input$cubo, "e modelos ajustados"),
           x = "Tempo", y = input$cubo, colour = NULL) +
      theme_minimal()
  })

  output$tabela <- renderTable(resultado()$tab, digits = 4)
}

shinyApp(ui, server)

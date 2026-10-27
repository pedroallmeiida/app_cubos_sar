# RVI - G0 

rm(list=ls())

## definindo caminho do repositorio
repo_root <- "C:/Users/pedro/Documents/github/artigo_kb_arma/codigo_cubos"
repo_cubos <- repo_root


## cubos do HD externo
#repo_cubos <- "E:/Snap/dados_area_preservacao_brasilia_spk"


# C:/Users/Elize/Downloads/Projeto - Series Temporais/Funcoes - Processamento/

## diretorio
modules_dir <- file.path(repo_root, "modulos")

## carregando os modulos
invisible(lapply(list.files(modules_dir, pattern = "\\.R$", full.names = TRUE), source))

#### Funcoes para o processo Gama
library(maxLik)
library(invgamma)
library(pso)
library(ggplot2)

# Carregando os cubos de dados
library(terra)



cube_vh <- rast( file.path(repo_cubos,"cubos/cube_vh.tif") )
cube_vv <- rast( file.path(repo_cubos,"cubos/cube_vv.tif"))
rvi <- rast( file.path(repo_cubos,"cubos/cube_rvi.tif"))


##### ANALISE INTENSIDADES 

##### Plotar a imagem - VH
cube_vh_mean <- app(rvi, mean, na.rm = TRUE)
cube_vh_mean

if (identical(names(grDevices::dev.cur()), "RStudioGD")) {
  grDevices::dev.new(noRStudioGD = TRUE)
}

plot(cube_vh_mean,
     main   = "Intensidade canal VH",
     col    =  colorRampPalette(c("black", "white"))(100),
     legend = TRUE)


pt <- tryCatch(click(cube_vh_mean, n = 1, xy = TRUE, show = FALSE),
               error = function(e) NULL)

points(pt$x, pt$y, col = "red", pch = 15, cex = 1.5)
ts_vh  <- terra::extract(cube_vh, pt[, c("x", "y")])
z_vh <- as.numeric(ts_vh[1, -1])

hist(z_vh, breaks = 15, main = "Histograma da Intensidade VH", xlab = "Intensidade VH", ylab = "Frequência", col = "lightblue", border = "black")

## Modelos KB-ARMA
auto_k_sqrt= auto.K_sqrt(z_vh, L=1, criterio = 'AIC')
auto_k_W = auto.K_W(z_vh, L=1, criterio = 'AIC')
auto_k_log = auto.K_log(z_vh, L=1, criterio = 'AIC')


## Modelos G0-ARMA
auto_g0_sqrt= auto.g0_sqrt(z_vh, L=1, criterio = 'AIC')
auto_g0_W = auto.g0_W(z_vh, L=1, criterio = 'AIC')
auto_g0_log = auto.g0_log(z_vh, L=1, criterio = 'AIC')


## Modelos GAMA-ARMA
auto_gama_sqrt= auto.gamma_sqrt(z_vh, L=1, criterio = 'AIC')
auto_gama_W = auto.gamma_W(z_vh, L=1, criterio = 'AIC')
auto_gama_log = auto.gamma_log(z_vh, L=1, criterio = 'AIC')


## Classicos ARIMA e ETS
auto_arima = auto.arima( z_vh )
ajuste_ets = ets(z_vh)


# Grafico da serie - Intensidade VH
z_hat_vh_auto_k_sqrt = auto_k_sqrt$mu.hat
z_hat_vh_auto_k_W = auto_k_W$mu.hat
z_hat_vh_auto_k_log = auto_k_log$mu.hat
z_hat_vh_auto_gama_sqrt = auto_gama_sqrt$mu.hat
z_hat_vh_auto_gama_W = auto_gama_W$mu.hat
z_hat_vh_auto_gama_log = auto_gama_log$mu.hat
z_hat_vh_auto_g0_sqrt = auto_g0_sqrt$mu.hat
z_hat_vh_auto_g0_W = auto_g0_W$mu.hat
z_hat_vh_auto_g0_log = auto_g0_log$mu.hat
z_hat_arima = fitted(auto_arima)
z_hat_ets = ajuste_ets$fitted


df_vh <- data.frame(
  t      = seq_along(z_vh),
  z      = z_vh,
  z_hat_vh_auto_k_sqrt  = z_hat_vh_auto_k_sqrt,
  z_hat_vh_auto_k_W  = z_hat_vh_auto_k_W,
  z_hat_vh_auto_k_log  = z_hat_vh_auto_k_log,
  z_hat_vh_auto_g0_sqrt  = z_hat_vh_auto_g0_sqrt,
  z_hat_vh_auto_g0_W  = z_hat_vh_auto_g0_W,
  z_hat_vh_auto_g0_log  = z_hat_vh_auto_g0_log,
  z_hat_vh_auto_gama_sqrt  = z_hat_vh_auto_gama_sqrt,
  z_hat_vh_auto_gama_W  = z_hat_vh_auto_gama_W,
  z_hat_vh_auto_gama_log  = z_hat_vh_auto_gama_log,
  z_hat_arima  = z_hat_arima,
  z_hat_ets  = z_hat_ets
)

ggplot(df_vh, aes(x = t)) +
  geom_line(aes(y = z,     colour = "Observado"), linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_k_sqrt, colour = "KB ARMA (sqrt)"),  linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_k_W, colour = "KB ARMA (W)"),  linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_k_log, colour = "KB ARMA (log)"),  linewidth = 0.7) +
  #geom_line(aes(y = z_hat_vh_auto_g0_sqrt, colour = "auto karma (sqrt)"),  linewidth = 0.7) +
  #geom_line(aes(y = z_hat_vh_auto_g0_W, colour = "auto karma (W)"),  linewidth = 0.7) +
  #geom_line(aes(y = z_hat_vh_auto_g0_log, colour = "auto karma (log)"),  linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_gama_sqrt, colour = "Gama ARMA (sqrt)"),  linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_gama_W, colour = "Gama ARMA (W)"),  linewidth = 0.7) +
  geom_line(aes(y = z_hat_vh_auto_gama_log, colour = "Gama ARMA (log)"),  linewidth = 0.7) +
  scale_colour_manual(values = c("Observado" = "black", "KB ARMA (sqrt)" = "red", "KB ARMA (W)" = "blue", "KB ARMA (log)" = "green",  "Gama ARMA (sqrt)" = "lightblue", "Gama ARMA (W)" = "pink", "Gama ARMA (log)" = "orange")) +
  labs(title   = "Série Intensidade VH - Ajuste KARMA",
       x       = "Tempo",
       y       = "Intensidade VH",
       colour  = NULL) +
  theme_minimal()



### Acuracia 

library(forecast)
library(dplyr)
library(purrr)
library(tibble)

# Modelos para avaliação
modelos <- c(
  "K-SQRT" = "z_hat_vh_auto_k_sqrt",
  "K-W"    = "z_hat_vh_auto_k_W",
  "K-LOG"  = "z_hat_vh_auto_k_log",
  #"G0-SQRT" = "z_hat_vh_auto_g0_sqrt",
  #"G0-W"    = "z_hat_vh_auto_g0_W",
  #"G0-LOG"  = "z_hat_vh_auto_g0_log",
  "Gama-SQRT" = "z_hat_vh_auto_gama_sqrt",
  "Gama-W"    = "z_hat_vh_auto_gama_W",
  "Gama-LOG"  = "z_hat_vh_auto_gama_log"
  #"ARIMA"  = "z_hat_arima",
  #"ETS"  = "z_hat_ets"
)

# Calcular métricas para todos os modelos
tabela_accuracy <- imap_dfr(modelos, function(coluna, modelo) {
  require(dplyr)
  observado <- df_vh$z
  previsto <- df_vh[[coluna]]
  metricas <- accuracy(
    previsto,
    observado
  )

  validos <- is.finite(observado) & is.finite(previsto)
  erro_absoluto <- abs(observado[validos] - previsto[validos])
  denominador_smape <- abs(observado[validos]) + abs(previsto[validos])
  smape <- if (length(erro_absoluto) > 0) {
    mean(ifelse(denominador_smape == 0, 0, 2 * erro_absoluto / denominador_smape))
  } else {
    NA_real_
  }

  escala_mase <- mean(abs(diff(observado)), na.rm = TRUE)
  mase <- if (is.finite(escala_mase) && escala_mase > 0 && length(erro_absoluto) > 0) {
    mean(erro_absoluto) / escala_mase
  } else {
    NA_real_
  }
  
  as_tibble(metricas) %>%
    mutate(sMAPE = smape, MASE = mase, Modelo = modelo, .before = 1)
})

# Visualizar resultados
tabela_accuracy %>% select( Modelo, MAPE, sMAPE, MASE )








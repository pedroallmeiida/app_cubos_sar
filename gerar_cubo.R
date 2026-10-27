# RVI - G0 

rm(list=ls())

## definindo caminho do repositorio
repo_root <- "C:/Users/pedro/Downloads/arquivos/arquivos/codigos"

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
library(terra)


### Carregamento os dados das imagens 
#folder_filtro <- "C:/Users/pedro/Downloads/arquivos/arquivos/codigos/cena_area_de_preservacao_BRASILIA"

folder_filtro <- "E:/Snap/dados_area_preservacao_brasilia_spk"

files_tif <- list.files(folder_filtro, pattern = "\\.tif$", full.names = TRUE)

if (length(files_tif) == 0) stop("Nenhum arquivo '.tif' encontrado.")

cat("Carregando", length(files_tif), "imagens...\n")

ref <- rast(files_tif[1])

# Carrega banda com resample se os extents divergirem
load_band <- function(f, band, ref) {
  r <- rast(f, lyrs = band)
  if (!compareGeom(r, ref, stopOnError = FALSE))
    r <- resample(r, ref, method = "bilinear")
  return(r)
}

vh_list <- lapply(files_tif, load_band, band = 1, ref = ref[[1]])
vv_list <- lapply(files_tif, load_band, band = 2, ref = ref[[2]])

## Cria o cubo dos dados
cube_vh <- rast(vh_list)
cube_vv <- rast(vv_list)



# RVI = (4 * VH) / (VH + VV)
cat("Calculando cubo de dados: RVI...\n")
cube_rvi <- (4 * cube_vh) / (cube_vh + cube_vv)
#cube_rvi[cube_rvi < 0 | cube_rvi > 4] <- NA


# gerar_cubo.R
writeRaster(cube_vh, "E:/Snap/dados_area_preservacao_brasilia_spk/cubos/cube_vh.tif", overwrite = TRUE)
writeRaster(cube_vv, "E:/Snap/dados_area_preservacao_brasilia_spk/cubos/cube_vv.tif", overwrite = TRUE)
writeRaster(cube_rvi, "E:/Snap/dados_area_preservacao_brasilia_spk/cubos/cube_rvi.tif", overwrite = TRUE)


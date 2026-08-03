# ============================================================
# TP2 Macroeconometría
# ============================================================

library(ggplot2)
library(tseries)
library(kableExtra)
library(vars)
library(tidyverse)
library(readxl)
library(lubridate)
library(patchwork)
library(stargazer)
library(urca)

# Para TCN, IPC y Expectativas se usa la librería rbcb que interactúa directo con la api del BCB.
# https://wilsonfreitas.github.io/rbcb/

# install.packages("rbcb")

library(rbcb)

# La consigna sólo dice que tienen que ser mínimo 15 años, por las dudas tomo del 2000 al 2025

fecha_inicio <- "2000-01-01"
fecha_fin    <- "2025-12-01"


# Helper: lleva una serie diaria/mensual a mensual (promedio simple del mes).
# df debe tener columnas 'date' (Date) y la columna de valor indicada.
# Para series ya mensuales el promedio es inocuo (una obs por mes).
to_monthly <- function(df, value_col, from = fecha_inicio, to = fecha_fin) {
  value_col <- rlang::ensym(value_col)
  df %>%
    filter(!is.na(date), !is.na(!!value_col)) %>%
    mutate(date_m = floor_date(date, "month")) %>%
    group_by(date_m) %>%
    summarise(value = mean(!!value_col, na.rm = TRUE), .groups = "drop") %>%
    arrange(date_m) %>%
    filter(date_m >= as.Date(from), date_m <= as.Date(to))
}


# ------------------------------------------------------------
# 1) TCN — Tipo de cambio nominal R$/US$ (BCB série 3696, prom. mensual)
#    Transformación: 100 * log(TCN)
# ------------------------------------------------------------
tcn <- get_series(c(tcn = 3696), start_date = fecha_inicio)

tcn_m <- tcn %>%
  transmute(date = as.Date(date), value = tcn) %>%
  to_monthly(value) %>%
  transmute(date_m, tcn = 100 * log(value))


# ------------------------------------------------------------
# 2) IPC — IPCA desestacionalizado
# BCB: série 433 (IPCA, variación % mensual)
#    Transformación final: 100 * log(IPC_sa)
# ------------------------------------------------------------

ipca_var <- get_series(c(ipca = 433), start_date = "1999-07-01")   # tomo algunos meses de más para desestacionalizar bien


# la serie está en % mensual. Convierto a número índice para que sea más fácil trabajar.

ipc_idx <- ipca_var %>%
  arrange(date) %>%
  mutate(ipc_nivel = 100 * cumprod(1 + ipca / 100))                # número índice (base 100 = jun-1999)

# desestacionalizo con X13 (paquete seasonal) -> no encontré series ya desestacionalizadas en IBGE o BCB

ipc_ts <- ts(ipc_idx$ipc_nivel, 
             start = c(year(min(ipc_idx$date)), month(min(ipc_idx$date))), 
             frequency = 12)

fit <- seas(ipc_ts)

ipc_idx <- ipc_idx %>% 
  mutate(ipc_nivel_sa = as.numeric(final(fit)))

ipc_m <- ipc_idx %>%
  transmute(date = as.Date(date), value = ipc_nivel_sa) %>%
  to_monthly(value) %>%
  transmute(date_m, ipc = 100 * log(value))

# ------------------------------------------------------------
# 3) Índice de términos de intercambio de commodities
#    IMF Commodity Terms of Trade: serie NETA (export - import),
#    ponderación por comercio total, pesos móviles, mensual:
#    BRA.CEMPI_CTOTNX_TT.R_RW_IX.M
#    Transformación: 100 * log(índice)
# ------------------------------------------------------------

raw_comm <- read_csv("Indice de precios de commodities exportados.csv",
                     locale = locale(encoding = "UTF-8"))

cols_fecha <- names(raw_comm) %>% str_subset("^[0-9]{4}-M[0-9]{2}$")   # sólo columnas mensuales

ctot_m <- raw_comm %>%
  filter(SERIES_CODE == "BRA.CEMPI_CTOTNX_TT.R_RW_IX.M") %>%
  slice(1) %>%
  select(all_of(cols_fecha)) %>%
  pivot_longer(everything(), names_to = "fecha_str", values_to = "ctot") %>%
  mutate(
    date = as.Date(paste0(str_replace(fecha_str, "-M", "-"), "-01")),
    ctot = as.numeric(ctot)
  ) %>%
  transmute(date, value = ctot) %>%
  to_monthly(value) %>%
  transmute(date_m, ctot = 100 * log(value))

# ------------------------------------------------------------
# 4) EMBI Brasil (EMBI+ JP Morgan, IPEAdata, diario -> promedio mensual)
#    Se mantiene en NIVEL (pbs).
# ------------------------------------------------------------
embi_m <- read_csv("ipeadata[07-07-2026-07-09].csv",
                   locale = locale(encoding = "UTF-8")) %>%
  rename(fecha = 1, embi = 2) %>%
  transmute(date = dmy(fecha), value = as.numeric(embi)) %>%
  to_monthly(value) %>%
  transmute(date_m, embi = value)

# ------------------------------------------------------------
# 5) Expectativa de INFLACIÓN a 12 meses (Focus, BCB)
#    Se toma IPCA, serie suavizada, mediana. Diaria -> promedio mensual.
#    Transformación: un registro de x% -> 100 * log(1 + x/100)
# ------------------------------------------------------------
exp_inf <- get_market_expectations("inflation-12-months", start_date = fecha_inicio)

exp_inf_m <- exp_inf %>%
  filter(Indicador == "IPCA", Suavizada == "S", baseCalculo == 0) %>%
  transmute(date = as.Date(Data), value = Mediana) %>%
  to_monthly(value) %>%
  transmute(date_m, exp_inf = 100 * log(1 + value / 100))

# ------------------------------------------------------------
# 6) Expectativa de TIPO DE CAMBIO a 12 meses (Focus, BCB)  [OPCIONAL]
#    Expectativas mensuales de Câmbio: para cada fecha de encuesta se toma
#    la proyección cuyo mes de referencia está 12 meses adelante.
#    Diaria -> promedio mensual.  Transformación: 100 * log(valor)
# ------------------------------------------------------------
exp_cam <- get_monthly_market_expectations("Câmbio", start_date = fecha_inicio)

exp_cam_m <- exp_cam %>%
  filter(base == 0) %>%
  mutate(ref_objetivo = format(floor_date(date, "month") %m+% months(12), "%Y-%m")) %>%
  filter(as.character(reference_date) == ref_objetivo) %>%
  transmute(date, value = median) %>%
  to_monthly(value) %>%
  transmute(date_m, exp_cam = 100 * log(value))

# ------------------------------------------------------------
# 7) Shock de política monetaria de la Fed — Jarocinski & Karadi (MP_median)
#    Archivo mensual con columnas year, month.
# ------------------------------------------------------------
mp_fed_m <- read_csv("shocks_fed_jk_m.csv") %>%
  transmute(date = make_date(year, month, 1), value = MP_median) %>%
  to_monthly(value) %>%
  transmute(date_m, mp_fed = value)

# ------------------------------------------------------------
# 8) Shock de oferta de petróleo — Baumeister & Hamilton (mensual)
# ------------------------------------------------------------
oil_m <- read_excel("BH2_supply_shocks.xlsx", sheet = "oil supply shocks",
                    skip = 2, col_names = FALSE) %>%
  select(date = 1, value = 2) %>%
  mutate(date = as.Date(date)) %>%
  to_monthly(value) %>%
  transmute(date_m, oil_shock = value)

# ------------------------------------------------------------
# 9) Excess Bond Premium (EBP) — mensual, primer día de mes
# ------------------------------------------------------------
ebp_m <- read_csv("ebp_csv.csv") %>%
  transmute(date = as.Date(date), value = ebp) %>%
  to_monthly(value) %>%
  transmute(date_m, ebp = value)

# ============================================================
# BASE UNIFICADA (mensual)
# ============================================================
Base_Unificada <- list(
  tcn_m, ipc_m, ctot_m, embi_m, exp_inf_m, exp_cam_m,   # locales
  mp_fed_m, oil_m, ebp_m                                # globales
) %>%
  reduce(full_join, by = "date_m") %>%
  arrange(date_m) %>%
  filter(date_m >= as.Date(fecha_inicio), date_m <= as.Date(fecha_fin)) %>%
  mutate(mes = format(date_m, "%Y-%m")) %>%
  relocate(mes, date_m)

#saveRDS(Base_Unificada, "Base_Unificada.rds")   # desmutear si se cambió algo en la construcción

# ------------------------------------------------------------
# Cobertura de cada serie (primer y último mes con dato + N)
# Sirve para elegir la muestra común final.
# ------------------------------------------------------------
cobertura <- Base_Unificada %>%
  select(-mes) %>%
  pivot_longer(-date_m, names_to = "serie", values_to = "valor") %>%
  filter(!is.na(valor)) %>%
  group_by(serie) %>%
  summarise(obs   = n(),
            desde = format(min(date_m), "%Y-%m"),
            hasta = format(max(date_m), "%Y-%m"),
            .groups = "drop")
print(cobertura)

# las series de expectations arrancan en 2001-12. La de EMBI se corta en 2024-07


# correr desde acá para no generar la base otra vez!

Base_Unificada <- read_rds("Base_Unificada.rds")





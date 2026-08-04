# ============================================================
# TP2 Macroeconometría — PUNTOS 1 y 2
#
#   Punto 1 : (a) gráficos de dos ejes contra el TCN
#             (b) tests de raíz unitaria -> orden de integración
#   Punto 2 : VAR bivariado (D.TCN, D.IPC), Cholesky con TCN primero.
#             IRF acumuladas + ERPT incondicional.
#
# Nota: 'vars' importa MASS, que enmascara dplyr::select. Por eso se
# carga vars PRIMERO y aun así se usa el prefijo dplyr:: donde importa.
# ============================================================

library(vars)
library(tidyverse)
library(lubridate)
library(patchwork)
library(urca)
library(kableExtra)
library(sandwich) 
library(lmtest)
library(AER)
library(car)
base <- read_rds("Base_Unificada.rds") %>%
  mutate(lembi = 100 * log(embi))   # el nivel en pbs está dominado por 2002


# ------------------------------------------------------------
# Parámetros globales
# ------------------------------------------------------------

H     <- 24     # horizonte máximo de las IRF
NBOOT <- 2000   # réplicas bootstrap
ALPHA <- 0.05   # -> bandas al 95%. Poner 0.10 para bandas al 90%.

NIVEL_BANDAS <- sprintf("%.0f", (1 - ALPHA) * 100)   # se usa en las notas


# ------------------------------------------------------------
# Guardado de gráficos y estilo
# ------------------------------------------------------------

dir_salida <- "graficos_TP2"        # <-- editar a tu ruta
dir.create(dir_salida, showWarnings = FALSE)

guardar_graf <- function(g, nombre, width = 11, height = 3.6) {
  ggsave(file.path(dir_salida, paste0(nombre, ".pdf")), g,
         width = width, height = height, device = cairo_pdf)
  ggsave(file.path(dir_salida, paste0(nombre, ".png")), g,
         width = width, height = height, dpi = 300, bg = "white")
}

tema_paper <- theme_bw(base_size = 11, base_family = "serif") +
  theme(
    plot.title       = element_text(size = 12, face = "bold", hjust = 0, color = "black"),
    strip.background = element_rect(fill = "grey92", color = NA),
    strip.text       = element_text(size = 13, face = "bold", color = "black"),
    panel.grid       = element_blank(),
    panel.border     = element_blank(),
    axis.line        = element_line(color = "black", linewidth = 0.4),
    axis.ticks       = element_line(color = "black"),
    axis.text        = element_text(color = "black", size = 10),
    axis.title       = element_text(size = 11, color = "black")
  )


# ============================================================
# PUNTO 1 (a) — GRÁFICOS DE DOS EJES
#
# No pueden hacerse con facet_wrap porque cada panel necesita su
# propio eje secundario, así que se combinan con patchwork.
# ============================================================

otras <- c("ipc", "ctot", "lembi", "exp_inf", "exp_cam",
           "mp_fed", "oil_shock", "ebp")

etiquetas <- c(
  ipc       = "IPC desest. (100·log)",
  ctot      = "Términos de intercambio (100·log)",
  lembi     = "EMBI (100·log de pbs)",
  exp_inf   = "Exp. inflación 12m",
  exp_cam   = "Exp. TC 12m (100·log)",
  mp_fed    = "Shock MP Fed (JK)",
  oil_shock = "Shock oferta petróleo (BH)",
  ebp       = "Excess Bond Premium"
)

# Los shocks se grafican como barras: una línea de ruido blanco contra
# un TCN tendencial no se lee.
como_barra <- c("mp_fed", "oil_shock", "ebp")

graf_2ejes <- function(v) {
  
  c1 <- base$tcn
  c2 <- base[[v]]
  
  r1 <- range(c1, na.rm = TRUE)
  r2 <- range(c2, na.rm = TRUE)
  
  a <- diff(r1) / diff(r2)
  b <- r1[1] - a * r2[1]
  
  df <- data.frame(date_m = base$date_m, tcn = c1, otra_esc = a * c2 + b)
  
  # Se reportan las dos correlaciones: la de niveles es espuria cuando
  # ambas series son I(1), y el contraste entre ambas es informativo.
  r_niv <- cor(c1, c2, use = "complete.obs")
  r_dif <- cor(c(NA, diff(c1)), c(NA, diff(c2)), use = "complete.obs")
  
  capa2 <- if (v %in% como_barra) {
    geom_col(aes(y = otra_esc - b), fill = "firebrick",
             width = 20, alpha = 0.8, position = position_nudge(y = b))
  } else {
    geom_line(aes(y = otra_esc), color = "firebrick", linewidth = 0.8)
  }
  
  ggplot(df, aes(x = date_m)) +
    capa2 +
    geom_line(aes(y = tcn), color = "steelblue4", linewidth = 0.8) +
    scale_y_continuous(
      name     = "TCN (100·log R$/US$)",
      sec.axis = sec_axis(~ (. - b) / a, name = etiquetas[v])
    ) +
    scale_x_date(date_breaks = "5 years", date_labels = "%Y") +
    labs(x = NULL,
         title = sprintf("%s   (r niv. = %.2f ; r dif. = %.2f)",
                         etiquetas[v], r_niv, r_dif)) +
    tema_paper +
    theme(
      axis.title.y       = element_text(color = "steelblue4"),
      axis.text.y        = element_text(color = "steelblue4"),
      axis.title.y.right = element_text(color = "firebrick"),
      axis.text.y.right  = element_text(color = "firebrick")
    )
}

grupos <- split(otras, ceiling(seq_along(otras) / 2))

for (i in seq_along(grupos)) {
  g <- wrap_plots(lapply(grupos[[i]], graf_2ejes), ncol = 2)
  print(g)
  guardar_graf(g, paste0("p1_panel_", i), width = 11, height = 3.6)
}


# ============================================================
# PUNTO 1 (b) — ORDEN DE INTEGRACIÓN
#     ADF  (ur.df)  : H0 = raíz unitaria.
#     KPSS (ur.kpss): H0 = estacionariedad.
#     Conclusión: I(0) -> NIVELES ; I(1) -> DIFERENCIAS
# ============================================================

adf_rechaza <- function(s) {
  s <- na.omit(as.numeric(s))
  t <- ur.df(s, type = "drift", lags = 12, selectlags = "AIC")
  t@teststat[1] < t@cval[1, "5pct"]
}

kpss_rechaza <- function(s) {
  s <- na.omit(as.numeric(s))
  t <- ur.kpss(s, type = "mu", lags = "long")
  t@teststat[1] > t@cval[1, "5pct"]
}

orden_int <- function(x, nombre) {
  
  a0 <- adf_rechaza(x);  k0 <- kpss_rechaza(x)
  dx <- diff(na.omit(as.numeric(x)))
  a1 <- adf_rechaza(dx); k1 <- kpss_rechaza(dx)
  
  veredicto <- if (a0 && !k0) {
    "I(0)"                                 # ambos coinciden: estacionaria
  } else if (!a0 && k0) {
    if (a1 && !k1) "I(1)" else "I(2)?"     # ambos coinciden: raíz unitaria
  } else {
    "AMBIGUO"                              # tests en conflicto
  }
  
  tibble(
    Serie      = nombre,
    N          = sum(!is.na(x)),
    ADF_nivel  = ifelse(a0, "rechaza RU",     "no rechaza"),
    KPSS_nivel = ifelse(k0, "rechaza estac.", "no rechaza"),
    ADF_dif    = ifelse(a1, "rechaza RU",     "no rechaza"),
    KPSS_dif   = ifelse(k1, "rechaza estac.", "no rechaza"),
    Conclusion = veredicto
  )
}

vars_test <- c("tcn", "ipc", "ctot", "lembi", "exp_inf",
               "exp_cam", "mp_fed", "oil_shock", "ebp")

tabla_I <- map_dfr(vars_test, ~ orden_int(base[[.x]], .x))
print(tabla_I, n = 20)

# Regla de uso en el VAR
tabla_I %>%
  mutate(Especificacion = case_when(
    Conclusion == "I(0)" ~ "Niveles",
    Conclusion == "I(1)" ~ "Diferencias",
    TRUE                 ~ "Decidir por teoría"
  )) %>%
  dplyr::select(Serie, Conclusion, Especificacion) %>%
  print(n = 20)

# ---- Tabla de estadísticos crudos (anexo) ------------------

stat_ur <- function(x, nombre) {
  x    <- na.omit(as.numeric(x))
  adf  <- ur.df(x, type = "drift", lags = 12, selectlags = "AIC")
  kpss <- ur.kpss(x, type = "mu", lags = "long")
  tibble(
    Serie = nombre,
    N     = length(x),
    ADF   = sprintf("%.2f%s", adf@teststat[1],
                    ifelse(adf@teststat[1]  < adf@cval[1, "5pct"],  "*", "")),
    KPSS  = sprintf("%.3f%s", kpss@teststat[1],
                    ifelse(kpss@teststat[1] > kpss@cval[1, "5pct"], "*", ""))
  )
}

tabla_stats <- bind_rows(
  map_dfr(vars_test, ~ stat_ur(base[[.x]], .x)),
  map_dfr(vars_test, ~ stat_ur(diff(na.omit(base[[.x]])), paste0("D.", .x)))
)
print(tabla_stats, n = 20)

tabla_stats %>%
  kbl(booktabs = TRUE, format = "latex", align = "lccc",
      caption = "Tests de raíz unitaria") %>%
  kable_styling(latex_options = "hold_position") %>%
  pack_rows("Niveles", 1, length(vars_test)) %>%
  pack_rows("Primeras diferencias", length(vars_test) + 1, 2 * length(vars_test)) %>%
  footnote(
    general = paste("(*) rechazo al 5\\%. ADF: H0 = raíz unitaria (con drift,",
                    "rezagos por AIC, máx. 12). KPSS: H0 = estacionariedad (ventana larga)."),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_ur.tex"))


# ============================================================
# PUNTO 2 — VAR BIVARIADO, IRF ACUMULADAS Y ERPT
# ============================================================

# ---- 1) Datos: ambas series son I(1) -> primeras diferencias

d_var2 <- base %>%
  transmute(date_m,
            d_tcn = c(NA, diff(tcn)),
            d_ipc = c(NA, diff(ipc))) %>%
  drop_na()

Y <- d_var2 %>% dplyr::select(d_tcn, d_ipc) %>% as.matrix()

cat("\nMuestra:", format(min(d_var2$date_m), "%Y-%m"), "a",
    format(max(d_var2$date_m), "%Y-%m"), " N =", nrow(Y), "\n")

# ---- 2) Selección de rezagos y estimación

sel <- VARselect(Y, lag.max = 12, type = "const")
print(sel$selection)      # AIC/FPE/HQ -> 4 ; SC(BIC) -> 1
print(round(t(sel$criteria), 3))

p    <- 4
var2 <- VAR(Y, p = p, type = "const")

# Diagnósticos: p = 4 es la única especificación con residuos limpios
for (pp in 1:4) {
  st <- serial.test(VAR(Y, p = pp, type = "const"),
                    lags.pt = 12, type = "PT.adjusted")
  cat("p =", pp, " Portmanteau p-valor =", round(st$serial$p.value, 4), "\n")
}

cat("Módulo de las raíces (deben ser < 1):\n"); print(roots(var2))
cat("Correlación entre innovaciones:", round(cor(resid(var2))[1, 2], 3), "\n")

# ---- 3) IRF acumuladas y ERPT
#
#   El VAR está en diferencias, de modo que la respuesta del NIVEL de
#   cada variable es la suma acumulada de la respuesta de su diferencia.
#   Con orden de Cholesky (D.TCN, D.IPC) el shock 1 es el cambiario:
#   puede afectar al IPC dentro del mes, pero un shock de precios no
#   afecta al TCN contemporáneamente.
#
#   Se normaliza el shock a 1% de depreciación en impacto, de modo que
#   el ERPT se lee directamente como fracción del traspaso.

cirf_erpt <- function(fit, H) {
  Ps    <- Psi(fit, nstep = H)     # IRF ortogonalizadas: [variable, shock, h]
  c_tcn <- cumsum(Ps[1, 1, ])
  c_ipc <- cumsum(Ps[2, 1, ])
  n     <- c_tcn[1]                # tamaño del shock en el impacto
  list(tcn   = c_tcn / n,
       ipc   = c_ipc / n,
       erpt  = c_ipc / c_tcn,      # la normalización se cancela en el ratio
       size  = n)
}

est <- cirf_erpt(var2, H)

cat("\nTamaño del shock de 1 d.e.:", round(est$size, 3),
    "puntos log de depreciación\n")

# Respuesta NO acumulada, para documentar el rezago del traspaso
cat("Respuesta de la inflación mensual por horizonte (h = 0..6):\n")
print(round(Psi(var2, nstep = 6)[2, 1, ] / est$size, 4))

# ---- 4) Bandas por bootstrap recursivo
#
#   El ERPT es un cociente de dos estimadores correlacionados, así que
#   el delta method es poco confiable. Se recomputa el cociente completo
#   en cada réplica y se toman percentiles de su distribución.

boot_bandas <- function(fit, H, nboot = NBOOT, alpha = ALPHA, seed = 1234) {
  
  set.seed(seed)
  Yd <- as.matrix(fit$y)
  p  <- fit$p
  K  <- fit$K
  Tt <- nrow(Yd)
  
  U  <- resid(fit)
  B  <- Bcoef(fit)                          # K x (K*p + 1), const al final
  A  <- B[, 1:(K * p), drop = FALSE]
  cc <- B[, K * p + 1]
  
  M_tcn <- M_ipc <- M_erpt <- matrix(NA_real_, nboot, H + 1)
  
  for (b in seq_len(nboot)) {
    
    Us <- U[sample(nrow(U), replace = TRUE), , drop = FALSE]
    Ys <- matrix(NA_real_, Tt, K)
    Ys[1:p, ] <- Yd[1:p, ]
    
    for (t in (p + 1):Tt) {
      # Bcoef ordena los regresores y1.l1, y2.l1, y1.l2, y2.l2, ...
      lagvec  <- as.vector(t(Ys[(t - 1):(t - p), , drop = FALSE]))
      Ys[t, ] <- cc + A %*% lagvec + Us[t - p, ]
    }
    colnames(Ys) <- colnames(Yd)
    
    fb <- try(VAR(Ys, p = p, type = "const"), silent = TRUE)
    if (inherits(fb, "try-error")) next
    
    r <- cirf_erpt(fb, H)
    M_tcn[b, ]  <- r$tcn
    M_ipc[b, ]  <- r$ipc
    M_erpt[b, ] <- r$erpt
  }
  
  qs <- function(M) tibble(
    lo = apply(M, 2, quantile, probs = alpha / 2,     na.rm = TRUE),
    hi = apply(M, 2, quantile, probs = 1 - alpha / 2, na.rm = TRUE)
  )
  list(tcn = qs(M_tcn), ipc = qs(M_ipc), erpt = qs(M_erpt))
}

bd <- boot_bandas(var2, H)

# ---- 5) Resultados en formato tidy

res_p2 <- bind_rows(
  tibble(h = 0:H, serie = "CIRF TCN", est = est$tcn,  lo = bd$tcn$lo,  hi = bd$tcn$hi),
  tibble(h = 0:H, serie = "CIRF IPC", est = est$ipc,  lo = bd$ipc$lo,  hi = bd$ipc$hi),
  tibble(h = 0:H, serie = "ERPT",     est = est$erpt, lo = bd$erpt$lo, hi = bd$erpt$hi)
) %>%
  mutate(serie = factor(serie, levels = c("CIRF TCN", "CIRF IPC", "ERPT")))

print(res_p2 %>% filter(h %in% c(0, 3, 4, 5, 6, 12, 18, 24)), n = 30)

# ---- 6) Gráfico
#
#   Con facet_wrap en lugar de patchwork: las tiras de facetas
#   resuelven solas el problema de los títulos superpuestos.

df_irf <- res_p2 %>%
  mutate(response = factor(
    serie,
    levels = c("CIRF TCN", "CIRF IPC", "ERPT"),
    labels = c("Respuesta acumulada del TCN",
               "Respuesta acumulada del IPC",
               "ERPT incondicional")
  ))

g_irf <- ggplot(df_irf, aes(h, est)) +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  facet_wrap(~ response, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, H, 6)) +
  labs(x = "Meses desde el shock", y = "Respuesta acumulada (nivel)") +
  tema_paper

print(g_irf)
guardar_graf(g_irf, "p2_irf_erpt")

# ---- 7) Tabla para el informe

tabla_p2 <- res_p2 %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
  mutate(txt = sprintf("%.3f [%.3f, %.3f]", est, lo, hi)) %>%
  dplyr::select(h, serie, txt) %>%
  pivot_wider(names_from = serie, values_from = txt)

print(tabla_p2)

tabla_p2 %>%
  kbl(booktabs = TRUE, format = "latex", align = "cccc",
      col.names = c("$h$", "CIRF TCN", "CIRF IPC", "ERPT"),
      caption = "IRF acumuladas y ERPT incondicional", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(
    general = paste0(
      "VAR(", p, ") en primeras diferencias de TCN e IPC, orden de Cholesky ",
      "(TCN, IPC). Shock normalizado a 1\\% de depreciación en impacto. ",
      "Intervalos al ", NIVEL_BANDAS, "\\% por bootstrap recursivo con ",
      NBOOT, " réplicas. Muestra: ",
      format(min(d_var2$date_m), "%Ym%m"), "--",
      format(max(d_var2$date_m), "%Ym%m"), " (N = ", nrow(Y), ")."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto2.tex"))

# ---- 8) Robustez: p = 1 (elegido por BIC)

var2_p1 <- VAR(Y, p = 1, type = "const")
est_p1  <- cirf_erpt(var2_p1, H)

tibble(h = 0:H, erpt_p4 = est$erpt, erpt_p1 = est_p1$erpt) %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
  print()

# ============================================================
# PUNTO 3
#
# Esquema de dos pasos:
#   Paso 1: VAR(4) como el del punto 2 -> shock estructural cambiario
#   Paso 2: x_{t+h} - x_{t-1} = a_h + b_h * u_t + g_h(L) Y_{t-1} + w_{t+h}
#
# Requiere: el script de los puntos 1-2 ya corrido (base, tema_paper,
# guardar_graf, dir_salida, H, res_p2).
# ============================================================

# ------------------------------------------------------------
# 1) Datos y shock estructural
# ------------------------------------------------------------

lev_tcn <- base$tcn
lev_ipc <- base$ipc

D <- cbind(d_tcn = diff(lev_tcn), d_ipc = diff(lev_ipc))
colnames(D) <- c("d_tcn", "d_ipc")
N <- nrow(D)     # 311
p <- 4

var_p3 <- VAR(D, p = p, type = "const")

# El shock estructural cambiario. Con Cholesky (TCN, IPC) la primera
# fila de A0^{-1} es (b11, 0), de modo que la innovación de la ecuación
# del TCN es b11 veces el shock estructural. Se usa la innovación
# directamente: al estar medida en puntos log, b_h queda normalizado a
# una depreciación de 1% en el impacto.
u <- rep(NA_real_, N)
u[(p + 1):N] <- resid(var_p3)[, 1]

cat("Shocks inferidos:", sum(!is.na(u)), " sd =", round(sd(u, na.rm = TRUE), 3), "\n")

# ------------------------------------------------------------
# 2) Local projections
#
#    Indexación: D[j] = lev[j+1] - lev[j], de modo que el shock en la
#    fila j de D corresponde al mes t = j+1 del vector de niveles.
#    El lado izquierdo x_{t+h} - x_{t-1} es entonces lev[j+1+h] - lev[j].
#
#    Nota importante: como el lado izquierdo ya es una diferencia de
#    NIVELES entre t+h y t-1, el coeficiente b_h ES la respuesta
#    acumulada. No hace falta ningún cumsum, a diferencia del VAR.
# ------------------------------------------------------------

controles <- function(js) {
  X <- lapply(1:p, function(l) D[js - l, , drop = FALSE])
  X <- do.call(cbind, X)
  colnames(X) <- paste0(rep(c("d_tcn", "d_ipc"), p), ".l", rep(1:p, each = 2))
  X
}

lp_h <- function(h) {
  
  js <- (p + 1):(N - h)
  W  <- controles(js)
  
  y_tcn <- lev_tcn[js + 1 + h] - lev_tcn[js]
  y_ipc <- lev_ipc[js + 1 + h] - lev_ipc[js]
  shock <- u[js]
  
  dat <- data.frame(y_tcn, y_ipc, shock, W)
  
  # --- IRF: un OLS por variable, con errores HAC.
  #     Los residuos w_{t+h} están autocorrelacionados por construcción
  #     (ventanas superpuestas siguen un MA(h)), de modo que HAC no es
  #     opcional sino necesario.
  fit_tcn <- lm(y_tcn ~ . - y_ipc, data = dat)
  fit_ipc <- lm(y_ipc ~ . - y_tcn, data = dat)
  
  hac <- function(m) sqrt(diag(NeweyWest(m, lag = h + 1, prewhite = FALSE)))["shock"]
  
  # --- ERPT: el cociente b_ipc / b_tcn se obtiene directamente como el
  #     coeficiente de 2SLS de regresar el cambio del IPC sobre el del
  #     TCN, instrumentando con el shock y manteniendo los mismos
  #     controles. Es algebraicamente idéntico al cociente, pero además
  #     entrega inferencia sobre él sin necesidad de método delta.
  iv <- ivreg(y_ipc ~ y_tcn + W | shock + W, data = dat)
  se_iv <- sqrt(diag(NeweyWest(iv, lag = h + 1, prewhite = FALSE)))["y_tcn"]
  
  tibble(
    h        = h,
    N        = length(js),
    b_tcn    = coef(fit_tcn)["shock"],  se_tcn  = hac(fit_tcn),
    b_ipc    = coef(fit_ipc)["shock"],  se_ipc  = hac(fit_ipc),
    erpt     = coef(iv)["y_tcn"],       se_erpt = se_iv
  )
}

lp_res <- map_dfr(0:H, lp_h)

print(lp_res %>% filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# Chequeo: b_tcn en h = 0 debe dar exactamente 1 por Frisch-Waugh-Lovell
cat("\nb_tcn en h=0 (debe ser 1):", round(lp_res$b_tcn[1], 6), "\n")

# ------------------------------------------------------------
# 3) Tidy con bandas al 95% y comparación contra el VAR
# ------------------------------------------------------------

z <- qnorm(1 - ALPHA / 2)

lp_tidy <- bind_rows(
  lp_res %>% transmute(h, serie = "CIRF TCN", est = b_tcn,
                       lo = b_tcn - z * se_tcn,  hi = b_tcn + z * se_tcn),
  lp_res %>% transmute(h, serie = "CIRF IPC", est = b_ipc,
                       lo = b_ipc - z * se_ipc,  hi = b_ipc + z * se_ipc),
  lp_res %>% transmute(h, serie = "ERPT",     est = erpt,
                       lo = erpt  - z * se_erpt, hi = erpt  + z * se_erpt)
) %>%
  mutate(metodo = "Local projections")

comp <- bind_rows(
  lp_tidy,
  res_p2 %>% mutate(metodo = "VAR")            # del punto 2
) %>%
  mutate(
    serie  = factor(serie, levels = c("CIRF TCN", "CIRF IPC", "ERPT"),
                    labels = c("Respuesta acumulada del TCN",
                               "Respuesta acumulada del IPC",
                               "ERPT")),
    metodo = factor(metodo, levels = c("VAR", "Local projections"))
  )

# ------------------------------------------------------------
# 4) Gráfico comparativo
# ------------------------------------------------------------

g_comp <- ggplot(comp, aes(h, est, color = metodo, fill = metodo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ serie, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, H, 6)) +
  scale_color_manual(values = c("VAR" = "steelblue4", "Local projections" = "firebrick")) +
  scale_fill_manual(values  = c("VAR" = "steelblue",  "Local projections" = "firebrick")) +
  labs(x = "Meses desde el shock", y = "Respuesta acumulada (nivel)",
       color = NULL, fill = NULL) +
  tema_paper +
  theme(legend.position = "bottom")

print(g_comp)
guardar_graf(g_comp, "p3_var_vs_lp", width = 11, height = 4.2)

# ------------------------------------------------------------
# 5) Tabla comparativa para el informe
# ------------------------------------------------------------

tabla_p3 <- comp %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24), serie == "ERPT") %>%
  mutate(txt = sprintf("%.3f [%.3f, %.3f]", est, lo, hi)) %>%
  dplyr::select(h, metodo, txt) %>%
  pivot_wider(names_from = metodo, values_from = txt)

print(tabla_p3)

tabla_p3 %>%
  kbl(booktabs = TRUE, format = "latex", align = "ccc",
      col.names = c("$h$", "VAR", "Local projections"),
      caption = "ERPT: VAR frente a local projections", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(
    general = paste0(
      "Intervalos al ", NIVEL_BANDAS, "\\%. Para el VAR, bootstrap recursivo con ",
      NBOOT, " réplicas. Para las local projections, errores estándar de ",
      "Newey-West con $h+1$ rezagos, obtenidos de la estimación por variables ",
      "instrumentales descrita en el texto."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto3.tex"))

# ------------------------------------------------------------
# 6) ¿Son estadísticamente distintos VAR y LP?
#    Se verifica si el punto estimado del VAR cae dentro de la banda de LP.
# ------------------------------------------------------------

comp %>%
  filter(serie == "ERPT") %>%
  dplyr::select(h, metodo, est, lo, hi) %>%
  pivot_wider(names_from = metodo, values_from = c(est, lo, hi)) %>%
  mutate(var_dentro_banda_lp =
           `est_VAR` >= `lo_Local projections` & `est_VAR` <= `hi_Local projections`) %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
  print()

# ============================================================
# PUNTO 4
# LP asimétricas según el signo del shock cambiario
#
# Especificación:
#
#   D_t = 1{u_t > 0}                       (dummy de signo)
#   u_pos = D_t * u_t                      (dummy x shock)
#   u_neg = (1 - D_t) * u_t
#
#   x_{t+h} - x_{t-1} = a_h + b_h^+ u_pos + b_h^- u_neg
#                       + g_h(L) Y_{t-1} + w_{t+h}
#
# La dummy va INTERACTUADA, no sola: lo que debe cambiar según el
# signo es la pendiente (el traspaso por punto de shock), no la
# constante. La especificación es continua en cero y anida la
# simetría: imponer b^+ = b^- devuelve la ecuación del punto 3.
#
# Requiere el punto 3 corrido (D, u, N, p, H, lev_tcn, lev_ipc,
# controles(), tema_paper, guardar_graf, dir_salida, ALPHA).
# ============================================================


# ------------------------------------------------------------
# 1) Descripción de los shocks por signo
# ------------------------------------------------------------

uu <- u[!is.na(u)]
cat("Shocks:", length(uu),
    "| positivos:", sum(uu > 0), sprintf("(%.1f%%)", 100 * mean(uu > 0)),
    "| negativos:", sum(uu < 0), sprintf("(%.1f%%)", 100 * mean(uu < 0)), "\n")
cat("Media positivos:", round(mean(uu[uu > 0]), 3),
    "| media negativos:", round(mean(uu[uu < 0]), 3),
    "| sd:", round(sd(uu), 3), "\n\n")

# ------------------------------------------------------------
# 2) LP asimétrica
# ------------------------------------------------------------

lp_signo_h <- function(h) {
  
  js <- (p + 1):(N - h)
  W  <- as.data.frame(controles(js))
  ct <- names(W)                                  # nombres de los controles
  
  dum <- as.numeric(u[js] > 0)                    # D_t
  
  dat <- data.frame(
    y_tcn = lev_tcn[js + 1 + h] - lev_tcn[js],
    y_ipc = lev_ipc[js + 1 + h] - lev_ipc[js],
    u_pos = dum * u[js],
    u_neg = (1 - dum) * u[js],
    W
  )
  
  rhs <- paste("u_pos + u_neg +", paste(ct, collapse = " + "))
  
  # --- IRF por signo, con test de Wald de igualdad
  #
  #     En h = 0 la ecuación del TCN tiene ajuste perfecto: por
  #     Frisch-Waugh-Lovell ambos coeficientes valen exactamente 1 y la
  #     suma de cuadrados residual es cero. El test de Wald es
  #     degenerado en ese caso y se devuelve NA.
  est_signo <- function(dep) {
    
    f <- lm(as.formula(paste(dep, "~", rhs)), data = dat)
    
    V <- tryCatch(
      NeweyWest(f, lag = h + 1, prewhite = FALSE),
      error = function(e) {
        k <- names(coef(f))
        matrix(0, length(k), length(k), dimnames = list(k, k))
      }
    )
    
    pv <- tryCatch({
      lh <- linearHypothesis(f, "u_pos = u_neg", vcov. = V, test = "Chisq")
      lh$`Pr(>Chisq)`[2]
    }, error = function(e) NA_real_)
    
    c(bp   = unname(coef(f)["u_pos"]), sep = sqrt(V["u_pos", "u_pos"]),
      bn   = unname(coef(f)["u_neg"]), sen = sqrt(V["u_neg", "u_neg"]),
      pval = pv)
  }
  
  r_tcn <- est_signo("y_tcn")
  r_ipc <- est_signo("y_ipc")
  
  # --- ERPT por signo, por variables instrumentales.
  #
  #     Para el ERPT+ se instrumenta el cambio del TCN con u_pos,
  #     INCLUYENDO u_neg como control exógeno. Sin ese control el
  #     coeficiente de IV deja de coincidir con b^+_ipc / b^+_tcn,
  #     porque los controles W son compartidos entre ambos signos.
  #     El bloque de verificación de más abajo lo confirma.
  erpt_signo <- function(inst, otro) {
    
    fml <- as.formula(paste0(
      "y_ipc ~ y_tcn + ", otro, " + ", paste(ct, collapse = " + "),
      " | ",  inst,  " + ", otro, " + ", paste(ct, collapse = " + ")
    ))
    iv <- ivreg(fml, data = dat)
    V  <- tryCatch(NeweyWest(iv, lag = h + 1, prewhite = FALSE),
                   error = function(e) {
                     k <- names(coef(iv))
                     matrix(NA_real_, length(k), length(k), dimnames = list(k, k))
                   })
    
    # F de primera etapa del instrumento
    f1 <- lm(as.formula(paste0("y_tcn ~ ", inst, " + ", otro, " + ",
                               paste(ct, collapse = " + "))), data = dat)
    V1 <- tryCatch(NeweyWest(f1, lag = h + 1, prewhite = FALSE),
                   error = function(e) {
                     k <- names(coef(f1))
                     matrix(NA_real_, length(k), length(k), dimnames = list(k, k))
                   })
    
    c(est = unname(coef(iv)["y_tcn"]),
      se  = sqrt(V["y_tcn", "y_tcn"]),
      F1  = unname(coef(f1)[inst])^2 / V1[inst, inst])
  }
  
  e_pos <- erpt_signo("u_pos", "u_neg")
  e_neg <- erpt_signo("u_neg", "u_pos")
  
  tibble(
    h = h, N = length(js),
    bp_tcn = r_tcn["bp"], sep_tcn = r_tcn["sep"],
    bn_tcn = r_tcn["bn"], sen_tcn = r_tcn["sen"], p_tcn = r_tcn["pval"],
    bp_ipc = r_ipc["bp"], sep_ipc = r_ipc["sep"],
    bn_ipc = r_ipc["bn"], sen_ipc = r_ipc["sen"], p_ipc = r_ipc["pval"],
    erpt_p = e_pos["est"], se_erpt_p = e_pos["se"], F_p = e_pos["F1"],
    erpt_n = e_neg["est"], se_erpt_n = e_neg["se"], F_n = e_neg["F1"]
  )
}

lp4 <- map_dfr(0:H, lp_signo_h)

# --- Verificación: el IV debe reproducir el cociente exactamente
cat("\nVerificación IV = cociente:\n")
lp4 %>%
  filter(h %in% c(3, 6, 12, 24)) %>%
  transmute(h,
            iv_pos = erpt_p, ratio_pos = bp_ipc / bp_tcn,
            iv_neg = erpt_n, ratio_neg = bn_ipc / bn_tcn) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>%
  print()

cat("\nResultados principales:\n")
print(lp4 %>% filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
        dplyr::select(h, bp_ipc, bn_ipc, p_ipc, erpt_p, erpt_n, F_p, F_n) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

cat("\nRespuesta del TCN por signo:\n")
print(lp4 %>% filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
        dplyr::select(h, bp_tcn, sep_tcn, bn_tcn, sen_tcn, p_tcn) %>%
        mutate(across(where(is.numeric), ~ round(.x, 4))))

# ------------------------------------------------------------
# 3) Gráfico
# ------------------------------------------------------------

z <- qnorm(1 - ALPHA / 2)

df4 <- bind_rows(
  lp4 %>% transmute(h, serie = "Respuesta acumulada del TCN", signo = "Depreciación",
                    est = bp_tcn, lo = bp_tcn - z*sep_tcn, hi = bp_tcn + z*sep_tcn),
  lp4 %>% transmute(h, serie = "Respuesta acumulada del TCN", signo = "Apreciación",
                    est = bn_tcn, lo = bn_tcn - z*sen_tcn, hi = bn_tcn + z*sen_tcn),
  lp4 %>% transmute(h, serie = "Respuesta acumulada del IPC", signo = "Depreciación",
                    est = bp_ipc, lo = bp_ipc - z*sep_ipc, hi = bp_ipc + z*sep_ipc),
  lp4 %>% transmute(h, serie = "Respuesta acumulada del IPC", signo = "Apreciación",
                    est = bn_ipc, lo = bn_ipc - z*sen_ipc, hi = bn_ipc + z*sen_ipc),
  lp4 %>% transmute(h, serie = "ERPT", signo = "Depreciación",
                    est = erpt_p, lo = erpt_p - z*se_erpt_p, hi = erpt_p + z*se_erpt_p),
  lp4 %>% transmute(h, serie = "ERPT", signo = "Apreciación",
                    est = erpt_n, lo = erpt_n - z*se_erpt_n, hi = erpt_n + z*se_erpt_n)
) %>%
  mutate(serie = factor(serie, levels = c("Respuesta acumulada del TCN",
                                          "Respuesta acumulada del IPC", "ERPT")),
         signo = factor(signo, levels = c("Depreciación", "Apreciación")))

g4 <- ggplot(df4, aes(h, est, color = signo, fill = signo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ serie, scales = "free_y") +
  scale_x_continuous(breaks = seq(0, H, 6)) +
  scale_color_manual(values = c("Depreciación" = "firebrick",
                                "Apreciación"  = "steelblue4")) +
  scale_fill_manual(values  = c("Depreciación" = "firebrick",
                                "Apreciación"  = "steelblue")) +
  labs(x = "Meses desde el shock", y = "Respuesta acumulada (nivel)",
       color = NULL, fill = NULL) +
  tema_paper +
  theme(legend.position = "bottom")

print(g4)
guardar_graf(g4, "p4_asimetria_signo", width = 11, height = 4.2)

# ------------------------------------------------------------
# 4) Tabla para el informe
# ------------------------------------------------------------


tabla_p4 <- lp4 %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
  transmute(
    h,
    ERPT_pos = sprintf("%.3f [%.3f, %.3f]", erpt_p,
                       erpt_p - z*se_erpt_p, erpt_p + z*se_erpt_p),
    ERPT_neg = sprintf("%.3f [%.3f, %.3f]", erpt_n,
                       erpt_n - z*se_erpt_n, erpt_n + z*se_erpt_n),
    pval = ifelse(is.na(p_ipc), "---", sprintf("%.3f", p_ipc)),
    Fpos = ifelse(h == 0, "---", sprintf("%.1f", F_p)),
    Fneg = ifelse(h == 0, "---", sprintf("%.1f", F_n))
  )
print(tabla_p4)

tabla_p4 %>%
  kbl(booktabs = TRUE, format = "latex", align = "cccccc",
      col.names = c("$h$", "ERPT$^{+}$", "ERPT$^{-}$",
                    "$p$-valor", "$F^{+}$", "$F^{-}$"),
      caption = "ERPT según el signo del shock cambiario", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(
    general = paste0(
      "Intervalos al ", NIVEL_BANDAS, "\\% con errores de Newey-West y $h+1$ rezagos. ",
      "El $p$-valor corresponde al test de Wald de igualdad entre los coeficientes de ",
      "respuesta del IPC ante shocks positivos y negativos. $F^{+}$ y $F^{-}$ son los ",
      "estadísticos de primera etapa de cada instrumento."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto4.tex"))

# ============================================================
# PUNTO 5
# Histograma de los shocks cambiarios estimados en el inciso 3
# y evaluación de ERPT según el tamaño del shock.
#
# ============================================================

# ------------------------------------------------------------
# 1) Distribución del shock: momentos y clasificación por tamaño
# ------------------------------------------------------------

# Serie de shocks con su fecha. D[j] = lev[j+1]-lev[j] corresponde al mes
# base$date_m[j+1], de modo que u (alineado a D) se fecha con date_m[-1].

u_dat <- tibble(date_m = base$date_m[-1], u = u) %>% drop_na(u)

uu <- u_dat$u                            # 307 shocks estructurales cambiarios
m  <- mean(uu)
s  <- sd(uu)

# Calculo asimetría y curtosis: te dan una medida de por qué el tamaño de los shocks 
# puede estar rompiendo la identificación del erpt

sk <- mean((uu - m)^3) / s^3             # asimetría
ku <- mean((uu - m)^4) / s^4             # curtosis

# "Normal" = dentro de 1 sd ; las colas son los shocks "grandes".
u_dat <- u_dat %>%
  mutate(regimen = cut(u, breaks = c(-Inf, -s, s, Inf),
                       labels = c("Apreciación grande (< -1 sd)",
                                  "Normal (|u| <= 1 sd)",
                                  "Depreciación grande (> +1 sd)")))

tabla_regimen <- u_dat %>% count(regimen) %>% mutate(pct = 100 * n / sum(n))

cat(sprintf("\nShocks: %d | media: %.3f | sd: %.3f | asimetría: %.2f | curtosis: %.2f (exceso %.2f)\n",
            length(uu), m, s, sk, ku, ku - 3))
cat(sprintf("Mínimo: %.2f (%.2f sd)  Máximo: %.2f (%.2f sd)\n",
            min(uu), min(uu) / s, max(uu), max(uu) / s))
print(tabla_regimen)

# Episodios extremos (|u| > 2 sd) para nombrarlos en el texto
cat("\nEpisodios extremos (|u| > 2 sd):\n")
u_dat %>% filter(abs(u) > 2 * s) %>% arrange(desc(abs(u))) %>%
  mutate(u = round(u, 2), en_sd = round(u / s, 2)) %>% print(n = 40)

# --- Tabla de ANEXO: episodios cambiarios extremos (LaTeX) ---
# Descompone el movimiento total del mes en la parte sorpresiva (el shock)
# y la parte predecible por la historia del VAR:  Δtcn = û_t + predecible.
tabla_episodios <- tibble(date_m = base$date_m, d_tcn = c(NA, diff(base$tcn))) %>%
  inner_join(u_dat, by = "date_m") %>%
  filter(abs(u) > 2 * s) %>%
  arrange(desc(abs(u))) %>%
  transmute(
    Mes        = format(date_m, "%Y-%m"),
    u_sd       = sprintf("%+.2f", u / s),
    shock      = sprintf("%+.1f", u),
    predecible = sprintf("%+.1f", d_tcn - u),
    d_tcn_log  = sprintf("%+.1f", d_tcn),
    deva_pct   = sprintf("%+.1f", (exp(d_tcn / 100) - 1) * 100),
    tipo       = ifelse(u > 0, "Depreciación", "Apreciación"))
print(tabla_episodios, n = 40)

tabla_episodios %>%
  kbl(booktabs = TRUE, format = "latex", align = "lcccccl",
      col.names = c("Mes", "$\\hat u_t/\\sigma$", "$\\hat u_t$ (log)", "Predecible",
                    "$\\Delta tcn$ (log)", "Deva. real (\\%)", "Tipo"),
      caption = "Episodios cambiarios extremos: shocks con $|\\hat u_t| > 2$ desvíos estándar",
      label = "episodios", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position", font_size = 8) %>%
  footnote(
    general = paste0(
      "Meses con shock cambiario estructural mayor a 2 desvíos estándar (", sprintf("%.2f", s),
      " puntos log), ordenados por magnitud. $\\hat u_t$ es la depreciación no anticipada; ",
      "\\emph{Predecible} $= \\Delta tcn - \\hat u_t$ es la parte explicada por la historia del VAR; ",
      "ambas en puntos log ($\\times 100$). \"Deva. real\" convierte $\\Delta tcn$ a variación porcentual del mes."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_anexo_episodios.tex"))


# ------------------------------------------------------------
# 2) Histograma
# ------------------------------------------------------------

bw <- 1.5   # ancho de bin (sd ~ 4.6 -> ~25 bins); boundary = 0 => borde en cero

g5 <- ggplot(u_dat, aes(u)) +
  # franjas de cola: territorio de shocks "grandes" (|u| > 1 sd)
  annotate("rect", xmin = -Inf, xmax = -s, ymin = 0, ymax = Inf,
           fill = "grey50", alpha = 0.08) +
  annotate("rect", xmin =  s,   xmax = Inf, ymin = 0, ymax = Inf,
           fill = "grey50", alpha = 0.08) +
  geom_histogram(binwidth = bw, boundary = 0,
                 fill = "steelblue4", color = "white", linewidth = 0.2) +
  # normal de referencia (misma media y sd), escalada a frecuencias
  stat_function(fun = function(x) dnorm(x, m, s) * nrow(u_dat) * bw,
                color = "firebrick", linewidth = 0.8) +
  geom_rug(sides = "b", alpha = 0.25) +
  geom_vline(xintercept = 0, color = "grey40", linewidth = 0.4) +
  geom_vline(xintercept = c(-s, s), linetype = "dashed",
             color = "grey20", linewidth = 0.5) +
  geom_vline(xintercept = c(-2 * s, 2 * s), linetype = "dotted",
             color = "grey50", linewidth = 0.4) +
  annotate("text", x = c(-s, s), y = Inf, label = c("-1 sd", "+1 sd"),
           vjust = 1.6, hjust = c(1.15, -0.15), size = 3, color = "grey20") +
  annotate("text", x = c(-2 * s, 2 * s), y = Inf, label = c("-2 sd", "+2 sd"),
           vjust = 1.6, hjust = c(1.15, -0.15), size = 3, color = "grey50") +
  labs(
    x = "Shock cambiario no anticipado û (puntos log ≈ % de depreciación)",
    y = "Frecuencia (meses)",
    subtitle = sprintf("N = %d ; sd = %.2f ; asimetría = %.2f ; curtosis = %.1f (normal = 3)",
                       nrow(u_dat), s, sk, ku)
    ) +
  tema_paper +
  theme(plot.subtitle = element_text(size = 10, color = "grey30"))

print(g5)
guardar_graf(g5, "p5_histograma_shocks", width = 8, height = 4.5)


# ------------------------------------------------------------
# 3) ERPT según el TAMAÑO del shock (3 regímenes, umbral ±1 sd)
#
#    Se descompone el shock en tres piezas que suman u_t:
#      u_neg = u * 1{u < -s}   (apreciación grande)
#      u_nrm = u * 1{|u| <= s} (normal)
#      u_pos = u * 1{u >  s}   (depreciación grande)
# ------------------------------------------------------------

# Error estándar HAC (Newey-West, h+1 rezagos) robusto a fallos numéricos
nwse <- function(m, hh, nm) {
  V <- tryCatch(NeweyWest(m, lag = hh + 1, prewhite = FALSE), error = function(e) NULL)
  if (is.null(V) || !(nm %in% rownames(V))) return(NA_real_)
  sqrt(V[nm, nm])
}

lp_size_h <- function(h) {
  js <- (p + 1):(N - h)
  W  <- as.data.frame(controles(js)); ct <- names(W)
  ush <- u[js]
  dat <- data.frame(
    y_tcn = lev_tcn[js + 1 + h] - lev_tcn[js],
    y_ipc = lev_ipc[js + 1 + h] - lev_ipc[js],
    u_neg = ush * (ush < -s),
    u_nrm = ush * (abs(ush) <= s),
    u_pos = ush * (ush >  s),
    W)
  rhs <- paste("u_neg + u_nrm + u_pos +", paste(ct, collapse = " + "))
  
  # IRF por régimen + test de Wald de igualdad de los tres coeficientes
  est_reg <- function(dep) {
    f <- lm(as.formula(paste(dep, "~", rhs)), data = dat)
    V <- tryCatch(NeweyWest(f, lag = h + 1, prewhite = FALSE),
                  error = function(e) { k <- names(coef(f)); matrix(NA, length(k), length(k), dimnames = list(k, k)) })
    pv <- tryCatch(linearHypothesis(f, c("u_neg = u_nrm", "u_nrm = u_pos"),
                                    vcov. = V, test = "Chisq")$`Pr(>Chisq)`[2],
                   error = function(e) NA_real_)
    c(bneg = unname(coef(f)["u_neg"]), seneg = unname(sqrt(V["u_neg", "u_neg"])),
      bnrm = unname(coef(f)["u_nrm"]), senrm = unname(sqrt(V["u_nrm", "u_nrm"])),
      bpos = unname(coef(f)["u_pos"]), sepos = unname(sqrt(V["u_pos", "u_pos"])), pval = unname(pv))
  }
  rt <- est_reg("y_tcn"); ri <- est_reg("y_ipc")
  
  # ERPT del régimen: IV con la pieza como instrumento, las otras dos como controles
  erpt_reg <- function(inst, otras) {
    fml <- as.formula(paste0("y_ipc ~ y_tcn + ", paste(otras, collapse = " + "), " + ", paste(ct, collapse = " + "),
                             " | ", inst, " + ", paste(otras, collapse = " + "), " + ", paste(ct, collapse = " + ")))
    iv <- ivreg(fml, data = dat)
    f1 <- lm(as.formula(paste0("y_tcn ~ ", inst, " + ", paste(otras, collapse = " + "), " + ", paste(ct, collapse = " + "))), data = dat)
    V1 <- tryCatch(NeweyWest(f1, lag = h + 1, prewhite = FALSE), error = function(e) NULL)
    F1 <- if (is.null(V1)) NA_real_ else unname(coef(f1)[inst])^2 / V1[inst, inst]
    c(est = unname(coef(iv)["y_tcn"]), se = nwse(iv, h, "y_tcn"), F1 = F1)
  }
  en <- erpt_reg("u_neg", c("u_nrm", "u_pos"))
  e0 <- erpt_reg("u_nrm", c("u_neg", "u_pos"))
  ep <- erpt_reg("u_pos", c("u_neg", "u_nrm"))
  
  tibble(h = h, N = length(js),
         bneg_ipc = ri["bneg"], seneg_ipc = ri["seneg"], bnrm_ipc = ri["bnrm"], senrm_ipc = ri["senrm"],
         bpos_ipc = ri["bpos"], sepos_ipc = ri["sepos"], p_ipc = ri["pval"],
         bneg_tcn = rt["bneg"], seneg_tcn = rt["seneg"],
         bnrm_tcn = rt["bnrm"], senrm_tcn = rt["senrm"],
         bpos_tcn = rt["bpos"], sepos_tcn = rt["sepos"],
         erpt_neg = en["est"], se_erpt_neg = en["se"], F_neg = en["F1"],
         erpt_nrm = e0["est"], se_erpt_nrm = e0["se"], F_nrm = e0["F1"],
         erpt_pos = ep["est"], se_erpt_pos = ep["se"], F_pos = ep["F1"])
}

lp5b <- map_dfr(0:H, lp_size_h)

# Verificación: IV = cociente b_ipc/b_tcn por régimen
cat("\nVerificación IV = cociente (por régimen):\n")
lp5b %>% filter(h %in% c(3, 6, 12, 24)) %>%
  transmute(h, iv_pos = erpt_pos, ratio_pos = bpos_ipc / bpos_tcn,
            iv_nrm = erpt_nrm, ratio_nrm = bnrm_ipc / bnrm_tcn) %>%
  mutate(across(where(is.numeric), ~ round(.x, 5))) %>% print()

cat("\nERPT por tamaño del shock:\n")
print(lp5b %>% filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
        transmute(h, erpt_neg = round(erpt_neg, 3), erpt_nrm = round(erpt_nrm, 3),
                  erpt_pos = round(erpt_pos, 3), p_ipc = round(p_ipc, 3),
                  F_pos = round(F_pos, 1)))

# --- Gráfico: ERPT y respuesta del IPC por régimen ---
df5b <- bind_rows(
  lp5b %>% transmute(h, serie = "Respuesta acumulada del TCN", reg = "Apreciación grande (< -1 sd)",  est = bneg_tcn, lo = bneg_tcn - z*seneg_tcn, hi = bneg_tcn + z*seneg_tcn),
  lp5b %>% transmute(h, serie = "Respuesta acumulada del TCN", reg = "Normal (|u| <= 1 sd)",          est = bnrm_tcn, lo = bnrm_tcn - z*senrm_tcn, hi = bnrm_tcn + z*senrm_tcn),
  lp5b %>% transmute(h, serie = "Respuesta acumulada del TCN", reg = "Depreciación grande (> +1 sd)", est = bpos_tcn, lo = bpos_tcn - z*sepos_tcn, hi = bpos_tcn + z*sepos_tcn),
  lp5b %>% transmute(h, serie = "Respuesta acumulada del IPC", reg = "Apreciación grande (< -1 sd)",  est = bneg_ipc, lo = bneg_ipc - z*seneg_ipc, hi = bneg_ipc + z*seneg_ipc),
  lp5b %>% transmute(h, serie = "Respuesta acumulada del IPC", reg = "Normal (|u| <= 1 sd)",          est = bnrm_ipc, lo = bnrm_ipc - z*senrm_ipc, hi = bnrm_ipc + z*senrm_ipc),
  lp5b %>% transmute(h, serie = "Respuesta acumulada del IPC", reg = "Depreciación grande (> +1 sd)", est = bpos_ipc, lo = bpos_ipc - z*sepos_ipc, hi = bpos_ipc + z*sepos_ipc),
  lp5b %>% transmute(h, serie = "ERPT", reg = "Apreciación grande (< -1 sd)",  est = erpt_neg, lo = erpt_neg - z*se_erpt_neg, hi = erpt_neg + z*se_erpt_neg),
  lp5b %>% transmute(h, serie = "ERPT", reg = "Normal (|u| <= 1 sd)",          est = erpt_nrm, lo = erpt_nrm - z*se_erpt_nrm, hi = erpt_nrm + z*se_erpt_nrm),
  lp5b %>% transmute(h, serie = "ERPT", reg = "Depreciación grande (> +1 sd)", est = erpt_pos, lo = erpt_pos - z*se_erpt_pos, hi = erpt_pos + z*se_erpt_pos)
) %>% mutate(
  serie = factor(serie, levels = c("Respuesta acumulada del TCN", "Respuesta acumulada del IPC", "ERPT")),
  reg   = factor(reg, levels = c("Apreciación grande (< -1 sd)", "Normal (|u| <= 1 sd)", "Depreciación grande (> +1 sd)")))

g5b <- ggplot(df5b, aes(h, est, color = reg, fill = reg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.9) +
  facet_wrap(~ serie, scales = "free_y") + scale_x_continuous(breaks = seq(0, H, 6)) +
  scale_color_manual(values = c("steelblue4", "grey40", "firebrick")) +
  scale_fill_manual(values  = c("steelblue4", "grey60", "firebrick")) +
  labs(x = "Meses desde el shock", y = "Respuesta acumulada (nivel)", color = NULL, fill = NULL) +
  tema_paper + theme(legend.position = "bottom")

print(g5b)
guardar_graf(g5b, "p5b_erpt_tamano", width = 11, height = 4.4)

# --- Tabla para el informe ---
tabla_p5b <- lp5b %>%
  filter(h %in% c(0, 3, 6, 12, 18, 24)) %>%
  transmute(
    h,
    `ERPT$^{-}$ (apr.)`  = sprintf("%.3f [%.3f, %.3f]", erpt_neg, erpt_neg - z*se_erpt_neg, erpt_neg + z*se_erpt_neg),
    `ERPT$^{0}$ (norm.)` = sprintf("%.3f [%.3f, %.3f]", erpt_nrm, erpt_nrm - z*se_erpt_nrm, erpt_nrm + z*se_erpt_nrm),
    `ERPT$^{+}$ (depr.)` = sprintf("%.3f [%.3f, %.3f]", erpt_pos, erpt_pos - z*se_erpt_pos, erpt_pos + z*se_erpt_pos),
    `$p$-valor` = ifelse(is.na(p_ipc), "---", sprintf("%.3f", p_ipc)),
    `$F^{-}$`   = ifelse(h == 0, "---", sprintf("%.1f", F_neg)),
    `$F^{0}$`   = ifelse(h == 0, "---", sprintf("%.1f", F_nrm)),
    `$F^{+}$`   = ifelse(h == 0, "---", sprintf("%.1f", F_pos))
  )
print(tabla_p5b)

tabla_p5b %>%
  kbl(booktabs = TRUE, format = "latex", align = "cccccccc",
      col.names = c("$h$", "ERPT$^{-}$ (apr.)", "ERPT$^{0}$ (norm.)", "ERPT$^{+}$ (depr.)",
                    "$p$-valor", "$F^{-}$", "$F^{0}$", "$F^{+}$"),
      caption = "ERPT según el tamaño del shock cambiario", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position", font_size = 8) %>%
  footnote(
    general = paste0(
      "Regímenes definidos por el umbral de $\\pm 1$ desvío estándar del shock. Intervalos al ",
      NIVEL_BANDAS, "\\% con errores de Newey-West y $h+1$ rezagos. El $p$-valor corresponde al ",
      "test de Wald de igualdad de los tres coeficientes de respuesta del IPC. $F^{-}$, $F^{0}$ y ",
      "$F^{+}$ son los estadísticos de primera etapa de los regímenes de apreciación grande, normal ",
      "y depreciación grande, respectivamente; el umbral convencional para descartar instrumento débil es 10."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto5.tex"))

# ============================================================
# PUNTO 6
# ERPT CONDICIONAL — VAR de 7 variables, 6 shocks estructurales
#
# Paso 1: VAR sobre (Fed, petróleo, EBP, TI, EMBI, TCN, IPC) con ese
#         orden como identificación recursiva (Cholesky). Se infieren los
#         primeros 6 shocks estructurales (se excluye el del IPC).
# Paso 2: para cada shock, ERPT por LP+IV como en el punto 3.
#
# Acá se recuperan los shocks estructurales ortogonalizando. En el p2 no hacía falta
# porque tenías sólo dos variables
# ============================================================

# ------------------------------------------------------------
# 1) Sistema de 7 variables (mezcla niveles/diferencias por orden de int.)
#    I(0): Fed, petróleo, EBP -> nivel. I(1): TI, EMBI, TCN, IPC -> diff.
#    La muestra la ata el EMBI (termina 2024m7).
# ------------------------------------------------------------
dat7 <- base %>%
  transmute(date_m,
            fed = mp_fed, oil = oil_shock, ebp = ebp,
            dctot = ctot - lag(ctot), dembi = lembi - lag(lembi),
            dtcn = tcn - lag(tcn), dipc = ipc - lag(ipc)) %>%
  drop_na()

Y7 <- dat7 %>% dplyr::select(fed, oil, ebp, dctot, dembi, dtcn, dipc) %>% as.matrix()

cat("\nPunto 6 - muestra VAR7:", format(min(dat7$date_m), "%Y-%m"), "a",
    format(max(dat7$date_m), "%Y-%m"), " N =", nrow(Y7), "\n")

# Selección de rezagos + diagnóstico de residuos
selv7 <- VARselect(Y7, lag.max = 8, type = "const")
print(selv7$selection)   # AIC/HQ/BIC -> 1
for (pp in 1:4) {
  st <- serial.test(VAR(Y7, p = pp, type = "const"), lags.pt = 18, type = "PT.adjusted")
  cat("p7 =", pp, " Portmanteau p =", round(st$serial$p.value, 4), "\n")
}

# Portmanteau (H0: no autocorrelación de residuos) para cada rezago
pv7 <- sapply(1:8, function(pp)
  serial.test(VAR(Y7, p = pp, type = "const"), lags.pt = 18, type = "PT.adjusted")$serial$p.value)
for (pp in 1:8) cat("p7 =", pp, " Portmanteau p =", round(pv7[pp], 4), "\n")

# Los criterios eligen 1, pero con 1 rezago los residuos quedan
# autocorrelacionados (p < 0.01). Ningún rezago los blanquea del todo
# (habitual en VAR mensuales de alta dimensión); 2 es el que más reduce
# la autocorrelación. Para la LP, el VAR es sólo un dispositivo que genera
# instrumentos, no el objeto de inferencia.
p7 <- 2L

crit7 <- as.data.frame(t(selv7$criteria))[, c("AIC(n)", "HQ(n)", "SC(n)")]
fmt_min <- function(x) { s <- sprintf("%.3f", x); s[which.min(x)] <- paste0("\\textbf{", s[which.min(x)], "}"); s }
tabla_lags7 <- tibble(
  Rezagos = 1:8,
  AIC = fmt_min(crit7[["AIC(n)"]]),
  HQ  = fmt_min(crit7[["HQ(n)"]]),
  BIC = fmt_min(crit7[["SC(n)"]]),
  `Portmanteau (18)` = sprintf("%.3f", pv7)
)
print(tabla_lags7)

tabla_lags7 %>%
  kbl(booktabs = TRUE, format = "latex", align = "ccccc",
      col.names = c("Rezagos", "AIC", "HQ", "BIC", "Portmanteau (18)"),
      caption = "Selección de rezagos del VAR de 7 variables", label = "varselect7", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position") %>%
  footnote(
    general = paste("Menor valor = rezago óptimo (en negrita). BIC = criterio de Schwarz,",
                    "HQ = Hannan-Quinn. La última columna reporta el $p$-valor del test de Portmanteau",
                    "ajustado a 18 rezagos (H0: ausencia de autocorrelación). Los criterios seleccionan",
                    "un rezago, pero éste deja residuos autocorrelacionados; se adopta $p=2$, la",
                    "especificación más parsimoniosa que mitiga la autocorrelación."),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_varselect_p6.tex"))

var7 <- VAR(Y7, p = p7, type = "const")
cat("¿raíces del VAR7 < 1?", all(roots(var7) < 1),
    " | máx:", round(max(roots(var7)), 3), "\n")

# ------------------------------------------------------------
# 2) Shocks estructurales (Cholesky en el orden dado)
# ------------------------------------------------------------
E7   <- resid(var7)
Sig7 <- cov(E7)
P7   <- t(chol(Sig7))                        # triangular inferior: Σ = P7 P7'
U7   <- E7 %*% t(solve(P7))                  # u_t = P7^{-1} e_t (var. unitaria)
colnames(U7) <- c("fed", "oil", "ebp", "ctot", "embi", "tcn", "ipc")
cat("cov(U7) ~ I ? diag:", round(diag(cov(U7)), 3),
    "| máx |fuera diag|:", round(max(abs(cov(U7) - diag(7))), 4), "\n")

shocks7 <- as_tibble(U7) %>% mutate(date_m = dat7$date_m[(p7 + 1):nrow(dat7)])
shocks7$bi <- match(shocks7$date_m, base$date_m)   # fila en 'base' de cada shock

# Controles del paso 2: rezagos del vector del VAR de 7 variables
rows7 <- (p7 + 1):nrow(dat7)
W7 <- do.call(cbind, lapply(1:p7, function(l) Y7[rows7 - l, , drop = FALSE]))
colnames(W7) <- paste0(rep(colnames(Y7), p7), ".l", rep(1:p7, each = ncol(Y7)))

shock_names <- c("fed", "oil", "ebp", "ctot", "embi", "tcn")
lab_shock   <- c(fed = "Shock MP Fed", oil = "Shock petróleo", ebp = "EBP",
                 ctot = "Términos interc.", embi = "EMBI", tcn = "Cambiario (TCN)")

# ------------------------------------------------------------
# 3) LP+IV del ERPT para cada shock
#    Instrumento = shock estructural k ; controles = rezagos del VAR7.
#    El primer estadístico F mide la relevancia: sólo los shocks que
#    mueven fuerte al TCN identifican bien el ERPT.
# ------------------------------------------------------------

lp6_hk <- function(h, k) {
  bi <- shocks7$bi
  ok <- (bi + h) <= nrow(base) & (bi - 1) >= 1
  dd <- data.frame(y_tcn = base$tcn[bi + h] - base$tcn[bi - 1],
                   y_ipc = base$ipc[bi + h] - base$ipc[bi - 1],
                   zshk  = shocks7[[k]], W7)[ok, , drop = FALSE]
  dd  <- dd[complete.cases(dd), ]
  ctl <- setdiff(names(dd), c("y_tcn", "y_ipc", "zshk"))
  ft  <- lm(reformulate(c("zshk", ctl), "y_tcn"), data = dd)
  fp  <- lm(reformulate(c("zshk", ctl), "y_ipc"), data = dd)
  iv  <- ivreg(as.formula(paste0("y_ipc ~ y_tcn + ", paste(ctl, collapse = " + "),
                                 " | zshk + ", paste(ctl, collapse = " + "))), data = dd)
  V1  <- tryCatch(NeweyWest(ft, lag = h + 1, prewhite = FALSE), error = function(e) NULL)
  F1  <- if (is.null(V1)) NA_real_ else unname(coef(ft)["zshk"])^2 / V1["zshk", "zshk"]
  tibble(h = h, shock = k, N = nrow(dd),
         b_tcn = unname(coef(ft)["zshk"]), se_tcn = nwse(ft, h, "zshk"),
         b_ipc = unname(coef(fp)["zshk"]), se_ipc = nwse(fp, h, "zshk"),
         erpt  = unname(coef(iv)["y_tcn"]), se_erpt = nwse(iv, h, "y_tcn"), F1 = F1)
}

lp6 <- map_dfr(shock_names, function(k) map_dfr(0:H, ~ lp6_hk(.x, k)))

cat("\nERPT condicional por shock (h=12) y F de primera etapa:\n")
print(lp6 %>% filter(h == 12) %>%
        transmute(shock = lab_shock[shock], erpt = round(erpt, 3),
                  ic = sprintf("[%.2f, %.2f]", erpt - z*se_erpt, erpt + z*se_erpt),
                  F1 = round(F1, 1), N))

cat("\nERPT del shock cambiario en el VAR7 (comparar con el punto 3):\n")
print(lp6 %>% filter(shock == "tcn", h %in% c(0, 3, 6, 12, 18, 24)) %>%
        transmute(h, erpt = round(erpt, 3), F1 = round(F1, 1)))

# --- Gráfico (escala fija; Fed y petróleo, F≈0, exceden el rango) ---

df6 <- lp6 %>% transmute(h, shock = factor(lab_shock[shock], levels = lab_shock),
                         est = erpt, lo = erpt - z*se_erpt, hi = erpt + z*se_erpt)
g6 <- ggplot(df6, aes(h, est)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), fill = "steelblue", alpha = 0.2) +
  geom_line(color = "steelblue4", linewidth = 0.9) +
  facet_wrap(~ shock, ncol = 3) + scale_x_continuous(breaks = seq(0, H, 6)) +
  coord_cartesian(ylim = c(-0.4, 0.55)) +
  labs(x = "Meses desde el shock", y = "ERPT",
       title = "ERPT condicional por tipo de shock (VAR de 7 variables)",
       subtitle = "Escala fija; las bandas de Fed y petróleo (F≈0) exceden el rango graficado") +
  tema_paper

print(g6)
guardar_graf(g6, "p6_erpt_condicional", width = 11, height = 6)

# --- Tabla para el informe ---
tabla_p6 <- lp6 %>%
  filter(h %in% c(6, 12, 24)) %>%
  mutate(txt = sprintf("%.3f [%.2f, %.2f]", erpt, erpt - z*se_erpt, erpt + z*se_erpt)) %>%
  dplyr::select(shock, h, txt) %>%
  pivot_wider(names_from = h, values_from = txt, names_prefix = "h") %>%
  left_join(lp6 %>% filter(h == 12) %>% transmute(shock, F12 = sprintf("%.1f", F1)), by = "shock") %>%
  mutate(shock = lab_shock[shock]) %>%
  dplyr::select(shock, h6, h12, h24, F12)
print(tabla_p6)

tabla_p6 %>%
  kbl(booktabs = TRUE, format = "latex", align = "lcccc",
      col.names = c("Shock", "$h=6$", "$h=12$", "$h=24$", "$F$ ($h{=}12$)"),
      caption = "ERPT condicional por tipo de shock estructural", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position", font_size = 9) %>%
  footnote(
    general = paste0(
      "Cada shock se infiere del VAR de 7 variables (orden Fed, petróleo, EBP, TI, EMBI, TCN, IPC) ",
      "y se usa como instrumento del cambio acumulado del TCN en la LP. Intervalos al ", NIVEL_BANDAS,
      "\\% (Newey-West, $h+1$ rezagos). El $F$ de primera etapa mide la relevancia del instrumento: ",
      "un ERPT es interpretable sólo si el shock mueve significativamente al TCN. Muestra: 2000m2--2024m7."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto6.tex"))

# ============================================================
# PUNTO 7
# ERPT por SIGNO, para cada uno de los 6 shocks del punto 6
#
# Se parte cada shock estructural en su parte positiva y negativa y se
# estima el ERPT de cada una por IV.
# ============================================================

lp7_hk <- function(h, k) {
  bi <- shocks7$bi
  ok <- (bi + h) <= nrow(base) & (bi - 1) >= 1
  zk <- shocks7[[k]]
  dd <- data.frame(y_tcn = base$tcn[bi + h] - base$tcn[bi - 1],
                   y_ipc = base$ipc[bi + h] - base$ipc[bi - 1],
                   zpos = zk * (zk > 0), zneg = zk * (zk < 0), W7)[ok, , drop = FALSE]
  dd  <- dd[complete.cases(dd), ]
  ctl <- setdiff(names(dd), c("y_tcn", "y_ipc", "zpos", "zneg"))
  side <- function(inst, other) {
    iv <- ivreg(as.formula(paste0("y_ipc ~ y_tcn + ", other, " + ", paste(ctl, collapse = " + "),
                                  " | ", inst, " + ", other, " + ", paste(ctl, collapse = " + "))), data = dd)
    f1 <- lm(reformulate(c(inst, other, ctl), "y_tcn"), data = dd)
    V1 <- tryCatch(NeweyWest(f1, lag = h + 1, prewhite = FALSE), error = function(e) NULL)
    F1 <- if (is.null(V1)) NA_real_ else unname(coef(f1)[inst])^2 / V1[inst, inst]
    c(est = unname(coef(iv)["y_tcn"]), se = nwse(iv, h, "y_tcn"), F1 = F1)
  }
  ep <- side("zpos", "zneg"); en <- side("zneg", "zpos")
  tibble(h = h, shock = k,
         erpt_pos = ep["est"], se_pos = ep["se"], F_pos = ep["F1"],
         erpt_neg = en["est"], se_neg = en["se"], F_neg = en["F1"])
}

lp7 <- map_dfr(shock_names, function(k) map_dfr(0:H, ~ lp7_hk(.x, k)))

cat("\nERPT por signo del shock (h=12):\n")
print(lp7 %>% filter(h == 12) %>%
        transmute(shock = lab_shock[shock],
                  erpt_pos = round(erpt_pos, 3), F_pos = round(F_pos, 1),
                  erpt_neg = round(erpt_neg, 3), F_neg = round(F_neg, 1)))

# --- Gráfico (escala fija; sólo el shock cambiario retiene señal) ---
df7 <- bind_rows(
  lp7 %>% transmute(h, shock = factor(lab_shock[shock], levels = lab_shock), signo = "Positivo",
                    est = erpt_pos, lo = erpt_pos - z*se_pos, hi = erpt_pos + z*se_pos),
  lp7 %>% transmute(h, shock = factor(lab_shock[shock], levels = lab_shock), signo = "Negativo",
                    est = erpt_neg, lo = erpt_neg - z*se_neg, hi = erpt_neg + z*se_neg))
g7 <- ggplot(df7, aes(h, est, color = signo, fill = signo)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.12, color = NA) +
  geom_line(linewidth = 0.85) +
  facet_wrap(~ shock, ncol = 3) + scale_x_continuous(breaks = seq(0, H, 6)) +
  coord_cartesian(ylim = c(-0.75, 1)) +
  scale_color_manual(values = c("Positivo" = "firebrick", "Negativo" = "steelblue4")) +
  scale_fill_manual(values  = c("Positivo" = "firebrick", "Negativo" = "steelblue4")) +
  labs(x = "Meses desde el shock", y = "ERPT", color = NULL, fill = NULL,
       title = "ERPT por signo del shock, para cada shock estructural",
       subtitle = "Escala fija; salvo el shock cambiario, la primera etapa débil (F<10) hace estallar las bandas") +
  tema_paper + theme(legend.position = "bottom")

print(g7)
guardar_graf(g7, "p7_erpt_signo", width = 11, height = 6)

# --- Tabla para el informe (ERPT por signo a h=12, con F de primera etapa) ---
tabla_p7 <- lp7 %>%
  filter(h == 12) %>%
  transmute(
    Shock = lab_shock[shock],
    `ERPT$^{+}$` = sprintf("%.3f", erpt_pos), `$F^{+}$` = sprintf("%.1f", F_pos),
    `ERPT$^{-}$` = sprintf("%.3f", erpt_neg), `$F^{-}$` = sprintf("%.1f", F_neg)
  )
print(tabla_p7)

tabla_p7 %>%
  kbl(booktabs = TRUE, format = "latex", align = "lcccc",
      col.names = c("Shock", "ERPT$^{+}$", "$F^{+}$", "ERPT$^{-}$", "$F^{-}$"),
      caption = "ERPT por signo del shock, para cada shock estructural ($h=12$)", escape = FALSE) %>%
  kable_styling(latex_options = "hold_position", font_size = 9) %>%
  footnote(
    general = paste0(
      "Cada shock del punto 6 se parte en pieza positiva y negativa; el ERPT de cada signo se estima ",
      "por IV con esa pieza como instrumento y la otra como control. Los $F$ de primera etapa, casi ",
      "todos inferiores a 10, muestran que al condicionar por shock y además partir por signo se agota ",
      "la variación identificante: salvo el shock cambiario propio, los ERPT por signo no están identificados."
    ),
    escape = FALSE, threeparttable = TRUE
  ) %>%
  cat(file = file.path(dir_salida, "tabla_punto7.tex"))

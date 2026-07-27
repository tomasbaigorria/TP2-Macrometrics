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

# ========================================
#  PROJET : GRACE TWS + EGMS (Caen)
#  CRS : Détection automatique
#  Sauvegarde GRACE CSV + EGMS CSV
#  Corrélation mensuelle robuste (2019–2023)
# ========================================

# --- Packages ---
library(ncdf4)
library(CFtime)
library(RColorBrewer)
library(fields)     # image.plot
library(sf)
library(leaflet)
library(lubridate)

# -----------------------------
# 1) Fonctions utilitaires
# -----------------------------

# Détection simple du CRS EGMS (heuristique)
detect_crs_egms <- function(df) {
  
  if (!all(c("easting", "northing") %in% names(df))) {
    stop("Colonnes 'easting' et 'northing' manquantes.")
  }
  
  e <- df$easting
  n <- df$northing
  
  max_e <- max(e, na.rm = TRUE)
  max_n <- max(n, na.rm = TRUE)
  min_e <- min(e, na.rm = TRUE)
  min_n <- min(n, na.rm = TRUE)
  
  # --- Détection WGS84 (EPSG:4326) ---
  if (max_e <= 180 && max_e >= -180 && max_n <= 90 && max_n >= -90) {
    cat("CRS détecté : WGS84 (EPSG:4326)")
    return(4326)
  }
  
  # --- Détection LAEA Europe (EPSG:3035) ---
  # Domaine EGMS 100km typique : 2e6 à 7e6
  if (max_e > 2e6 && max_e < 7e6 && max_n > 1e6 && max_n < 7e6) {
    cat(" CRS détecté : ETRS89 / LAEA Europe (EPSG:3035)")
    return(3035)
  }
  
  # --- Détection UTM ---
  # Eastings UTM : 166 000 – 833 000
  # Northings : 0 – 10 000 000
  if (max_e > 160000 && max_e < 900000 && max_n > 0 && max_n < 10000000) {
    # Estimation de la zone UTM basée sur l'easting
    # Zone 31N typique pour la France : easting ~200k-700k
    mean_e <- mean(e, na.rm = TRUE)
    
    # France = Zone 30-32
    if (mean_e > 150000 && mean_e < 400000) {
      zone <- 30
    } else if (mean_e >= 400000 && mean_e < 700000) {
      zone <- 31
    } else {
      zone <- 32
    }
    
    epsg <- 32600 + zone
    cat(paste0("CRS détecté : UTM Zone ", zone, epsg, ")"))
    return(epsg)
  }
  
  # --- Si aucune détection ---
  cat("CRS non reconnu automatiquement. Valeurs easting/northing inhabituelles.")
  cat("  max_e =", max_e, "| max_n =", max_n, "\n")
  return(NA)
}

# Récupérer le TWS au pixel le plus proche (sans trier lon/lat)
extract_grace_pixel <- function(lon, lat, tws_array, target_lon, target_lat) {
  lon_idx <- which.min(abs(lon - target_lon))
  lat_idx <- which.min(abs(lat - target_lat))
  list(
    lon_idx = lon_idx,
    lat_idx = lat_idx,
    lon_val = lon[lon_idx],
    lat_val = lat[lat_idx],
    ts      = tws_array[lon_idx, lat_idx, ]
  )
}

# -----------------------------
# 2) GRACE : lecture NetCDF
# -----------------------------

ncpath <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/G3P_v1.12/TWS/"
ncname <- "tws_g3p_200204_202309_monthly_05_degree_costg_gravis.nc"
ncfname <- file.path(ncpath, ncname)
if (!file.exists(ncfname)) stop("Fichier GRACE introuvable : ", ncfname)

ncin <- nc_open(ncfname)
cat("Fichier GRACE ouvert.\n")

lon       <- ncvar_get(ncin, "lon")
lat       <- ncvar_get(ncin, "lat")
time      <- ncvar_get(ncin, "time")
tunits    <- ncatt_get(ncin, "time", "units")
tws_array <- ncvar_get(ncin, "tws")

# Remplacer FillValue
fillvalue <- ncatt_get(ncin, "tws", "_FillValue")
if (!fillvalue$hasatt) fillvalue <- ncatt_get(ncin, "tws", "missing_value")
if (fillvalue$hasatt) {
  tws_array[tws_array == fillvalue$value] <- NA
  cat("Valeurs de remplissage remplacées.\n")
}

# Décodage du temps
calendar_att <- ncatt_get(ncin, "time", "calendar")
calendar <- if (calendar_att$hasatt) calendar_att$value else "proleptic_gregorian"

timestamps <- NULL
tryCatch({
  cf <- CFtime(tunits$value, calendar = calendar, time)
  timestamps <- format(cf)
}, error = function(e) {
  origin <- as.Date(sub(".* since ", "", tunits$value))
  timestamps <- format(origin + time, "%Y-%m-%d")
})

time_cf <- as.POSIXct(timestamps, tz = "UTC")
time_components <- data.frame(
  year  = as.numeric(format(time_cf, "%Y")),
  month = as.numeric(format(time_cf, "%m")),
  day   = as.numeric(format(time_cf, "%d"))
)

# -----------------------------
# 3) Caen : extraction pixel GRACE
# -----------------------------

target_lon <- -0.37
target_lat <- 49.18

px <- extract_grace_pixel(lon, lat, tws_array, target_lon, target_lat)
cat("Point GRACE le plus proche de Caen :", px$lon_val, px$lat_val, "\n")

tws_timeseries <- px$ts

# Sauvegarde GRACE Caen
grace_csv_file <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/grace_caen.csv"
grace_caen_df <- data.frame(
  longitude = px$lon_val,
  latitude  = px$lat_val,
  date      = format(time_cf, "%Y-%m-%d"),
  tws_mm    = tws_timeseries
)
write.csv(grace_caen_df, grace_csv_file, row.names = FALSE)
cat("Fichier GRACE Caen sauvegardé :", grace_csv_file, "\n")

# Plot TWS global
plot(time_cf, tws_timeseries, type = "l", col = "blue", lwd = 2,
     main = paste("TWS à Caen (", round(px$lon_val, 3), ",", round(px$lat_val, 3), ")"),
     xlab = "Temps", ylab = "TWS (mm)")
points(time_cf, tws_timeseries, pch = 20, cex = 0.5, col = "blue")

# -----------------------------
# 4) EGMS : lecture + reprojection
# -----------------------------

egms_csv_path <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/EGMS_L3_E35N29_100km_U_2019_2023_1/EGMS_L3_E35N29_100km_U_2019_2023_1.csv"
if (!file.exists(egms_csv_path)) {
  nc_close(ncin)
  stop("Fichier EGMS introuvable : ", egms_csv_path)
}

egms_data <- read.csv(egms_csv_path, header = TRUE)
cat("EGMS chargé :", nrow(egms_data), "points\n")

crs_code <- detect_crs_egms(egms_data)

sf_egms  <- st_as_sf(egms_data, coords = c("easting", "northing"), crs = crs_code, remove = FALSE)
sf_wgs84 <- st_transform(sf_egms, 4326)
coords_e <- st_coordinates(sf_wgs84)

egms_data$longitude <- coords_e[, 1]
egms_data$latitude  <- coords_e[, 2]

# Sauvegarde EGMS reprojeté (chemin Mac)
egms_csv_out <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/egms_caen_wgs84.csv"
write.csv(egms_data, egms_csv_out, row.names = FALSE)
cat("EGMS sauvegardé (WGS84) :", egms_csv_out, "\n")

# -----------------------------
# 5) Filtrage EGMS autour de Caen
# -----------------------------

egms_data$dist_deg <- sqrt((egms_data$longitude - target_lon)^2 + (egms_data$latitude - target_lat)^2)

# Rayon progressif : 0.5° -> 1° -> 2°
egms_caen <- egms_data[egms_data$dist_deg <= 0.5, ]
if (nrow(egms_caen) == 0) egms_caen <- egms_data[egms_data$dist_deg <= 1.0, ]
if (nrow(egms_caen) == 0) egms_caen <- egms_data[egms_data$dist_deg <= 2.0, ]

if (nrow(egms_caen) == 0) {
  nc_close(ncin)
  stop("Aucun point EGMS trouvé près de Caen (même avec 2°).")
}

cat("Points EGMS retenus autour de Caen :", nrow(egms_caen), "\n")

# Statistiques
if (!all(c("mean_velocity", "rmse") %in% names(egms_caen))) {
  nc_close(ncin)
  stop("EGMS : colonnes 'mean_velocity' et/ou 'rmse' absentes.")
}

mean_vel <- mean(egms_caen$mean_velocity, na.rm = TRUE)
rmse_val <- mean(egms_caen$rmse, na.rm = TRUE)
cat("EGMS - Vitesse moyenne :", round(mean_vel, 3), "mm/an | RMSE :", round(rmse_val, 3), "\n")

# -----------------------------
# 6) Corrélation mensuelle robuste (2019–2023)
# -----------------------------

# GRACE 2019-2023
idx_grace <- which(time_components$year >= 2019 & time_components$year <= 2023)
tws_19_23  <- tws_timeseries[idx_grace]
date_grace <- time_cf[idx_grace]
key_grace  <- format(date_grace, "%Y-%m")

# EGMS dates
date_cols <- grep("^X\\d{8}$", names(egms_caen), value = TRUE)
if (length(date_cols) == 0) {
  nc_close(ncin)
  stop("Aucune colonne EGMS de type XYYYYMMDD trouvée.")
}

egms_dates_raw <- as.Date(sub("^X", "", date_cols), format = "%Y%m%d")
egms_month     <- floor_date(egms_dates_raw, "month")
key_egms       <- format(egms_month, "%Y-%m")

egms_matrix <- as.matrix(egms_caen[, date_cols, drop = FALSE])
storage.mode(egms_matrix) <- "numeric"

# Moyenne mensuelle sur TOUS les points
n_points <- nrow(egms_matrix)
egms_monthly_mean <- tapply(
  X     = as.vector(egms_matrix),                # empile par colonnes
  INDEX = rep(key_egms, each = n_points),        # chaque date répétée n_points fois
  FUN   = mean,
  na.rm = TRUE
)

df_egms_monthly <- data.frame(
  year_month = names(egms_monthly_mean),
  egms_mm    = as.numeric(egms_monthly_mean),
  stringsAsFactors = FALSE
)

df_grace <- data.frame(
  year_month = key_grace,
  date       = as.Date(format(date_grace, "%Y-%m-01")),
  tws_mm     = tws_19_23,
  stringsAsFactors = FALSE
)

df_merged <- merge(df_grace, df_egms_monthly, by = "year_month", all.x = TRUE)
df_merged <- df_merged[order(df_merged$date), ]

valid <- complete.cases(df_merged$tws_mm, df_merged$egms_mm)
df_valid <- df_merged[valid, ]

cat("\n=== RÉSULTATS AGRÉGATION MENSUELLE ===\n")
cat("Mois totaux GRACE (2019–2023) :", nrow(df_grace), "\n")
cat("Mois valides communs :", nrow(df_valid), "\n")

if (nrow(df_valid) < 10) {
  warning("Pas assez de mois communs (>=10 conseillé) pour une corrélation robuste.")
} else {
  correl <- cor(df_valid$tws_mm, df_valid$egms_mm, use = "complete.obs")
  cat("Corrélation mensuelle TWS–EGMS (2019–2023) :", round(correl, 4), "\n")
  
  # Normalisation pour visu
  tws_norm  <- as.numeric(scale(df_valid$tws_mm))
  egms_norm <- as.numeric(scale(df_valid$egms_mm))
  
  plot(df_valid$date, tws_norm, type = "l", col = "blue", lwd = 2,
       ylab = "Valeur normalisée", xlab = "Année",
       main = paste0("TWS GRACE vs EGMS Caen (mensuel 2019–2023) | r = ", round(correl, 3)))
  lines(df_valid$date, egms_norm, col = "red", lwd = 2)
  abline(h = 0, lty = 3, col = "gray")
  legend("topright", legend = c("TWS (GRACE)", "EGMS (moyenne mensuelle)"),
         col = c("blue", "red"), lwd = 2, cex = 0.9)

}


# Sauvegarde du tableau mensuel commun
out_monthly <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/caen_tws_egms_mensuel_2019_2023.csv"
write.csv(df_valid, out_monthly, row.names = FALSE)
cat("Tableau mensuel commun sauvegardé :", out_monthly, "\n")

# -----------------------------
# --- Carte TWS autour de Caen (mois m) ---
m <- 1
if (m >= 1 && m <= length(time)) {
  tws_slice <- tws_array[, , m]
  
  lon_range <- c(target_lon - 5, target_lon + 5)
  lat_range <- c(target_lat - 5, target_lat + 5)
  
  lon_sub <- which(lon >= lon_range[1] & lon <= lon_range[2])
  lat_sub <- which(lat >= lat_range[1] & lat <= lat_range[2])
  
  if (length(lon_sub) > 1 && length(lat_sub) > 1) {
    
    # extraire
    lon_r <- lon[lon_sub]
    lat_r <- lat[lat_sub]
    tws_r <- tws_slice[lon_sub, lat_sub]
    
    # trier lon/lat + réordonner la matrice pour que image() accepte
    o_lon <- order(lon_r)
    o_lat <- order(lat_r)
    
    lon_r2 <- lon_r[o_lon]
    lat_r2 <- lat_r[o_lat]
    tws_r2 <- tws_r[o_lon, o_lat]
    
    image(lon_r2, lat_r2, tws_r2,
          col = rev(brewer.pal(10, "RdBu")),
          main = paste("TWS près de Caen -", format(time_cf[m], "%Y-%m-%d")),
          xlab = "Longitude", ylab = "Latitude")
    
    points(target_lon, target_lat, pch = 19, col = "black", cex = 1.5)
    text(target_lon, target_lat, "Caen", pos = 3, col = "black", cex = 1.2)
    
    image.plot(legend.only = TRUE,
               col = rev(brewer.pal(10, "RdBu")),
               zlim = range(tws_r2, na.rm = TRUE))
    
  } else {
    warning("Sous-grille trop petite pour image() (pas assez de points lon/lat).")
  }
} else {
  warning("m est hors bornes.")
}

# -----------------------------
# 8) Carte interactive Leaflet
# -----------------------------

leaflet(data = egms_caen) %>%
  addTiles() %>%
  addCircleMarkers(
    lng = ~longitude, lat = ~latitude,
    radius = 3, color = "red", fillOpacity = 0.8,
    popup = ~paste0("Vélocité : ", round(mean_velocity, 2), " mm/an<br>",
                    "RMSE : ", round(rmse, 2))
  ) %>%
  addMarkers(lng = target_lon, lat = target_lat,
             popup = "<b>Caen</b><br>Pixel GRACE")


# --- Nettoyage ---
nc_close(ncin)


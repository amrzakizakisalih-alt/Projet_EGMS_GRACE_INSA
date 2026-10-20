# ========================================
#  PROJET : GRACE TWS + EGMS (Caen)
#  Pixel GRACE (0.5° × 0.5°) et Pixel EGMS (1° x 1°)
#  CRS : Détection automatique
#  Sauvegarde GRACE CSV + EGMS CSV
#  Corrélation mensuelle (2019–2023)
# ========================================

# --- Packages ---
library(ncdf4)
library(CFtime)
library(RColorBrewer)
library(fields)
library(sf)
library(leaflet)
library(lubridate)

# =============================================
# 1) FONCTIONS
# =============================================

detecter_crs_egms <- function(df) {
  
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

# Extraction du pixel GRACE le plus proche
extraire_pixel_grace <- function(lon, lat, tws_array, lon_cible, lat_cible) {
  idx_lon <- which.min(abs(lon - lon_cible))
  idx_lat <- which.min(abs(lat - lat_cible))
  list(
    idx_lon = idx_lon,
    idx_lat = idx_lat,
    lon_val = lon[idx_lon],
    lat_val = lat[idx_lat],
    serie_tws = tws_array[idx_lon, idx_lat, ]
  )
}

# =============================================
# 2) LECTURE GRACE (NetCDF)
# =============================================

cat("\n========== ÉTAPE 1 : CHARGEMENT GRACE ==========\n")

chemin_nc  <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/G3P_v1.12/TWS/"
nom_nc     <- "tws_g3p_200204_202309_monthly_05_degree_costg_gravis.nc"
fichier_nc <- file.path(chemin_nc, nom_nc)

if (!file.exists(fichier_nc)) {
  stop("Fichier GRACE introuvable : ", fichier_nc)
}

nc_in <- nc_open(fichier_nc)
cat("Fichier GRACE ouvert\n")

# Variables
lon <- ncvar_get(nc_in, "lon")
lat <- ncvar_get(nc_in, "lat")
temps <- ncvar_get(nc_in, "time")
unite_temps <- ncatt_get(nc_in, "time", "units")
tws_array <- ncvar_get(nc_in, "tws")

cat("Dimensions : Lon =", length(lon), "| Lat =", length(lat), "| Temps =", length(temps), "\n")

# Valeurs manquantes
valeur_manquante <- ncatt_get(nc_in, "tws", "_FillValue")
if (!valeur_manquante$hasatt) {
  valeur_manquante <- ncatt_get(nc_in, "tws", "missing_value")
}
if (valeur_manquante$hasatt) {
  tws_array[tws_array == valeur_manquante$value] <- NA
  cat("Valeurs manquantes remplacées\n")
}

# Décodage des dates
calendrier_att <- ncatt_get(nc_in, "time", "calendar")
calendrier <- if (calendrier_att$hasatt) calendrier_att$value else "proleptic_gregorian"

horodatages <- NULL
tryCatch({
  cf <- CFtime(unite_temps$value, calendar = calendrier, temps)
  horodatages <- format(cf)
}, error = function(e) {
  origine <- as.Date(sub(".* since ", "", unite_temps$value))
  horodatages <- format(origine + temps, "%Y-%m-%d")
})

temps_cf <- as.POSIXct(horodatages, tz = "UTC")
composantes_temps <- data.frame(
  annee = as.numeric(format(temps_cf, "%Y")),
  mois  = as.numeric(format(temps_cf, "%m")),
  jour  = as.numeric(format(temps_cf, "%d"))
)

cat("Dates : ", format(range(temps_cf), "%Y-%m-%d"), "\n")

# =============================================
# 3) EXTRACTION PIXEL GRACE CAEN
# =============================================

cat("\n========== ÉTAPE 2 : PIXEL GRACE CAEN ==========\n")

lon_caen <- -0.37
lat_caen <- 49.18

pixel <- extraire_pixel_grace(lon, lat, tws_array, lon_caen, lat_caen)

cat("Pixel GRACE :\n")
cat("Longitude : ", pixel$lon_val, "°\n")
cat("Latitude  : ", pixel$lat_val, "°\n")

serie_tws <- pixel$serie_tws

# Sauvegarde GRACE
fichier_grace_csv <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/grace_caen.csv"
df_grace_caen <- data.frame(
  longitude = pixel$lon_val,
  latitude = pixel$lat_val,
  date = format(temps_cf, "%Y-%m-%d"),
  tws_mm = serie_tws
)
write.csv(df_grace_caen, fichier_grace_csv, row.names = FALSE)
cat("GRACE sauvegardé :", fichier_grace_csv, "\n")

# =============================================
# 4) LECTURE ET REPROJECTION EGMS
# =============================================

cat("\n========== ÉTAPE 3 : CHARGEMENT EGMS ==========\n")

chemin_egms_csv <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/EGMS_L3_E35N29_100km_U_2019_2023_1/EGMS_L3_E35N29_100km_U_2019_2023_1.csv"

if (!file.exists(chemin_egms_csv)) {
  nc_close(nc_in)
  stop("Fichier EGMS introuvable : ", chemin_egms_csv)
}

donnees_egms <- read.csv(chemin_egms_csv, header = TRUE)
cat("EGMS chargé :", nrow(donnees_egms), "points\n")

# Détection et reprojection CRS
code_crs <- detecter_crs_egms(donnees_egms)
sf_egms <- st_as_sf(donnees_egms, coords = c("easting", "northing"), crs = code_crs, remove = FALSE)
sf_wgs84 <- st_transform(sf_egms, 4326)
coords_wgs <- st_coordinates(sf_wgs84)
donnees_egms$longitude <- coords_wgs[, 1]
donnees_egms$latitude  <- coords_wgs[, 2]

cat("Reprojection WGS84 effectuée\n")

# Sauvegarde EGMS reprojeté
fichier_egms_out <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/egms_caen_wgs84.csv"
write.csv(donnees_egms, fichier_egms_out, row.names = FALSE)
cat("EGMS sauvegardé :", fichier_egms_out, "\n")

# =============================================
# 5) FILTRAGE PAR PIXEL GRACE EXACT
# =============================================

cat("\n========== ÉTAPE 4 : FILTRAGE PIXEL GRACE ==========\n")

# Résolution GRACE : 0.5°
resolution_grace <- 0.5
demi_resolution <- resolution_grace / 2

# Centre du pixel
lon_centre_grace <- pixel$lon_val
lat_centre_grace <- pixel$lat_val

# Limites du pixel
bbox_grace <- list(lon_min = lon_centre_grace - demi_resolution, lon_max = lon_centre_grace + demi_resolution,lat_min = lat_centre_grace - demi_resolution,lat_max = lat_centre_grace + demi_resolution)

cat("Limites pixel GRACE (0.5° × 0.5°) :\n")
cat("Longitude : [", round(bbox_grace$lon_min, 4), ", ", round(bbox_grace$lon_max, 4), "]\n")
cat("Latitude  : [", round(bbox_grace$lat_min, 4), ", ", round(bbox_grace$lat_max, 4), "]\n")

# Filtrage des points EGMS dans le pixel
egms_dans_pixel <- donnees_egms[donnees_egms$longitude >= bbox_grace$lon_min & donnees_egms$longitude <= bbox_grace$lon_max & donnees_egms$latitude  >= bbox_grace$lat_min & donnees_egms$latitude  <= bbox_grace$lat_max,]

cat("Points EGMS dans le pixel : ", nrow(egms_dans_pixel), "\n")

if (nrow(egms_dans_pixel) == 0) {
  cat("Aucun point, extension à 1° × 1°...\n") 
  egms_dans_pixel <- donnees_egms[donnees_egms$longitude >= (lon_centre_grace - 0.5) & donnees_egms$longitude <= (lon_centre_grace + 0.5) & donnees_egms$latitude  >= (lat_centre_grace - 0.5) & donnees_egms$latitude  <= (lat_centre_grace + 0.5),]
  cat("Points dans pixel élargi : ", nrow(egms_dans_pixel), "\n")
}

if (nrow(egms_dans_pixel) == 0) {
  nc_close(nc_in)
  stop("Aucun point EGMS trouvé.")
}

# =============================================
# 6) CORRÉLATION MENSUELLE (2019–2023)
# =============================================

cat("\n========== ÉTAPE 5 : CORRÉLATION MENSUELLE ==========\n")

# Filtrer GRACE 2019-2023
idx_grace_19_23 <- which(composantes_temps$annee >= 2019 & composantes_temps$annee <= 2023)
tws_19_23  <- serie_tws[idx_grace_19_23]
date_grace_19_23 <- temps_cf[idx_grace_19_23]
cle_grace  <- format(date_grace_19_23, "%Y-%m")

cat("Mois GRACE (2019-2023) : ", length(tws_19_23), "\n")

# Dates EGMS
colonnes_dates <- grep("^X\\d{8}$", names(egms_dans_pixel), value = TRUE)

if (length(colonnes_dates) == 0) {
  nc_close(nc_in)
  stop("Aucune colonne EGMS XYYYYMMDD trouvée.")
}

dates_egms_brut <- as.Date(sub("^X", "", colonnes_dates), format = "%Y%m%d")
mois_egms <- floor_date(dates_egms_brut, "month")
cle_egms  <- format(mois_egms, "%Y-%m")

# Matrice EGMS
matrice_egms <- as.matrix(egms_dans_pixel[, colonnes_dates, drop = FALSE])
storage.mode(matrice_egms) <- "numeric"

# Moyenne mensuelle
n_points <- nrow(matrice_egms)
moyenne_mensuelle_egms <- tapply(X = as.vector(matrice_egms), INDEX = rep(cle_egms, each = n_points),FUN = mean, na.rm = TRUE)

df_egms_mensuel <- data.frame(annee_mois = names(moyenne_mensuelle_egms), egms_mm    = as.numeric(moyenne_mensuelle_egms),stringsAsFactors = FALSE)
df_grace_mensuel <- data.frame( annee_mois = cle_grace,date = as.Date(format(date_grace_19_23, "%Y-%m-01")),tws_mm = tws_19_23, stringsAsFactors = FALSE)

# Fusion
df_fusionne <- merge(df_grace_mensuel, df_egms_mensuel, by = "annee_mois", all.x = TRUE)
df_fusionne <- df_fusionne[order(df_fusionne$date), ]

valides <- complete.cases(df_fusionne$tws_mm, df_fusionne$egms_mm)
df_valide <- df_fusionne[valides, ]

cat("Mois communs valides : ", nrow(df_valide), "\n")

if (nrow(df_valide) >= 10) {
  correlation <- cor(df_valide$tws_mm, df_valide$egms_mm, use = "complete.obs")
  cat("\n CORRÉLATION TWS–EGMS : ", round(correlation, 4), "\n\n")
  
  # Normalisation
  tws_norm  <- as.numeric(scale(df_valide$tws_mm))
  egms_norm <- as.numeric(scale(df_valide$egms_mm))
  
  # Graphique comparatif
 
  
  plot(df_valide$date, tws_norm, type = "l", col = "dodgerblue", lwd = 2.5, ylab = "Valeur normalisée", xlab = "", main = paste0("TWS GRACE vs EGMS à Caen (2019-2023) | r = ", round(correlation, 3)),las = 1)
  lines(df_valide$date, egms_norm, col = "red3", lwd = 2.5)
  abline(h = 0, lty = 3, col = "gray40")
  grid()
  legend("topright", legend = c("TWS (GRACE)", "EGMS (moyenne)"), col = c("dodgerblue", "red3"), lwd = 2.5, cex = 0.9)
  
} else {
  warning("Moins de 10 mois communs.")
}

# Sauvegarde
fichier_mensuel <- "/home/mouad/Bureau/ProjetS7MSRO_SHERADJDRISSIMouad_ZAKISALIHAmr_SECKCheikh/data/caen_tws_egms_mensuel_2019_2023.csv"
write.csv(df_valide, fichier_mensuel, row.names = FALSE)
cat("Tableau mensuel sauvegardé :", fichier_mensuel, "\n")



# =============================================
# 7) CARTE INTERACTIVE LEAFLET
# =============================================
cat("\n========== ÉTAPE 6 : CARTE INTERACTIVE ==========\n")

carte_interactive <- leaflet() %>%
  addTiles(group = "OpenStreetMap") %>%
  addProviderTiles("Esri.WorldImagery", group = "Satellite") %>%
  
  # Rectangle du pixel GRACE
  addRectangles(
    lng1 = bbox_grace$lon_min, lat1 = bbox_grace$lat_min,
    lng2 = bbox_grace$lon_max, lat2 = bbox_grace$lat_max,
    color = "blue", fillColor = "lightblue", fillOpacity = 0.25,
    weight = 3, 
    popup = paste0("<b>Pixel GRACE</b><br>",
                   "Résolution : 0.5° × 0.5° (~50 km)<br>",
                   "Centre : (", round(lon_centre_grace, 4), ", ", 
                   round(lat_centre_grace, 4), ")")
  ) %>%
  
  # Ville de Caen
  addMarkers(
    lng = lon_caen, 
    lat = lat_caen,
    popup = paste0("<b>CAEN</b><br>",
                   "Coordonnées : (", lon_caen, ", ", lat_caen, ")")
  ) %>%
  
  # Centre du pixel GRACE
  addMarkers(
    lng = lon_centre_grace, 
    lat = lat_centre_grace,
    popup = paste0("<b>Centre pixel GRACE</b><br>",
                   "TWS moyen : ", round(mean(serie_tws, na.rm = TRUE), 2), " mm")
  ) %>%
  
  # Points EGMS
  addCircleMarkers(
    data = egms_dans_pixel,
    lng = ~longitude, 
    lat = ~latitude,
    radius = 5, 
    color = "red",
    fillOpacity = 0.8,
    stroke = TRUE,
    weight = 1,
    popup = ~paste0("<b>Point EGMS</b><br>",
                    "Position : (", round(longitude, 4), ", ", 
                    round(latitude, 4), ")")
  ) %>%
  
  # Contrôles
  addLayersControl(
    baseGroups = c("OpenStreetMap", "Satellite"),
    options = layersControlOptions(collapsed = FALSE)
  )

print(carte_interactive)
cat("Carte interactive générée\n")
# =============================================
# 8) NETTOYAGE
# =============================================

nc_close(nc_in)
cat("\n==========================================================\n")
cat("\n RÉSUMÉ :\n")
cat("Pixel GRACE   : (", round(pixel$lon_val, 4), ", ", round(pixel$lat_val, 4), ")\n")
cat("Points EGMS   : ", nrow(egms_dans_pixel), "\n")
cat("Corrélation   : r = ", round(correlation, 4), "\n")

cat("\n Fichiers :\n")
cat(fichier_grace_csv, "\n")
cat(fichier_egms_out, "\n")
cat(fichier_mensuel, "\n")
cat("\n==========================================================\n")
cat("\nANALYSE TERMINÉE AVEC SUCCÈS !\n")


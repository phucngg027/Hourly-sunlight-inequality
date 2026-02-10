suppressPackageStartupMessages({ 
  library(sf)
  library(terra)
  library(suncalc)
  library(osmdata)
})

sf::sf_use_s2(FALSE)
sf::sf_proj_network(TRUE)

dir0 <- "your path"
bld_path   <- file.path(dir0, "your path")
roads_path <- file.path(dir0, "your path")

LAT   <- 41.3166411809927
LON   <- -72.92433810053593
date0 <- as.Date("2025-06-21")     
tz0   <- "America/New_York"

zoom_km   <- 0.60
res_m     <- 2.0
h_default <- 12

hours_vec <- c(8, 9, 10, 11, 13, 15, 17)
ncol_pan  <- 4

# Roads
road_col_base   <- "#9AA6B2"
road_alpha_base <- 0.70
road_blur_cells <- 3

# Buildings (outline only, crisp)
bld_edge <- "#768391"
bld_lwd  <- 1.10

# Parks (STRONGER)
add_osm_green <- TRUE
park_fill     <- "#BFE3CF"
park_alpha    <- 1.00
park_border   <- NA
park_edge_soft_m <- 6     

# optional hatch
park_hatch_on <- TRUE
park_hatch_col   <- "#A9C8B6"
park_hatch_alpha <- 0.22
park_hatch_lwd   <- 1.0

# Shadow (crisper)
shadow_q       <- 0.94
shadow_gamma   <- 0.65
sigma_m_soft   <- 3.5
shadow_alpha   <- 0.78

# Buffers
bld_buf_cells <- 7                
park_mask_thresh <- 0.03          

# Output
out_dir <- file.path(dir0, sprintf("_hourly_panels_zoom_%dm", as.integer(round(zoom_km*1000))))
if (!dir.exists(out_dir)) dir.create(out_dir, recursive = TRUE)

out_png_grid <- file.path(out_dir, sprintf("Yale_GRID_SUMMER_PARKS_NYT_%s_zoom_%dm.png",
                                           date0, as.integer(round(zoom_km*1000))))
px_panel <- 1500
dpi0 <- 300

# CRS
crs_bld   <- "EPSG:6434"
crs_roads <- "EPSG:4269"
crs_utm   <- "EPSG:32618"

odd_int <- function(x, default = 3L) {
  x <- suppressWarnings(as.integer(x))
  if (length(x) != 1 || is.na(x) || x < 1) x <- default
  if (x %% 2 == 0) x <- x + 1
  x
}

gauss_kernel <- function(size, sigma) {
  size <- odd_int(size, 31L)
  ax <- seq(-(size-1)/2, (size-1)/2)
  K <- outer(ax, ax, function(x,y) exp(-(x^2 + y^2)/(2*sigma^2)))
  K / sum(K)
}

clamp01_rast <- function(x) {
  x <- terra::ifel(x < 0, 0, x)
  x <- terra::ifel(x > 1, 1, x)
  x
}

utm_bbox_from_lonlat <- function(lon, lat, zoom_km, crs_out) {
  pt_ll  <- sf::st_sfc(sf::st_point(c(lon, lat)), crs = 4326)
  pt_utm <- suppressWarnings(sf::st_transform(pt_ll, crs_out))
  xy <- sf::st_coordinates(pt_utm)[1, ]
  if (!all(is.finite(xy))) stop("UTM transform returned NA/Inf. Check PROJ network.")
  half <- zoom_km * 1000
  list(
    xmin = as.numeric(xy[1] - half),
    xmax = as.numeric(xy[1] + half),
    ymin = as.numeric(xy[2] - half),
    ymax = as.numeric(xy[2] + half),
    center = xy
  )
}

bbox_poly_sf <- function(bbox_list, crs) {
  bb <- sf::st_bbox(c(
    xmin = bbox_list$xmin,
    ymin = bbox_list$ymin,
    xmax = bbox_list$xmax,
    ymax = bbox_list$ymax
  ), crs = sf::st_crs(crs))
  sf::st_as_sfc(bb)
}

sf_clean_polygons <- function(x) {
  x <- sf::st_make_valid(x)
  x <- x[!sf::st_is_empty(x), ]
  if (nrow(x) == 0) return(x)
  gt <- unique(as.character(sf::st_geometry_type(x)))
  if (any(gt %in% "GEOMETRYCOLLECTION")) {
    x <- sf::st_collection_extract(x, "POLYGON", warn = FALSE)
  }
  gt <- unique(as.character(sf::st_geometry_type(x)))
  if (any(gt %in% "MULTIPOLYGON")) {
    x <- suppressWarnings(sf::st_cast(x, "POLYGON"))
  }
  x <- x[!sf::st_is_empty(x), ]
  x
}

alpha_ramp <- function(hex_col, alpha_max = 1.0, n = 256L, pow = 0.35) {
  x <- seq(0, 1, length.out = n)
  x <- x^pow
  vapply(x, function(a) grDevices::adjustcolor(hex_col, alpha.f = alpha_max * a), character(1))
}

get_osm_green_cached <- function(lon, lat, zoom_km, crs_out, cache_path, timeout = 240) {
  if (file.exists(cache_path)) {
    x <- suppressWarnings(sf::st_read(cache_path, quiet = TRUE))
    if (!inherits(x, "try-error") && nrow(x) > 0) return(sf::st_transform(x, crs_out))
  }
  
  half_lat <- zoom_km / 111.0
  half_lon <- zoom_km / (111.0 * cos(lat * pi/180))
  bb <- c(lon - half_lon, lat - half_lat, lon + half_lon, lat + half_lat)
  
  fetch <- function(key, value) {
    q <- osmdata::opq(bbox = bb, timeout = timeout) |>
      osmdata::add_osm_feature(key = key, value = value)
    tryCatch(osmdata::osmdata_sf(q), error = function(e) NULL)
  }
  
  union_geom <- function(res) {
    if (is.null(res)) return(NULL)
    geoms <- list()
    if (!is.null(res$osm_polygons) && nrow(res$osm_polygons) > 0) {
      geoms[[length(geoms)+1]] <- sf::st_geometry(sf::st_as_sf(res$osm_polygons))
    }
    if (!is.null(res$osm_multipolygons) && nrow(res$osm_multipolygons) > 0) {
      geoms[[length(geoms)+1]] <- sf::st_geometry(sf::st_as_sf(res$osm_multipolygons))
    }
    if (length(geoms) == 0) return(NULL)
    g <- do.call(c, geoms)
    g <- g[!sf::st_is_empty(g)]
    if (length(g) == 0) return(NULL)
    g <- suppressWarnings(sf::st_make_valid(g))
    g <- g[!sf::st_is_empty(g)]
    if (length(g) == 0) return(NULL)
    suppressWarnings(sf::st_union(g))
  }
  
  inc1 <- fetch("leisure", c("park","garden","common"))
  inc2 <- fetch("landuse", c("grass","recreation_ground","village_green","park"))
  inc3 <- fetch("natural", c("wood","scrub","grassland"))
  
  g1 <- union_geom(inc1); g2 <- union_geom(inc2); g3 <- union_geom(inc3)
  inc_all <- sf::st_sfc(crs = 4326)
  if (!is.null(g1) && length(g1) > 0) inc_all <- c(inc_all, g1)
  if (!is.null(g2) && length(g2) > 0) inc_all <- c(inc_all, g2)
  if (!is.null(g3) && length(g3) > 0) inc_all <- c(inc_all, g3)
  inc_all <- inc_all[!sf::st_is_empty(inc_all)]
  if (length(inc_all) == 0) return(NULL)
  
  inc_u <- suppressWarnings(sf::st_union(inc_all))
  if (length(inc_u) == 0 || sf::st_is_empty(inc_u)) return(NULL)
  
  # exclude rinks/fields
  ex1 <- fetch("leisure", c("pitch","track","stadium","sports_centre","playground","ice_rink"))
  ex2 <- fetch("sport", c("soccer","baseball","basketball","tennis","ice_hockey","hockey"))
  
  e1 <- union_geom(ex1); e2 <- union_geom(ex2)
  ex_all <- sf::st_sfc(crs = 4326)
  if (!is.null(e1) && length(e1) > 0) ex_all <- c(ex_all, e1)
  if (!is.null(e2) && length(e2) > 0) ex_all <- c(ex_all, e2)
  ex_all <- ex_all[!sf::st_is_empty(ex_all)]
  
  if (length(ex_all) > 0) {
    ex_u <- suppressWarnings(sf::st_union(ex_all))
    if (!sf::st_is_empty(ex_u)) inc_u <- suppressWarnings(sf::st_difference(inc_u, ex_u))
  }
  
  if (length(inc_u) == 0 || sf::st_is_empty(inc_u)) return(NULL)
  
  out <- sf::st_as_sf(sf::st_sfc(inc_u, crs = 4326))
  suppressWarnings(sf::st_write(out, cache_path, delete_dsn = TRUE, quiet = TRUE))
  sf::st_transform(out, crs_out)
}

add_park_hatch <- function(parks_sf, spacing = 30, angle = 35,
                           col = "#A9C8B6", alpha = 0.22, lwd = 1.0) {
  if (is.null(parks_sf) || nrow(parks_sf) == 0) return(invisible(NULL))
  parks_u <- suppressWarnings(sf::st_union(parks_sf))
  if (sf::st_is_empty(parks_u)) return(invisible(NULL))
  
  bb <- sf::st_bbox(parks_u)
  cx <- (bb["xmin"] + bb["xmax"]) / 2
  cy <- (bb["ymin"] + bb["ymax"]) / 2
  
  xs <- seq(bb["xmin"] - 2000, bb["xmax"] + 2000, by = spacing)
  lines <- lapply(xs, function(x0) {
    sf::st_linestring(matrix(c(x0, bb["ymin"]-2000, x0, bb["ymax"]+2000), ncol=2, byrow=TRUE))
  })
  L <- sf::st_sfc(lines, crs = sf::st_crs(parks_sf))
  
  ang <- angle * pi/180
  rot <- function(xy){
    x <- xy[,1]-cx; y <- xy[,2]-cy
    xr <- x*cos(ang) - y*sin(ang)
    yr <- x*sin(ang) + y*cos(ang)
    cbind(xr+cx, yr+cy)
  }
  
  L_rot <- sf::st_sfc(lapply(L, function(g){
    sf::st_linestring(rot(sf::st_coordinates(g)))
  }), crs = sf::st_crs(parks_sf))
  
  hatch <- suppressWarnings(sf::st_intersection(sf::st_as_sf(L_rot), parks_u))
  if (!inherits(hatch, "try-error") && nrow(hatch) > 0) {
    plot(sf::st_geometry(hatch), add = TRUE,
         col = grDevices::adjustcolor(col, alpha.f = alpha),
         lwd = as.numeric(lwd))
  }
  invisible(NULL)
}

make_shadow_depth_shift <- function(B, H, az_sun_deg, alt_deg) {
  if (!terra::compareGeom(B, H, stopOnError = FALSE)) {
    H <- terra::resample(H, B, method = "near")
  }
  
  alt_rad <- alt_deg * pi / 180
  if (!is.finite(alt_rad) || alt_rad <= 0) {
    Sd <- B; Sd[] <- 0
    return(Sd)
  }
  
  az_shadow_deg <- (az_sun_deg + 180) %% 360
  azr <- az_shadow_deg * pi/180
  
  Hn <- H
  Hn[is.na(Hn)] <- 0
  
  hcap <- suppressWarnings(as.numeric(
    terra::global(Hn, fun = function(v) stats::quantile(v, 0.98, na.rm = TRUE, names = FALSE))[1,1]
  ))
  if (!is.finite(hcap) || hcap <= 0) hcap <- 10
  
  cell <- terra::res(B)[1]
  Lmax <- hcap / tan(alt_rad)
  nsteps <- max(1L, min(1600L, as.integer(ceiling(Lmax / cell))))
  stride <- if (nsteps > 900) 3L else if (nsteps > 450) 2L else 1L
  
  ux <- sin(azr)
  uy <- cos(azr)
  
  Occ <- B
  Occ[] <- 0
  
  for (k in seq(1L, nsteps, by = stride)) {
    drop_h <- k * cell * tan(alt_rad)
    dxm <- ux * k * cell
    dym <- uy * k * cell
    
    Sh <- terra::shift(Hn, dx = dxm, dy = dym)
    Sh <- terra::resample(Sh, Occ, method = "near")
    Sh[is.na(Sh)] <- 0
    
    cand <- Sh - drop_h
    cand <- terra::ifel(cand > 0, cand, 0)
    
    Occ <- terra::ifel(cand > Occ, cand, Occ)
  }
  
  Sd <- terra::ifel(B == 1, 0, Occ)
  Sd[is.na(Sd)] <- 0
  Sd
}

bld0   <- sf::st_read(bld_path, quiet = TRUE)
roads0 <- sf::st_read(roads_path, quiet = TRUE)

sf::st_crs(bld0)   <- sf::st_crs(crs_bld)
sf::st_crs(roads0) <- sf::st_crs(crs_roads)

bld0   <- sf_clean_polygons(bld0)
roads0 <- sf::st_make_valid(roads0)
roads0 <- roads0[!sf::st_is_empty(roads0), ]

if (nrow(bld0) == 0) stop("Buildings empty.")
if (nrow(roads0) == 0) stop("Roads empty.")

bb_utm <- utm_bbox_from_lonlat(LON, LAT, zoom_km, crs_out = crs_utm)
aoi_utm_poly <- bbox_poly_sf(bb_utm, crs = crs_utm)

bld_utm_all   <- suppressWarnings(sf::st_transform(bld0,   crs_utm))
roads_utm_all <- suppressWarnings(sf::st_transform(roads0, crs_utm))
rm(bld0, roads0); gc()

bld_utm_all   <- sf_clean_polygons(bld_utm_all)
roads_utm_all <- sf::st_make_valid(roads_utm_all)
roads_utm_all <- roads_utm_all[!sf::st_is_empty(roads_utm_all), ]

bld_utm   <- suppressWarnings(sf::st_intersection(bld_utm_all,   aoi_utm_poly))
roads_utm <- suppressWarnings(sf::st_intersection(roads_utm_all, aoi_utm_poly))

bld_utm   <- bld_utm[!sf::st_is_empty(bld_utm), ]
roads_utm <- roads_utm[!sf::st_is_empty(roads_utm), ]

if (nrow(bld_utm) == 0) stop("0 buildings after crop. Increase zoom_km or check CRS.")

bld_utm   <- suppressWarnings(sf::st_simplify(bld_utm,   dTolerance = 0.15, preserveTopology = TRUE))
roads_utm <- suppressWarnings(sf::st_simplify(roads_utm, dTolerance = 0.60, preserveTopology = TRUE))

cand <- c("height","HEIGHT","H","h","elev","ELEV","building_h","bld_h","Z","z","H_m")
hcol <- cand[cand %in% names(bld_utm)][1]

if (!is.na(hcol) && length(hcol) == 1) {
  hm <- suppressWarnings(as.numeric(bld_utm[[hcol]]))
  hm[!is.finite(hm)] <- NA
} else {
  hm <- rep(NA_real_, nrow(bld_utm))
}

if (all(is.na(hm))) {
  a <- as.numeric(sf::st_area(bld_utm))
  hm <- pmax(6, pmin(50, (sqrt(pmax(a, 1)) / 3.2)))
}
bld_utm$H_m <- pmax(3, ifelse(is.finite(hm), hm, h_default))

r_tmpl <- terra::rast(
  terra::ext(bb_utm$xmin, bb_utm$xmax, bb_utm$ymin, bb_utm$ymax),
  resolution = res_m,
  crs = crs_utm
)

B <- terra::rasterize(terra::vect(bld_utm),   r_tmpl, field = 1, background = 0)
H <- terra::rasterize(terra::vect(bld_utm),   r_tmpl, field = "H_m", fun = "max", background = 0)
R <- terra::rasterize(terra::vect(roads_utm), r_tmpl, field = 1, background = 0)

B[is.na(B)] <- 0
H[is.na(H)] <- 0
R[is.na(R)] <- 0

Bpoly <- terra::as.polygons(B, dissolve = TRUE, values = TRUE, na.rm = TRUE)
Bpoly <- Bpoly[Bpoly$layer == 1, ]
Bpoly_sf <- sf::st_make_valid(sf::st_as_sf(Bpoly))

parks_utm <- NULL
Pmask_soft <- NULL
P_cols <- NULL
P_breaks <- seq(0, 1, length.out = 257)

if (isTRUE(add_osm_green)) {
  cache_path <- file.path(out_dir, sprintf("OSM_green_SUBTRACT_zoom_%dm.gpkg",
                                           as.integer(round(zoom_km*1000))))
  try(osmdata::set_overpass_url("https://overpass.kumi.systems/api/interpreter"), silent = TRUE)
  if (is.null(osmdata::get_overpass_url())) {
    try(osmdata::set_overpass_url("https://overpass-api.de/api/interpreter"), silent = TRUE)
  }
  
  parks_utm <- get_osm_green_cached(LON, LAT, zoom_km = zoom_km, crs_out = crs_utm,
                                    cache_path = cache_path, timeout = 240)
  
  if (!is.null(parks_utm) && nrow(parks_utm) > 0) {
    Pmask <- terra::rasterize(terra::vect(parks_utm), r_tmpl, field = 1, background = 0)
    Pmask[is.na(Pmask)] <- 0
    Pmask <- terra::ifel(Pmask > 0, 1, 0)
    
    sigmaP <- max(1.0, park_edge_soft_m / res_m)
    kP <- odd_int(as.integer(round(sigmaP * 6)), 21L)
    KP <- gauss_kernel(kP, sigmaP)
    
    Pmask_soft <- terra::focal(Pmask, w = KP, fun = "mean", na.rm = TRUE, fillvalue = 0)
    Pmask_soft[is.na(Pmask_soft)] <- 0
    Pmask_soft <- clamp01_rast(Pmask_soft)
    
    P_cols <- alpha_ramp(park_fill, alpha_max = park_alpha, n = 256L, pow = 0.35)
    rm(Pmask, sigmaP, kP, KP)
  }
}

bld_buf_cells <- odd_int(bld_buf_cells, 7L)
Kb <- matrix(1, nrow = bld_buf_cells, ncol = bld_buf_cells)
Bbuf <- terra::focal(B, w = Kb, fun = "max", na.rm = TRUE, fillvalue = 0)
Bbuf[is.na(Bbuf)] <- 0
Bbuf <- terra::ifel(Bbuf > 0, 1, 0)

road_blur_cells <- odd_int(road_blur_cells, 3L)
Kr <- matrix(1, nrow = road_blur_cells, ncol = road_blur_cells)
R_base <- terra::focal(R, w = Kr, fun = "max", na.rm = TRUE, fillvalue = 0)
R_base[is.na(R_base)] <- 0
R_base <- terra::ifel(R_base > 0, 1, NA)

sigma_cells2 <- max(1.0, sigma_m_soft / res_m)
k2 <- odd_int(as.integer(round(sigma_cells2 * 6)), 15L)
K2 <- gauss_kernel(k2, sigma_cells2)

shadow_pal <- grDevices::colorRampPalette(
  c("#FFF8F0", "#FEE9D6", "#FDD6B5", "#FDBB84", "#F08A4B", "#D85A1A")
)(256)
shadow_breaks <- seq(0, 1, length.out = 257)

draw_panel <- function(hour0) {
  t0 <- as.POSIXct(sprintf("%s %02d:00:00", as.character(date0), hour0), tz = tz0)
  pos <- suncalc::getSunlightPosition(date = t0, lat = LAT, lon = LON)
  
  alt_deg <- pos$altitude * 180 / pi
  az_deg  <- (pos$azimuth  * 180 / pi + 180) %% 360
  
  Sd <- make_shadow_depth_shift(B, H, az_deg, alt_deg)
  
  Sd <- terra::focal(Sd, w = matrix(1, 3, 3), fun = "max", na.rm = TRUE, fillvalue = 0)
  Sd[is.na(Sd)] <- 0
  
  S <- terra::focal(Sd, w = K2, fun = "mean", na.rm = TRUE, fillvalue = 0)
  S[is.na(S)] <- 0
  
  cap <- suppressWarnings(as.numeric(
    terra::global(S, fun = function(v) stats::quantile(v, shadow_q, na.rm = TRUE, names = FALSE))[1,1]
  ))
  if (!is.finite(cap) || cap <= 0) cap <- 1e-6
  
  S <- clamp01_rast(S / cap)
  S <- S^0.85
  S <- terra::ifel(S == 0, 0, S^shadow_gamma)
  
  if (!is.null(Pmask_soft)) {
    S_mask <- terra::ifel((B == 1) | (Pmask_soft > park_mask_thresh), NA, S)
  } else {
    S_mask <- terra::ifel(B == 1, NA, S)
  }
  
  base <- B; base[] <- 1
  plot(base, col = "white", legend = FALSE, axes = FALSE, box = FALSE)
  
  if (!is.null(Pmask_soft) && !is.null(P_cols)) {
    plot(Pmask_soft, add = TRUE, col = P_cols, breaks = P_breaks, legend = FALSE)
    if (isTRUE(park_hatch_on) && !is.null(parks_utm) && nrow(parks_utm) > 0) {
      add_park_hatch(parks_utm, spacing = 30, angle = 35,
                     col = park_hatch_col, alpha = park_hatch_alpha, lwd = park_hatch_lwd)
    }
  }
  
  pal_alpha <- grDevices::adjustcolor(shadow_pal, alpha.f = shadow_alpha)
  plot(S_mask, add = TRUE, col = pal_alpha, breaks = shadow_breaks, legend = FALSE)
  
  plot(R_base, add = TRUE,
       col = grDevices::adjustcolor(road_col_base, alpha.f = road_alpha_base),
       legend = FALSE)
  
  plot(sf::st_geometry(Bpoly_sf), add = TRUE, border = bld_edge, lwd = bld_lwd)
  
  mtext(sprintf("%02d:00", hour0), side = 1, line = 0.6, adj = 0.5, cex = 0.95)
  
  rm(Sd, S, S_mask, cap, base, pos, t0, alt_deg, az_deg)
  invisible(NULL)
}

n <- length(hours_vec)
nrow_pan <- ceiling(n / ncol_pan)

png(out_png_grid,
    width  = px_panel * ncol_pan,
    height = px_panel * nrow_pan,
    res = dpi0, bg = "white")

par(mfrow = c(nrow_pan, ncol_pan),
    mar = c(2.6, 0.2, 0.4, 0.2),
    oma = c(0, 0, 2.6, 0))

for (i in seq_len(nrow_pan * ncol_pan)) {
  if (i <= n) draw_panel(hours_vec[i]) else plot.new()
  gc()
}

mtext(sprintf("Hourly summer sunlight inequality — | %s | zoom=%.2f km",
              date0, zoom_km),
      outer = TRUE, cex = 1.15, font = 2)

dev.off()
cat("SAVED ->", out_png_grid, "\n")

rm(bld_utm_all, roads_utm_all, bld_utm, roads_utm, aoi_utm_poly, bb_utm,
   r_tmpl, B, H, R, Bbuf, R_base, K2, Kb, Kr, parks_utm, Pmask_soft,
   shadow_pal, shadow_breaks, Bpoly, Bpoly_sf, P_cols, P_breaks, hcol, hm)
gc()
terra::tmpFiles(remove = TRUE)


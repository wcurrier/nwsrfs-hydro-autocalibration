#Verify rsnwelev works standalone
library(rfchydromodels)
library(data.table)

pars <- fread("runs/2zone/SAKW1/pars_default.csv")
ae_tbl <- read.csv("runs/2zone/SAKW1/area_elev_curve.csv")

# Check pxtemp and talr exist
print(pars[name %in% c("pxtemp", "talr", "elev")])

# Load forcing
forcing <- list()
forcing[["SAKW1-1"]] <- fread("runs/2zone/SAKW1/forcing_por_SAKW1-1.csv")
forcing[["SAKW1-2"]] <- fread("runs/2zone/SAKW1/forcing_por_SAKW1-2.csv")

# Run rsnwelev
forcing_adj <- rsnwelev(forcing, pars, ae_tbl)

# Check ptps was replaced with sensible values
cat("Original ptps range zone 1:", range(forcing[[1]]$ptps), "\n")
cat("rsnwelev ptps range zone 1:", range(forcing_adj[[1]]$ptps), "\n")
cat("Original ptps range zone 2:", range(forcing[[2]]$ptps), "\n")
cat("rsnwelev ptps range zone 2:", range(forcing_adj[[2]]$ptps), "\n")

# Spot check: cold day should have ptps near 1, warm day near 0
cat("\nSample values (first 20 steps):\n")
print(data.frame(
  mat = forcing[[1]]$mat_degc[1:20],
  ptps_orig = forcing[[1]]$ptps[1:20],
  ptps_rsnwelev = forcing_adj[[1]]$ptps[1:20]
))

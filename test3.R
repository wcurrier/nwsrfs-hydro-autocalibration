# Vertify pxtemp sensitivity
library(rfchydromodels)
library(data.table)
library(hydroGOF)
source("wrappers.R")
source("obj_fun.R")

basin_dir <- "runs/2zone/SAKW1"
pars <- fread(file.path(basin_dir, "pars_default.csv"))
limits <- fread(file.path(basin_dir, "pars_limits.csv"))
ae_tbl <- read.csv(file.path(basin_dir, "area_elev_curve.csv"))
forcing <- list()
for (z in sort(grep("SAKW1-", unique(pars$zone), value=TRUE))) {
  forcing[[z]] <- fread(file.path(basin_dir, paste0("forcing_por_", z, ".csv")))
}

# Run with pxtemp = -1 (more snow) vs pxtemp = 3 (less snow)
p_mid <- (limits$lower + limits$upper) / 2
names(p_mid) <- limits$name

p_cold <- p_mid
p_cold["pxtemp_SAKW1-1"] <- -1; p_cold["pxtemp_SAKW1-2"] <- -1

p_warm <- p_mid
p_warm["pxtemp_SAKW1-1"] <- 3; p_warm["pxtemp_SAKW1-2"] <- 3

obs_daily <- fread(file.path(basin_dir, "flow_daily_SAKW1.csv"))[, -"Source"]
obs_daily <- obs_daily |> dplyr::right_join(
  forcing[[1]] |> as.data.table() |> dplyr::group_by(year,month,day) |>
  dplyr::summarise(.groups="drop"), by=c("year","month","day")) |> as.data.table()

run_cold <- model_wrapper(p_cold, names(p_mid), 6, pars, obs_daily, NULL,
                          forcing, NULL, "nselognse_NULL", 2, character(0),
                          ae_tbl=ae_tbl, return_flow=TRUE)
run_warm <- model_wrapper(p_warm, names(p_mid), 6, pars, obs_daily, NULL,
                          forcing, NULL, "nselognse_NULL", 2, character(0),
                          ae_tbl=ae_tbl, return_flow=TRUE)

cat("pxtemp=-1 mean flow:", mean(run_cold$sim_flow_cfs, na.rm=TRUE), "\n")
cat("pxtemp= 3 mean flow:", mean(run_warm$sim_flow_cfs, na.rm=TRUE), "\n")
cat("Difference:", mean(run_cold$sim_flow_cfs - run_warm$sim_flow_cfs, na.rm=TRUE), "\n")

if (abs(mean(run_cold$sim_flow_cfs - run_warm$sim_flow_cfs, na.rm=TRUE)) > 1) {
  cat("PASS: pxtemp affects model output\n")
} else {
  cat("FAIL: pxtemp has no effect - check rsnwelev integration\n")
}

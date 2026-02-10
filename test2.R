# Single model_wrapper call (no optimization)
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

obs_daily <- fread(file.path(basin_dir, "flow_daily_SAKW1.csv"))[, -"Source"]
obs_daily <- obs_daily |>
  dplyr::right_join(
    forcing[[1]] |> as.data.table() |>
      dplyr::group_by(year, month, day) |>
      dplyr::summarise(.groups="drop") ,
    by = c("year", "month", "day")) |>
  as.data.table()

# Use midpoint of all parameter bounds as a test
lower <- limits$lower; names(lower) <- limits$name
upper <- limits$upper; names(upper) <- limits$name
p_test <- (lower + upper) / 2

result <- model_wrapper(
  p = p_test, p_names = names(lower),
  dt_hours = 6, default_pars = pars,
  obs_daily = obs_daily, obs_inst = NULL,
  forcing_raw = forcing, upflow = NULL,
  obj_fun = "nselognse_NULL", n_zones = 2,
  cu_zones = character(0),
  ae_tbl = ae_tbl,
  return_flow = TRUE
)
cat("Success! Flow range:", range(result$sim_flow_cfs, na.rm=TRUE), "\n")
cat("Rows:", nrow(result), "\n")


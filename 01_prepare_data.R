# ------------------------------------------------------------------
# 01_prepare_data.R  -  Unemployment in KwaZulu-Natal municipalities (Census 2011)
# Joins the census table to the municipal boundaries. Run ONCE.
# Creates data/municipalities.rds, which app.R reads.
# Boundaries: Municipal Demarcation Board 2011 local municipalities (KZN only, 51)
# Census:     Stats SA Census 2011 Municipal Fact Sheet
# ------------------------------------------------------------------
library(sf)
library(dplyr)

muni   <- st_read("data/kzn_boundaries.geojson")           # muni_code, muni_name, district
census <- read.csv("data/census2011_kzn_municipal.csv")    # one row per municipality

muni_data <- muni %>% inner_join(census, by = "muni_code")  # join on the CODE, not the name

nrow(muni_data)                                  # must be 51
setdiff(census$muni_code, muni_data$muni_code)   # must be character(0)

muni_data <- muni_data %>% st_make_valid() %>% st_transform(4326)
saveRDS(muni_data, "data/municipalities.rds")

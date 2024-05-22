### Load libraries
library(tidyverse)
library(rstan)
library(reshape2)
library(xlsx)

here::i_am("R/002_display_reg_noid.R")
outfile <- here::here("outputs", "regnoid_workbook_v2.xlsx")
alpha <- 1/40

load(here::here("working",
                "mrp_samples_bivreg_2024-04-25.RData"))

### Post-process ward counts
reg_with_id_per_ward <- rstan::extract(samples,
                                 pars = "consty_reg_with_id",
                                 permuted = TRUE)[[1]]

### Melt
reg_with_id_per_ward <- reshape2::melt(reg_with_id_per_ward)
names(reg_with_id_per_ward) <- c("iter", "ward", "reg_with_id_count")
reg_with_id_per_ward$ward <- consty_levels[reg_with_id_per_ward$ward]

### Add on the actual counts
psf <- readRDS(here::here("data",
                          "ward_joint.rds")) |>
    filter(Age != "16 to 17")

tots <- psf |>
    group_by(geogcode) |>
    summarize(tot = sum(count))

reg_with_id_per_ward <- left_join(reg_with_id_per_ward,
                                    tots,
                                    by = join_by(ward == geogcode))

reg_noid_per_ward <- rstan::extract(samples,
                                 pars = "consty_reg_no_id",
                                 permuted = TRUE)[[1]]

### Melt
reg_noid_per_ward <- reshape2::melt(reg_noid_per_ward)
names(reg_noid_per_ward) <- c("iter", "ward", "reg_noid_count")
reg_noid_per_ward$ward <- consty_levels[reg_noid_per_ward$ward]

comb <- left_join(reg_with_id_per_ward, reg_noid_per_ward,
                  by = join_by(iter, ward))


noid_per_ward <- rstan::extract(samples,
                                pars = "consty_noid",
                                permuted = TRUE)[[1]]

### Melt
noid_per_ward <- reshape2::melt(noid_per_ward)
names(noid_per_ward) <- c("iter", "ward", "noid_count")
noid_per_ward$ward <- consty_levels[noid_per_ward$ward]


comb <- left_join(comb, noid_per_ward,
                  by = join_by(iter, ward))

### Probably need to add the LAD on to this
aux <- readRDS(here::here("working", "tidy_aux_data.rds"))

comb <- left_join(comb,
                  aux |> dplyr::select(geogcode,
                                       wardname, LAD23NM),
                  by = join_by(ward == geogcode))

### 
### Calculate proportions
comb <- comb |>
    mutate(Pct_Registered = (reg_with_id_count + reg_noid_count) / tot,
           Reg_NoID_prop = reg_noid_count / (reg_noid_count + reg_with_id_count))

comb <- comb |>
    dplyr::select(Borough = LAD23NM,
                  Ward = wardname,
                  `Ward code` = ward,
                  everything())

### What do we want
t1 <- comb |>
    group_by(Ward, `Ward code`, Borough) |>
    summarize(`Avg Pct registered` = round(100 * mean(Pct_Registered), 1),
              `VAP lo` = floor(100 * quantile(Pct_Registered, alpha)),
              `VAP hi` = ceiling(100 * quantile(Pct_Registered, 1 - alpha)),
              `Avg registered w/o ID` = round(100 * mean(Reg_NoID_prop), 1),
              `Registered w/o ID, lo` = round(100 * quantile(Reg_NoID_prop, alpha), 1),
              `Registered w/o ID, hi` = round(100 * quantile(Reg_NoID_prop, 1 - alpha), 1),
              .groups = "drop")


class(t1) <- "data.frame"

write.xlsx(t1,
           outfile,
           sheetName = "Proportion lacking ID by ward",
           append = FALSE,
           row.names = FALSE)


## vars <- c("Avg VAP turnout",
##           "Avg voters w/o ID",
##           "Avg no ID turnout")
## pairs(t1[,vars])

### What's the correlation between the two sets?
t0 <- read.xlsx(file = here::here("outputs",
                                  "regnoid_workbook.xlsx"),
                sheetIndex = 1,
                startRow = 2)

comb <- full_join(t1 |> dplyr::select(Ward, `Ward code`, `Avg registered w/o ID`),
                  t0 |> dplyr::select(Ward, Ward.code, Best.estimate.1),
                  by = join_by(Ward, `Ward code` == Ward.code))

cor(comb[,3], comb[,4])

ggplot(comb, aes(x = `Best.estimate.1`,
                 y = `Avg registered w/o ID`)) +
    geom_point() +
    geom_smooth(method = "lm", se = FALSE) +
    geom_abline(slope = 1, intercept = 0) +
    theme_bw()

summary(lm(`Avg registered w/o ID` ~ Best.estimate.1, data = comb))

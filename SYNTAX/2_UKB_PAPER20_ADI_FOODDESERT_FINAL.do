
capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\DATA_MANAGEMENT.smcl",replace


**Food desert data**
**https://data.cdrc.ac.uk/dataset/e-food-desert-index**
**https://data.hasp.ac.uk/browser/dataset/5347/0**


**Index of multiple deprivation-imd**
**https://data.cdrc.ac.uk/dataset/index-multiple-deprivation-imd**
**https://data.geods.ac.uk/dataset/index-of-multiple-deprivation-imd


**Shapefile: Great Britain 2011**
**https://statistics.ukdataservice.ac.uk/dataset/2011-census-geography-boundaries-lower-layer-super-output-areas-and-data-zones**

**Safeguarded data**
**https://data.cdrc.ac.uk/dataset/uk-lsoa-dz-sdz-classification-20212-lsoac/resource/variable-dictionary-classification-codes#{currentView:!grid,view-grid:{columnsWidth:[{column:!Classification++Name,width:722}]}}

**https://data.cdrc.ac.uk/dataset/uk-lsoa-dz-sdz-classification-20212-lsoac?utm_source

**https://data.geods.ac.uk/dataset/lsoac



*-----------------------------*
* Step 0: Set working directory
*-----------------------------*
cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA"

*-----------------------------*
* Step 1: Load Food Desert Data
*-----------------------------*
use UK_fooddesertdatafinalized, clear

* Clean LSOA codes (ensure uppercase, no trailing spaces)
replace lsoa11 = upper(strtrim(lsoa11))

* Save cleaned file
save UK_fooddesertdatafinalized_clean, replace

*-----------------------------*
* Step 2: Load IMD Data
*-----------------------------*
import delimited "Dataset1_CDRC_Harmonised_IMD2019.csv", clear

* Create `lsoa11` variable for merging
gen lsoa11 = upper(strtrim(lsoa))  // Clean LSOA codes
drop lsoa
save IMD2019_clean, replace

*-----------------------------*
* Step 3: Merge IMD with Food Desert Data
*-----------------------------*
use IMD2019_clean, clear
merge 1:1 lsoa11 using UK_fooddesertdatafinalized_clean

* Keep matched data
keep if _merge == 3
drop _merge

*-----------------------------*
* Step 4: Standardize and Transform Scores
*-----------------------------*
egen zfd = std(fd)
egen zsoa_pct = std(soa_pct)
gen zsoa_pctinv = -1 * zsoa_pct   // Invert so high = worse for deprivation

keep if zfd~=. & zsoa_pctinv~=. 


save UK_fooddesertdatafinalizedIMD, replace

su zfd zsoa_pctinv

histogram zfd
histogram zsoa_pctinv


save UK_fooddesertdatafinalizedIMD, replace


*-----------------------------*
* Step 5: Create Cluster Categories
*-----------------------------*
xtile zsoa_pctinv_cluster = zsoa_pctinv, nq(4)
xtile zfd_cluster = zfd, nq(4)

gen cluster_soa_pctinvfd = .
replace cluster_soa_pctinvfd = 1 if zsoa_pctinv_cluster >= 3 & zfd_cluster >= 3  &  zsoa_pctinv_cluster~=. &  zfd_cluster~=. // High-High
replace cluster_soa_pctinvfd = 2 if zsoa_pctinv_cluster <= 2 & zfd_cluster <= 2  &  zsoa_pctinv_cluster~=. &  zfd_cluster~=.   // Low-Low
replace cluster_soa_pctinvfd = 3 if zsoa_pctinv_cluster > 2 & zfd_cluster <= 2   &  zsoa_pctinv_cluster~=. &  zfd_cluster~=.   // Mixed High
replace cluster_soa_pctinvfd = 4 if zsoa_pctinv_cluster <= 2 & zfd_cluster > 2   &  zsoa_pctinv_cluster~=. &  zfd_cluster~=.   // Mixed Low
replace cluster_soa_pctinvfd = 0 if missing(cluster_soa_pctinvfd)                   // Non-significant

tab cluster_soa_pctinvfd

label define cluster_soa_pctinvfd 0 "Non-significant" 1 "High-High" 2 "Low-Low" 3 "Mixed High" 4 "Mixed Low"
label values cluster_soa_pctinvfd cluster_soa_pctinvfd

tab cluster_soa_pctinvfd

save UK_fooddesertdatafinalizedIMD, replace


*-----------------------------*
* Step 6: Convert SHP to Stata format (if not already done)
*-----------------------------*
* This step requires that infuse_gb_2011_clipped.shp and associated .dbf and .shx files are present in this directory

shp2dta using infuse_lsoa_lyr_2011_clipped.shp, database(uk_db) coordinates(uk_coords) replace


*-----------------------------*
* Step 7: Prepare shapefile attribute data for merging
*-----------------------------*
* Load the attribute data
use uk_db, clear

* Rename and clean LSOA code
rename geo_code lsoa11
replace lsoa11 = upper(strtrim(lsoa11))

save uk_db_str, replace

*-----------------------------*
* Step 8: Merge cluster data with shapefile data
*-----------------------------*

use UK_fooddesertdatafinalizedIMD, clear
replace lsoa11 = upper(strtrim(lsoa11))
merge 1:1 lsoa11 using uk_db_str
tab _merge


*-----------------------------*
* Step 9: Check cluster distribution
*-----------------------------*
tab cluster_soa_pctinvfd, missing

*-----------------------------*
* Step 10: Map the Clusters
*-----------------------------*
* Note: ID variable must be numeric, use _ID from shapefile

**0 "Non-significant" 1 "High-High" 2 "Low-Low" 3 "Mixed High" 4 "Mixed Low"

capture drop cluster_ordered
gen cluster_ordered = .
replace cluster_ordered = 2 if cluster_soa_pctinvfd == 1   // High-High
replace cluster_ordered = 0 if cluster_soa_pctinvfd == 2   // Low-Low
replace cluster_ordered = 1 if cluster_soa_pctinvfd == 3   // Mixed-High
replace cluster_ordered = 1 if cluster_soa_pctinvfd == 4   // Mixed-Low
replace cluster_ordered = . if cluster_soa_pctinvfd == 0   // Mixed-Low




capture label drop cluster_lbl
label define cluster_lbl 0 "Low-Low" 1 "M-H, M-L" 2 "High-High", replace
label values cluster_ordered cluster_lbl

tab cluster_ordered

save UK_fooddesertdatafinalizedIMD, replace



clear
clear matrix
clear mata
set maxvar 20000
set memory 2g

use UK_fooddesertdatafinalizedIMD,clear

spmap cluster_ordered using uk_coords, id(_ID) ///
    fcolor("0 102 204" "211 211 211" "204 0 0") clmethod(unique)

		
*-----------------------------*
* Step 11: Save Map Output
*-----------------------------*
graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1A.gph", replace
graph export "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1A.png", width(2000) replace


*-----------------------------*
* Step 12: Create individual zfd and zsoa_pctinv maps and Save Map Output
*-----------------------------*

clear
clear matrix
clear mata
set maxvar 20000
set memory 2g

use UK_fooddesertdatafinalizedIMD,clear

spmap zfd using uk_coords, id(_ID) ///
    fcolor(Rainbow) clmethod(kmeans)


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1B.gph", replace
graph export "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1B.png", width(2000) replace

clear
clear matrix
clear mata
set maxvar 20000
set memory 2g

use UK_fooddesertdatafinalizedIMD,clear


spmap zsoa_pctinv using uk_coords, id(_ID) ///
    fcolor(Rainbow) clmethod(kmeans)


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1C.gph", replace
graph export "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE1C.png", width(2000) replace

	
	*-----------------------------*
* Step 13: Save Final Dataset
*-----------------------------*
save UK_fooddesertdatafinalizedIMD_merged, replace


use UK_fooddesertdatafinalizedIMD_merged,clear

tab zsoa_pctinv_cluster zfd_cluster, row col chi 

mlogit  zfd_cluster zsoa_pctinv_cluster, rrr

mlogit  zfd_cluster i.zsoa_pctinv_cluster, rrr

bysort  zsoa_pctinv_cluster: su zsoa_pctinv

bysort  zfd_cluster: su zfd

reg zfd zsoa_pctinv

lpoly zfd zsoa_pctinv

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE2A.gph",replace

lpoly zfd zsoa_pctinv if Country=="England and Wales"

reg zfd zsoa_pctinv if Country=="England and Wales"

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE2B.gph",replace

reg zfd zsoa_pctinv if Country=="Scottland"

lpoly zfd zsoa_pctinv if Country=="Scottland"

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FIGURES\FIGURE2C.gph",replace



capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\DATA_MANAGEMENT_LSOA_LIST.txt", text replace


cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA"

use UK_fooddesertdatafinalizedIMD_merged,clear

bysort cluster_ordered: list lsoa11


**Get Area in R**

**library(sf)
**library(dplyr)
**library(haven)

**shp_path <- "C:/Users/jmw596/Downloads/infuse_lsoa_lyr_2011_clipped/infuse_lsoa_lyr_2011_clipped.shp"   # or full path to the .shp

**gdf <- st_read(shp_path, quiet = FALSE)
**print(names(gdf))
**gdf <- st_transform(gdf, 27700)

**gdf <- as.data.frame(gdf %>%
**  mutate(area_km2 = as.numeric(st_area(geometry)) / 1e6)) %>%
**  select(-c(geometry))

**write_dta(gdf, "lsoa2011_area_km2.dta")

use lsoa2011_area_km2.dta,clear

capture drop lsoa11
gen lsoa11=geo_code

sort lsoa11

save lsoa2011_area_km2,replace 

use UK_fooddesertdatafinalizedIMD_merged,clear

sort lsoa11

save UK_fooddesertdatafinalizedIMD_merged,replace
capture drop _merge

merge lsoa11 using lsoa2011_area_km2

save UK_fooddesertdatafinalizedIMD_mergedfin,replace

bysort cluster_ordered: su area_km2

**----------------------------------------------------------------------------------------------------------------------------------------
**-> cluster_ordered = Low-Low
**
**    Variable |        Obs        Mean    Std. dev.       Min        Max
**-------------+---------------------------------------------------------
**    area_km2 |     13,487    1.221803    3.283567   .0181287   183.0458
**
**----------------------------------------------------------------------------------------------------------------------------------------
**-> cluster_ordered = M-H, M-L
**
**    Variable |        Obs        Mean    Std. dev.       Min        Max
**-------------+---------------------------------------------------------
**    area_km2 |     14,966    9.196109    34.76696    .015549   1004.007
**
**----------------------------------------------------------------------------------------------------------------------------------------
**-> cluster_ordered = High-High
**
**    Variable |        Obs        Mean    Std. dev.       Min        Max
**-------------+---------------------------------------------------------
**    area_km2 |     13,276    5.652416    29.41317   .0093673    1162.51


save UK_fooddesertdatafinalizedIMD_mergedfin,replace


 graph bar (mean) area_km2, over(nation) over(cluster_ordered)
 
 graph save "FIGURE1_PANELD.gph",replace

capture log close

**R SCRIPT FOR MORE DETAILED PANEL D**

###############################################################################
# FIGURE: 3 panels of mean area (km^2) by Nation ×
#   Panel D.1: EFDI quartiles (zfd_cluster -> zfd_q)
#   Panel D.2: IMD quartiles  (zsoa_pctinv_cluster -> imd_q)
#   Panel D.3: Ordered EFDI/IMD clusters (cluster_ordered -> cl_ord)
#
# Adds Bonferroni-corrected p-values and significance stars:
#   *   p < 0.05
#   **  p < 0.010
#   *** p < 0.001
#
# FIX for your rstatix error:
# - Do NOT use rstatix::anova_test() output directly inside mutate/transmute + p.adjust
# - Use base R summary(aov()) to get a plain numeric p-value per stratum
###############################################################################

# ---- 0) Packages ----
pkgs <- c("haven","dplyr","tidyr","stringr","forcats","ggplot2","patchwork","readr")
to_install <- pkgs[!pkgs %in% rownames(installed.packages())]
if (length(to_install)) install.packages(to_install)
invisible(lapply(pkgs, library, character.only = TRUE))

# ---- 1) Helpers ----
p_to_stars <- function(p){
  dplyr::case_when(
    is.na(p)  ~ "",
    p < 0.001 ~ "***",
    p < 0.010 ~ "**",
    p < 0.05  ~ "*",
    TRUE      ~ ""
  )
}

# safe ANOVA p-value extractor (returns numeric scalar or NA)
anova_p <- function(df, formula){
  # Need >=2 groups and >=2 non-missing rows
  d <- df %>% dplyr::filter(stats::complete.cases(stats::model.frame(formula, data = df)))
  if (nrow(d) < 3) return(NA_real_)
  mm <- stats::model.frame(formula, data = d)
  # response + one factor
  if (ncol(mm) < 2) return(NA_real_)
  g <- mm[[2]]
  if (!is.factor(g)) g <- factor(g)
  if (nlevels(g) < 2) return(NA_real_)
  fit <- stats::aov(formula, data = d)
  sm <- summary(fit)
  # first term p-value
  p <- tryCatch(as.numeric(sm[[1]][["Pr(>F)"]][1]), error = function(e) NA_real_)
  p
}

# ---- 2) File paths (EDIT THESE) ----
infile <- "E:/16GBBACKUPUSB/BACKUP_USB_SEPTEMBER2014/May Baydoun_folder/UK_BIOBANK_PROJECT/UKB_PAPER20_ADI_FOODESERT/DATA/UK_fooddesertdatafinalizedIMD_mergedfin.dta"

out_png <- "E:/16GBBACKUPUSB/BACKUP_USB_SEPTEMBER2014/May Baydoun_folder/UK_BIOBANK_PROJECT/UKB_PAPER20_ADI_FOODESERT/FIGURES/FIGURE_area_3panel_bonferroni.png"

out_csv <- "E:/16GBBACKUPUSB/BACKUP_USB_SEPTEMBER2014/May Baydoun_folder/UK_BIOBANK_PROJECT/UKB_PAPER20_ADI_FOODESERT/OUTPUT/area_bonferroni_pvalues.csv"

# ---- 3) Read data ----
df <- haven::read_dta(infile)

# ---- 4) Harmonize nation variable (SAFE) ----
if ("nation" %in% names(df)) {
  df <- df %>% mutate(nation = as.character(nation))
} else if ("Country" %in% names(df)) {
  df <- df %>% mutate(nation = as.character(Country))
} else {
  stop("No nation or Country variable found. Columns:\n", paste(names(df), collapse = ", "))
}
df <- df %>% mutate(nation = stringr::str_trim(nation))

# OPTIONAL: enforce desired order
# df <- df %>% mutate(nation = factor(nation, levels = c("England","Scotland","Wales")))

# ---- 5) Create factor versions of cluster variables (match Stata) ----
df <- df %>%
  mutate(
    zfd_cluster = as.integer(zfd_cluster),
    zsoa_pctinv_cluster = as.integer(zsoa_pctinv_cluster),
    cluster_ordered = as.integer(cluster_ordered),
    
    zfd_q = factor(zfd_cluster, levels = 1:4,
                   labels = c("EFDI Q1 (lowest)","EFDI Q2","EFDI Q3","EFDI Q4 (highest)")),
    imd_q = factor(zsoa_pctinv_cluster, levels = 1:4,
                   labels = c("IMD Q1 (lowest)","IMD Q2","IMD Q3","IMD Q4 (highest)")),
    cl_ord = factor(cluster_ordered, levels = c(0,1,2),
                    labels = c("Low-Low","M-H, M-L","High-High"))
  )

# ---- 6) Core function to build each panel + Bonferroni stars ----
make_panel <- function(dat, group_var, panel_title){
  
  d <- dat %>%
    filter(!is.na(area_km2), !is.na(nation), !is.na(.data[[group_var]])) %>%
    mutate(
      nation = factor(nation),
      grp    = factor(.data[[group_var]])
    )
  
  # Summary for bars
  sumdat <- d %>%
    group_by(nation, grp) %>%
    summarise(
      n = n(),
      mean = mean(area_km2, na.rm = TRUE),
      se = sd(area_km2, na.rm = TRUE) / sqrt(n),
      .groups = "drop"
    )
  
  # (1) Across nations within each cluster level: p(nation | grp)
  p_nation_by_grp <- d %>%
    group_by(grp) %>%
    summarise(
      p_raw = anova_p(cur_data(), area_km2 ~ nation),
      .groups = "drop"
    ) %>%
    mutate(
      p_bonf = p.adjust(p_raw, method = "bonferroni"),
      stars_nation = p_to_stars(p_bonf)
    )
  
  # Legend labels with stars appended
  grp_levels <- levels(d$grp)
  grp_label_df <- tibble(grp = factor(grp_levels, levels = grp_levels)) %>%
    left_join(p_nation_by_grp, by = "grp") %>%
    mutate(legend_lab = paste0(as.character(grp), stars_nation))
  
  legend_map <- setNames(grp_label_df$legend_lab, as.character(grp_label_df$grp))
  
  # (2) Within each nation across cluster levels: p(grp | nation)
  p_grp_by_nation <- d %>%
    group_by(nation) %>%
    summarise(
      p_raw = anova_p(cur_data(), area_km2 ~ grp),
      .groups = "drop"
    ) %>%
    mutate(
      p_bonf = p.adjust(p_raw, method = "bonferroni"),
      stars_within_nation = p_to_stars(p_bonf)
    )
  
  # y position for nation-level stars
  ystars <- sumdat %>%
    group_by(nation) %>%
    summarise(y = max(mean + se, na.rm = TRUE) * 1.08, .groups = "drop") %>%
    left_join(p_grp_by_nation, by = "nation")
  
  # Plot
  p <- ggplot(sumdat, aes(x = nation, y = mean, fill = grp)) +
    geom_col(position = position_dodge(width = 0.85), width = 0.75) +
    geom_errorbar(
      aes(ymin = mean - se, ymax = mean + se),
      position = position_dodge(width = 0.85),
      width = 0.20
    ) +
    geom_text(
      data = ystars,
      aes(x = nation, y = y, label = stars_within_nation),
      inherit.aes = FALSE,
      vjust = 0,
      size = 5
    ) +
    labs(
      title = panel_title,
      x = NULL,
      y = "Mean area (km^2)",
      fill = "Cluster level\n(legend stars = across nations,\nBonferroni)"
    ) +
    scale_fill_discrete(labels = legend_map) +
    theme_classic(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold"),
      legend.position = "right"
    )
  
  list(
    plot = p,
    nation_diff_within_cluster = p_nation_by_grp %>%
      mutate(panel = panel_title, test = "Across nations within cluster"),
    cluster_diff_within_nation = p_grp_by_nation %>%
      mutate(panel = panel_title, test = "Within nation across clusters")
  )
}

# ---- 7) Build three panels ----
res_A <- make_panel(df, "zfd_q",  "Panel D.1 Area by Nation × EFDI quartiles")
res_B <- make_panel(df, "imd_q",  "Panel D.2 Area by Nation × IMD quartiles")
res_C <- make_panel(df, "cl_ord", "Panel D.3 Area by Nation × Ordered EFDI/IMD clusters")

# ---- 8) Combine ----
fig_3panel <- (res_A$plot / res_B$plot / res_C$plot) +
  patchwork::plot_annotation(
    caption = paste(
      "Legend stars: Bonferroni-adjusted omnibus p for differences ACROSS NATIONS within each cluster level.",
      "Stars above nations: Bonferroni-adjusted omnibus p for differences ACROSS CLUSTER LEVELS within each nation.",
      "(* p<0.05; ** p<0.010; *** p<0.001).",
      sep = "\n"
    )
  )

print(fig_3panel)

# ---- 9) Export figure ----
ggsave(out_png, fig_3panel, width = 10, height = 14, dpi = 300)

# ---- 10) Export p-values ----
pvals_out <- bind_rows(
  res_A$nation_diff_within_cluster,
  res_B$nation_diff_within_cluster,
  res_C$nation_diff_within_cluster,
  res_A$cluster_diff_within_nation,
  res_B$cluster_diff_within_nation,
  res_C$cluster_diff_within_nation
) %>%
  mutate(stars = p_to_stars(p_bonf))

readr::write_csv(pvals_out, out_csv)

###############################################################################
# OPTIONAL: If you want a nonparametric alternative (safer for skewed area):
# - Replace anova_p(...) with kruskal.test(...) p-values in anova_p()
###############################################################################
# kruskal_p <- function(df, formula){
#   d <- df %>% dplyr::filter(stats::complete.cases(stats::model.frame(formula, data = df)))
#   if (nrow(d) < 3) return(NA_real_)
#   mm <- stats::model.frame(formula, data = d)
#   g <- mm[[2]]
#   if (!is.factor(g)) g <- factor(g)
#   if (nlevels(g) < 2) return(NA_real_)
#   out <- tryCatch(kruskal.test(formula, data = d)$p.value, error = function(e) NA_real_)
#   as.numeric(out)
# }
###############################################################################



*************************************

capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE1.smcl",replace

cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA"


use UK_fooddesertdatafinalizedIMD_merged,clear

*******Food desert risk*******
su fd
su fd if Country=="England and Wales"
su fd if Country=="Scottland"

su zfd
su zfd if Country=="England and Wales"
su zfd if Country=="Scottland"

bysort zfd_cluster: su zfd 
bysort zfd_cluster: su zfd if Country=="England and Wales"
bysort zfd_cluster: su zfd if Country=="Scottland"


*******IMD**************

su soa_pct
su soa_pct if Country=="England and Wales"
su soa_pct if Country=="Scottland"

su zsoa_pctinv
su zsoa_pctinv if Country=="England and Wales"
su zsoa_pctinv if Country=="Scottland"

bysort zsoa_pctinv_cluster: su zsoa_pctinv 
bysort zsoa_pctinv_cluster: su zsoa_pctinv if Country=="England and Wales"
bysort zsoa_pctinv_cluster: su zsoa_pctinv if Country=="Scottland"

******Clustering of IMD and Food desert index**************

tab cluster_soa_pctinvfd
tab cluster_soa_pctinvfd if Country=="England and Wales"
tab cluster_soa_pctinvfd if Country=="Scottland"

save UK_fooddesertdatafinalizedIMD_merged, replace

capture log close


********************************INDIVIDUAL COMPONENTS OF THE IMD*********************************

cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\IndexofMultipleDeprivation(IMD)_UK\"

capture log close

capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\DATA_MANAGEMENT2",replace

**STEP 1: OPEN CSV FILES AND COVERT TO DTA FILES**

**Dataset 1

import delimited Dataset1_CDRC_Harmonised_IMD2019.csv,clear

describe 

save Dataset1_CDRC_Harmonised_IMD2019, replace

**Dataset 2

import delimited Dataset2_CDRC_London_IMD2019.csv, clear

describe

save Dataset2_CDRC_London_IMD2019,replace

**Dataset 3

import delimited Dataset3_English_IMD2010_2011LOAs_2015.csv, clear

describe

save Dataset3_English_IMD2010_2011LOAs_2015,replace

**Dataset 4

import delimited Dataset4_English_IMD2004_2007_2010.csv, clear

describe

save Dataset4_English_IMD2004_2007_2010,replace


import delimited Dataset5_NorthernIrish_MDM2005_2010_2017.csv, clear

describe

save Dataset5_NorthernIrish_MDM2005_2010_2017,replace


import delimited Dataset6_NorthernIrish_MDM2001.csv, clear

describe

save Dataset6_NorthernIrish_MDM2001,replace


import delimited Dataset7_Scottish_IMD_SIMD_2016_2020.csv, clear

describe

save Dataset7_Scottish_IMD_SIMD_2016_2020,replace


import delimited Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012.csv, clear

describe

save Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012,replace


import delimited Dataset9_Welsh_IMD_WIMD_2014_2019.csv, clear

describe 

save Dataset9_Welsh_IMD_WIMD_2014_2019,replace


import delimited Dataset10_Welsh_IMD_WIMD_2005_2008_2011.csv, clear

describe

save Dataset10_Welsh_IMD_WIMD_2005_2008_2011,replace

save, replace


**STEP 2: OPEN DTA FILES**

use Dataset1_CDRC_Harmonised_IMD2019, clear 

use Dataset2_CDRC_London_IMD2019, clear 

use Dataset3_English_IMD2010_2011LOAs_2015, clear

use Dataset4_English_IMD2004_2007_2010, clear

use Dataset5_NorthernIrish_MDM2005_2010_2017, clear

use Dataset6_NorthernIrish_MDM2001, clear

use Dataset7_Scottish_IMD_SIMD_2016_2020, clear

use Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012, clear

use Dataset9_Welsh_IMD_WIMD_2014_2019, clear

use Dataset10_Welsh_IMD_WIMD_2005_2008_2011, clear

**Locate Common Varialbe**

* Install findname if not installed
ssc install findname

* Find variables with specific properties
findname, local(in1)
use Dataset2_CDRC_London_IMD2019.dta,clear
findname, local(in2)
use Dataset3_English_IMD2010_2011LOAs_2015.dta,clear
findname, local(in3)
use Dataset4_English_IMD2004_2007_2010.dta,clear
findname, local(in4)
use Dataset5_NorthernIrish_MDM2005_2010_2017.dta,clear
findname, local(in5)
use Dataset6_NorthernIrish_MDM2001.dta,clear
findname, local(in6)
use Dataset7_Scottish_IMD_SIMD_2016_2020.dta,clear
findname, local(in7)
use Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012.dta,clear
findname, local(in8)
use Dataset9_Welsh_IMD_WIMD_2014_2019.dta,clear
findname, local(in9)
use Dataset10_Welsh_IMD_WIMD_2005_2008_2011.dta,clear
findname, local(in10)

//NO COMMON VARIALBE AMONG ALL DATASETS//

***STEP 3: COMBINE DATASETS TOGETHER**

//Merge English Datasets, Wide//

use Dataset3_English_IMD2010_2011LOAs_2015, clear
capture rename lsoacode2011 lsoa 
capture rename v34 IDAC_Rank
capture rename v35 IDACI_Decile
capture rename v37 IDAOPI_Rank
capture rename v38 IDAOPI_Decile
capture rename incomedeprivationaffectingolderp IDAOPI_Score
capture  rename incomedeprivationaffectingchildr IDACI_Score
save, replace
sort lsoa 
save, replace

use Dataset2_CDRC_London_IMD2019, clear
capture rename ls11cd lsoa 
save, replace
sort lsoa 
save, replace

use Dataset4_English_IMD2004_2007_2010, clear
capture rename v22 IMD2007_SCORE
capture rename v23 RANK_OF_IMD2007
capture rename v24 INCOME2007_SCORE
capture rename v25 RANK_OF_INCOME2007_SCORE
capture rename v26 EMPLOYMENT2007_SCORE
capture rename v27 RANK_OF_EMPLOYMENT2007_SCORE
capture rename v28 HEALTH_DEPDIS2007_SCORE 
capture rename v29 RANK_OF_HEALTH_DEPDIS2007_SCORE
capture rename v30 EDUC_AND_TRAINING2007_SCORE
capture rename v31 RANK_EDUC_AND_TRAINING2007_SCORE
capture rename v32 BARRI_HOUSE_AND_SERV2007_SCORE
capture rename v33 RANK_BARRI_HOUSE_SERV2007_SCORE
capture rename v34 CRIME_AND_DISORDER2007_SCORE
capture rename v35 RANKCRIME_AND_DISORDER2007_SCORE
capture rename v36 LIVING_ENVIRONMENT2007_SCORE
capture rename v37 RANKLIVING_ENVIRONMENT2007_SCORE

sort lsoa
save, replace


use Dataset3_English_IMD2010_2011LOAs_2015, clear
merge lsoa using Dataset4_English_IMD2004_2007_2010
capture drop _merge
sort lsoa
merge lsoa using Dataset2_CDRC_London_IMD2019
capture drop _merge
sort lsoa
save English_Data, replace


//Merge Northern Irish Datasets, Wide//

use Dataset5_NorthernIrish_MDM2005_2010_2017, clear
capture rename v61 rankincomedeprivaffectingolder
capture rename v33 RankofHealthDepriDisabDomScore
capture rename v35 RankofEducSkillsTrainDomScore
sort lgdname 
save, replace

use Dataset6_NorthernIrish_MDM2001, clear
sort lgdname 
save, replace

use Dataset5_NorthernIrish_MDM2005_2010_2017, clear
merge lgdname using Dataset6_NorthernIrish_MDM2001
capture drop _merge
sort lgdname
save NorthernIrish_Dataset, replace

//Merge Scottish Datasets, Wide//

use Dataset7_Scottish_IMD_SIMD_2016_2020, clear
capture rename data_zone datazone
capture rename v48 Total_population2020
capture rename v62 income_rate2020
capture rename v63 income_count2020
capture rename v64 employment_rate2020
capture rename v65 employment_count2020
capture rename v66 CIF2020
capture rename v67 ALCOHOL2020
capture rename v68 DRUG2020
capture rename v69 SMR2020
capture rename v70 DEPRESS2020
capture rename v71 LBWT2020
capture rename v72 EMERG2020
capture rename v73 Attendance2020
capture rename v74 Attainment2020
capture rename v79 crime_rate2020
capture rename v80 overcrowded_count2020
capture rename v82 overcrowded_rate2020
capture renamev84 drive_petrol2020
capture rename v85 drive_GP2020
capture rename v87 drive_primary2020
capture rename v88 drive_retail2020
capture rename v89 drive_secondary2020
capture rename v90 PT_GP2020
capture rename v91 PT_post2020
capture rename v92 PT_retail2020
save, replace
sort datazone
save, replace


use Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012, clear
capture rename v29 Score2006
capture rename v30 Rank2006
capture rename v31 IncRate2006
capture rename v32 NumIncDep2006
capture rename v33 IncRank2006
capture rename v34 EmpRate2006
capture rename v35 NumEmpDep2006
capture rename v36 EmpRank2006
capture rename v37 HlthScore2006
capture rename v38 HlthRank2006
capture rename v39 EduScore2006
capture rename v40 EduRank2006
capture rename v41 HouseScore2006
capture rename v42 HouseRank2006
capture rename v43 GAccScore2006
capture rename v44 GAccRank2006
capture rename v47 Quintile2006
capture rename v48 Decile2006
capture rename v49 Vigintile2006
capture rename v50 Percentile2006
capture rename v53 Score2009
capture rename v54 Rank2009
capture rename v55 IncRate2009
capture rename v56 NumIncDep2009
capture rename v57 IncRank2009
capture rename v58 EmpRate2009
capture rename v59 NumEmpDep2009
capture rename v60 EmpRank2009
capture rename v61 HlthScore2009
capture rename v62 HlthRank2009
capture rename v63 EduScore2009
capture rename v64 EduRank2009
capture rename v65 HouseScore2009
capture rename v66 HouseRank2009
capture rename v67 GAccScore2009
capture rename v68 GAccRank2009
capture rename v69 CrimeScore2009
capture rename v70 CrimeRank2009
capture rename v71 Quintile2009
capture rename v72 Decile2009
capture rename v73 Vigintile2009
capture rename v74 Percentile2009
capture rename v77 Score2012
capture rename v78 Rank2012
capture rename v79 IncRate2012
capture rename v80 NumIncDep2012
capture rename v81 IncRank2012
capture rename v82 EmpRate2012
capture rename v83 NumEmpDep2012
capture rename v84 EmpRank2012
capture rename v85 HlthScore2012
capture rename v86 HlthRank2012
capture rename v87 EduScore2012
capture rename v88 EduRank2012
capture rename v89 HouseScore2012
capture rename v90 HouseRank2012
capture rename v91 GAccScore2012
capture rename v92 GAccRank2012
capture rename v93 CrimeScore2012
capture rename v94 CrimeRank2012
capture rename v95 Quintile2012
capture rename v96 Decile2012
capture rename v97 Vigintile2012
capture rename v98 Percentile2012

sort datazone
save, replace

use Dataset7_Scottish_IMD_SIMD_2016_2020, clear
merge datazone using Dataset8_Scottish_IMD_SIMD_2004_2006_2009_2012
capture drop _merge
sort datazone
save Scottish_Dataset, replace

//Merge Welsh Datasets, Wide//

use Dataset9_Welsh_IMD_WIMD_2014_2019, clear
capture rename v1 lsoa_code
capture rename v2 LSOA_Name_Eng
capture rename v3 LSOA_Name_ONS_classification
capture rename v4 Local_Authority_Name_Eng
capture rename v5 WIMD2014_r
capture rename v6 Income2014_r	
capture rename v7 Employment2014	
capture rename v8 Health2014	
capture rename v9 Education2014	
capture rename v10 Access_to_Services2014	
capture rename v11 Community_Safety2014	
capture rename v12 Physical_Environment2014	
capture rename v13 Housing2014	
capture rename v14 WIMD2014_Overall_Rank_r
capture rename v15 WIMD2014_Overall_Decile
capture rename v16 WIMD2014_Overall_Quintile
capture rename v17 WIMD2014_Overall_Quartile	
capture rename v18 WIMD_Score2019	
capture rename v19 Income_Score2019	
capture rename v20 Employment_Score2019	
capture rename v21 Health_Score2019		
capture rename v22 Education_Score2019	
capture rename v23 AccessServices_Score2019
capture rename v24 Housing_Score2019	
capture rename v25 CommunitySafety_Score2019	
capture rename v26 PhysicalEnvironment_Score2019	
capture rename v27 WIMD_Rank2019	
capture rename v28 Income_Rank2019	
capture rename v29 Employment_Rank2019	
capture rename v30 Health_Rank2019	
capture rename v31 Education_Rank2019	
capture rename v32 AccessServices_Rank2019	
capture rename v33 Housing_Rank2019	
capture rename v34 CommunitySafety_Rank2019		
capture rename v35 PhysicalEnvironment_Rank2019	
capture rename v36 WIMD_Decile2019
capture rename v37 WIMD_Quintile2019
capture rename v38 WIMD_Quartile2019
save, replace
sort lsoa_code
recast str lsoa_code
save, replace

use Dataset10_Welsh_IMD_WIMD_2005_2008_2011, clear
capture rename v17 HealthRank2005
capture rename v18 HealthScore2005
capture rename v42 IncomeRank2005
capture rename v46 AccesstoServicesRank2005
capture rename v47 HousingRank2005
capture rename v48 PhysicalEnvironmentRank2005
capture rename v68 EmploymentScore2011
capture rename v69 IncomeScore2011
capture rename v71 HealthScore2011
capture rename v75 HousingScore2011
capture rename v76 ChildIndex2011
capture rename v44 Health_Rank2005
capture rename lsoa01cd lsoa_code
save, replace
sort lsoa_code
save, replace

use Dataset9_Welsh_IMD_WIMD_2014_2019, clear
merge lsoa_code using Dataset10_Welsh_IMD_WIMD_2005_2008_2011
capture drop _merge
sort lsoa_code
save Welsh_Dataset, replace

//Merge Wide Datasets together, Long//

**find common variable**
use  English_Data, clear
findname, local(in1)
use NorthernIrish_Dataset, clear
findname, local(in2)

**common varialbe English lsoa and lgdname, change NorthernIrish_Dataset lgdname to lsoa***

use NorthernIrish_Dataset, clear
capture rename lgd2014code lsoa
sort lsoa
save, replace 

**Merge English & NorthernIrish Data, Long**
use  English_Data, clear 
merge lsoa using NorthernIrish_Dataset
capture drop _merge
sort lsoa
save IMD_UK_Final, replace

**find common variable**
use  Scottish_Dataset, clear
findname, local(in1)
use Welsh_Dataset, clear
findname, local(in2)

***common varialbe Scottish and Welsh is datazone and lsoa_code, change data_zone and lsoa_code to lsoa***

use Scottish_Dataset, clear
capture rename datazone lsoa
sort lsoa
save, replace

use Welsh_Dataset,clear
capture rename lsoa_code lsoa
sort lsoa
save, replace

***Merge Scottish and Welsh Data to IMD_UK_Final***
use IMD_UK_Final, clear
merge lsoa using Scottish_Dataset
capture drop _merge
sort lsoa
save, replace

use IMD_UK_Final, clear
merge lsoa using Welsh_Dataset
capture drop _merge
sort lsoa
save, replace

findname
use IMD_UK_Final.dta



* Create the 'country' variable based on the first character of lsoa
capture drop country
gen country = ""
replace country = "England" if substr(lsoa, 1, 1) == "E"
replace country = "Wales"   if substr(lsoa, 1, 1) == "W"
replace country = "Scotland" if substr(lsoa, 1, 1) == "S"

save, replace


codebook if country=="England"
codebook if country=="Wales"
codebook if country=="Scotland"


save, replace

capture log close




***********************MERGE THE DATASET WITH LSOA variable***************************

use IMD_UK_Final.dta,clear
sort lsoa



capture drop _merge

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\IMD_UK_Final.dta",replace


use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_merged.dta",clear

capture drop lsoa
gen lsoa=lsoa11
sort lsoa
capture drop _merge

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_merged.dta",replace

merge lsoa using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\IMD_UK_Final.dta"

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfin.dta",replace



capture log close

capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\LASSO_ENGLAND.smcl",replace

************************MACHINE LEARNING MODEL FOR PREDICTORS FOR FOOD DESERT FOR ENGLAND, WALES AND SCOTTLAND**********************
  
  
use   "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfin.dta",clear
  


**England potential predictor variables****

su imdadjusted nationalquintile1 nationaldecile2 imd_rank indexofmultipledeprivationimdsco indexofmultipledeprivationimdran indexofmultipledeprivationimddec incomescorerate incomerankwhere1ismostdeprived incomedecilewhere1ismostdeprived employmentscorerate employmentrankwhere1ismostdepriv employmentdecilewhere1ismostdepr educationskillsandtrainingscore educationskillsandtrainingrankwh educationskillsandtrainingdecile healthdeprivationanddisabilitysc healthdeprivationanddisabilityra healthdeprivationanddisabilityde crimescore crimerankwhere1ismostdeprived crimedecilewhere1ismostdeprived1 barrierstohousingandservicesscor barrierstohousingandservicesrank barrierstohousingandservicesdeci livingenvironmentscore livingenvironmentrankwhere1ismos livingenvironmentdecilewhere1ism IDACI_Score IDAC_Rank IDACI_Decile IDAOPI_Score IDAOPI_Rank IDAOPI_Decile childrenandyoungpeoplesubdomains childrenandyoungpeoplesubdomainr childrenandyoungpeoplesubdomaind adultskillssubdomainscore adultskillssubdomainrankwhere1is adultskillssubdomaindecilewhere1 geographicalbarrierssubdomainsco geographicalbarrierssubdomainran geographicalbarrierssubdomaindec widerbarrierssubdomainscore widerbarrierssubdomainrankwhere1 widerbarrierssubdomaindecilewher indoorssubdomainscore indoorssubdomainrankwhere1ismost indoorssubdomaindecilewhere1ismo outdoorssubdomainscore outdoorssubdomainrankwhere1ismos outdoorssubdomaindecilewhere1ism totalpopulationmid2012excludingp dependentchildrenaged015mid2012e populationaged1659mid2012excludi olderpopulationaged60andovermid2 workingagepopulation185964foruse

corr zsoa_pctinv imdadjusted nationalquintile1 nationaldecile2 imd_rank indexofmultipledeprivationimdsco indexofmultipledeprivationimdran indexofmultipledeprivationimddec incomescorerate incomerankwhere1ismostdeprived incomedecilewhere1ismostdeprived employmentscorerate employmentrankwhere1ismostdepriv employmentdecilewhere1ismostdepr educationskillsandtrainingscore educationskillsandtrainingrankwh educationskillsandtrainingdecile healthdeprivationanddisabilitysc healthdeprivationanddisabilityra healthdeprivationanddisabilityde crimescore crimerankwhere1ismostdeprived crimedecilewhere1ismostdeprived1 barrierstohousingandservicesscor barrierstohousingandservicesrank barrierstohousingandservicesdeci livingenvironmentscore livingenvironmentrankwhere1ismos livingenvironmentdecilewhere1ism IDACI_Score IDAC_Rank IDACI_Decile IDAOPI_Score IDAOPI_Rank IDAOPI_Decile childrenandyoungpeoplesubdomains childrenandyoungpeoplesubdomainr childrenandyoungpeoplesubdomaind adultskillssubdomainscore adultskillssubdomainrankwhere1is adultskillssubdomaindecilewhere1 geographicalbarrierssubdomainsco geographicalbarrierssubdomainran geographicalbarrierssubdomaindec widerbarrierssubdomainscore widerbarrierssubdomainrankwhere1 widerbarrierssubdomaindecilewher indoorssubdomainscore indoorssubdomainrankwhere1ismost indoorssubdomaindecilewhere1ismo outdoorssubdomainscore outdoorssubdomainrankwhere1ismos outdoorssubdomaindecilewhere1ism totalpopulationmid2012excludingp dependentchildrenaged015mid2012e populationaged1659mid2012excludi olderpopulationaged60andovermid2 workingagepopulation185964foruse


************List of domain and sub-domain scores computed for England**********************

reg zfd incomescorerate employmentscorerate educationskillsandtrainingscore healthdeprivationanddisabilitysc crimescore barrierstohousingandservicesscor livingenvironmentscore

capture drop healthdeprivationanddisabsc
capture rename healthdeprivationanddisabilitysc healthdeprivationanddisabsc

capture drop barrierstohousingandservscore
capture rename barrierstohousingandservicesscor barrierstohousingandservscore

capture drop childrenandyoungpeoplesubdom
capture rename childrenandyoungpeoplesubdomains childrenandyoungpeoplesubdom

capture drop geographicalbarrierssubdomscore
capture rename geographicalbarrierssubdomainsco geographicalbarrierssubdomscore


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


foreach x of varlist incomescorerate employmentscorerate educationskillsandtrainingscore healthdeprivationanddisabsc crimescore barrierstohousingandservscore  livingenvironmentscore {
capture drop z`x'
egen z`x'=std(`x')
	
}
	

su zincomescorerate zemploymentscorerate zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zcrimescore zbarrierstohousingandservscore zlivingenvironmentscore

corr zfd zsoa_pctinv zincomescorerate zemploymentscorerate zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zcrimescore zbarrierstohousingandservscore zlivingenvironmentscore

capture drop sample_LASSO_ENGLAND
sort nation lsoa
save, replace
splitsample if nation=="England" , generate(sample_LASSO_ENGLAND) nsplit(2) rseed(1234)
tab sample_LASSO_ENGLAND

save, replace

lasso linear zfd  zincomescorerate zemploymentscorerate zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zcrimescore zbarrierstohousingandservscore  zlivingenvironmentscore if sample_LASSO_ENGLAND==1 & nation=="England", rseed(1234)


cvplot 
graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_ENGLAND_CV.gph",replace

estimates store cvENGLAND

lassoknots, display(nonzero osr2 bic)

lassoselect id=74

cvplot

estimates store minBICENGLAND

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_ENGLAND_MINBIC.gph",replace

lasso linear zfd zincomescorerate zemploymentscorerate zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zcrimescore zbarrierstohousingandservscore zlivingenvironmentscore if sample_LASSO_ENGLAND==1 & nation=="England", rseed(1234) selection(adaptive)

cvplot

estimates store adaptiveENGLAND


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_ENGLAND_ADAPTIVE.gph",replace


lassocoef cvENGLAND minBICENGLAND adaptiveENGLAND, sort(coef, standardized) nofvlabel
 

**----------------------------------------------------------------------------
**                                 | cvENGLAND minBICENGLAND  adaptiveENGLAND 
**---------------------------------+------------------------------------------
**            zemploymentscorerate |     x           x               x        
**                zincomescorerate |     x           x               x        
**                     zcrimescore |     x           x               x        
**  zbarrierstohousingandservscore |     x           x               x        
**zeducationskillsandtrainingscore |     x           x               x        
**    zhealthdeprivationanddisabsc |     x           x               x        
**         zlivingenvironmentscore |     x           x               x        
**                           _cons |     x           x               x        
**----------------------------------------------------------------------------

				
				
 
lassogof cvENGLAND minBICENGLAND adaptiveENGLAND, over(sample_LASSO_ENGLAND) postselection


**Postselection coefficients
*-------------------------------------------------------------
**Name         sample_L~D |         MSE    R-squared        Obs
**------------------------+------------------------------------
**cvENGLAND               |
**                      1 |    .4744777       0.4198     16,422
**                      2 |    .4817145       0.4231     16,422
**------------------------+------------------------------------
**minBICENGLAND           |
**                      1 |    .4744777       0.4198     16,422
**                      2 |    .4817145       0.4231     16,422
**------------------------+------------------------------------
**adaptiveENGLAND         |
**                      1 |    .4744777       0.4198     16,422
**                      2 |    .4817145       0.4231     16,422
**-------------------------------------------------------------



save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


capture log close

capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\LASSO_SCOTTLAND.smcl",replace

****************************SCOTTLAND**************************************************

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear

su simd2020*

codebook simd2020*

su simd2020v2_income_domain_rank  simd2020_employment_domain_rank simd2020_education_domain_rank  simd2020_access_domain_rank  simd2020_crime_domain_rank    simd2020_housing_domain_rank 

foreach x of varlist simd2020v2_income_domain_rank  simd2020_employment_domain_rank simd2020_education_domain_rank  simd2020_access_domain_rank  simd2020_crime_domain_rank    simd2020_housing_domain_rank  {
	capture drop z`x'
	egen z`x'=std(`x')
	
}



capture drop sample_LASSO_SCOTLAND
sort nation lsoa
save, replace
splitsample if nation=="Scotland" , generate(sample_LASSO_SCOTLAND) nsplit(2) rseed(1234)
tab sample_LASSO_SCOTLAND

save, replace


su zsimd2020v2_income_domain_rank  zsimd2020_employment_domain_rank zsimd2020_education_domain_rank  zsimd2020_access_domain_rank  zsimd2020_crime_domain_rank    zsimd2020_housing_domain_rank 

corr zfd zsoa_pct zsimd2020v2_income_domain_rank  zsimd2020_employment_domain_rank zsimd2020_education_domain_rank  zsimd2020_access_domain_rank  zsimd2020_crime_domain_rank    zsimd2020_housing_domain_rank   
   
*******************


lasso linear zfd  zsimd2020v2_income_domain_rank  zsimd2020_employment_domain_rank zsimd2020_education_domain_rank  zsimd2020_access_domain_rank  zsimd2020_crime_domain_rank    zsimd2020_housing_domain_rank  if sample_LASSO_SCOTLAND==1 & nation=="Scotland", rseed(1234)


cvplot 
graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_SCOTLAND_CV.gph",replace

estimates store cvSCOTLAND

lassoknots, display(nonzero osr2 bic)

lassoselect id=61

cvplot

estimates store minBICSCOTLAND

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_SCOTLAND_MINBIC.gph",replace


lasso linear zsimd2020v2_income_domain_rank  zsimd2020_employment_domain_rank zsimd2020_education_domain_rank  zsimd2020_access_domain_rank  zsimd2020_crime_domain_rank    zsimd2020_housing_domain_rank  if sample_LASSO_SCOTLAND==1 & country=="Scotland", rseed(1234) selection(adaptive)

cvplot

estimates store adaptiveSCOTLAND


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_SCOTLAND_ADAPTIVE.gph",replace


lassocoef cvSCOTLAND minBICSCOTLAND adaptiveSCOTLAND, sort(coef, standardized) nofvlabel


**-------------------------------------------------------------------------------
**                                 | cvSCOTLAND minBICSCOTLAND  adaptiveSCOTLAND 
**---------------------------------+---------------------------------------------
**  zsimd2020v2_income_domain_rank |     x             x       
**    zsimd2020_access_domain_rank |     x             x       
** zsimd2020_employment_domain_rank |     x             x                x        
**     zsimd2020_crime_domain_rank |     x             x       
**   zsimd2020_housing_domain_rank |     x             x                x        
** zsimd2020_education_domain_rank |     x             x                x        
**                           _cons |     x             x                x        
**-------------------------------------------------------------------------------
**Legend:
**  b - base level
**  e - empty cell
**  o - omitted
**  x - estimated


 
lassogof cvSCOTLAND minBICSCOTLAND adaptiveSCOTLAND, over(sample_LASSO_SCOTLAND) postselection

**Postselection coefficients
**-------------------------------------------------------------
**Name         samp~TLAND |         MSE    R-squared        Obs
**------------------------+------------------------------------
**cvSCOTLAND              |
**                      1 |    .9659442       0.3761      3,488
**                      2 |    .9942468       0.3676      3,488
**------------------------+------------------------------------
**minBICSCOTLAND          |
**                      1 |    .9659442       0.3761      3,488
**                      2 |    .9942468       0.3676      3,488
**------------------------+------------------------------------
**adaptiveSCOTLAND        |
**                      1 |    .0559045       0.9442      3,488
**                      2 |    .0568349       0.9431      3,488
**-------------------------------------------------------------



save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\LASSO_WALES.smcl",replace


use  "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear


************************WALES*********************************

su Income_Score2019 Employment_Score2019 Health_Score2019 Education_Score2019 AccessServices_Score2019 Housing_Score2019 CommunitySafety_Score2019 PhysicalEnvironment_Score2019

 

foreach x of varlist Income_Score2019 Employment_Score2019 Health_Score2019 Education_Score2019 AccessServices_Score2019 Housing_Score2019 CommunitySafety_Score2019 PhysicalEnvironment_Score2019 {
	capture drop z`x'
	egen z`x'=std(`x')
	
}



capture drop sample_LASSO_WALES
sort nation lsoa
save, replace
splitsample if nation=="Wales" , generate(sample_LASSO_WALES) nsplit(2) rseed(1234)
tab sample_LASSO_WALES

save, replace


su zIncome_Score2019 zEmployment_Score2019 zHealth_Score2019 zEducation_Score2019 zAccessServices_Score2019 zHousing_Score2019 zCommunitySafety_Score2019 PhysicalEnvironment_Score2019

corr zfd zsoa_pct zIncome_Score2019 zEmployment_Score2019 zHealth_Score2019 zEducation_Score2019 zAccessServices_Score2019 zHousing_Score2019 zCommunitySafety_Score2019 PhysicalEnvironment_Score2019
   
*******************


lasso linear zfd  zIncome_Score2019 zEmployment_Score2019 zHealth_Score2019 zEducation_Score2019 zAccessServices_Score2019 zHousing_Score2019 zCommunitySafety_Score2019 zPhysicalEnvironment_Score2019 if sample_LASSO_WALES==1 & nation=="Wales", rseed(1234)


cvplot 
graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_WALES_CV.gph",replace

estimates store cvWALES

lassoknots, display(nonzero osr2 bic)

lassoselect id=64

cvplot

estimates store minBICWALES

graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_WALES_MINBIC.gph",replace


lasso linear zIncome_Score2019 zEmployment_Score2019 zHealth_Score2019 zEducation_Score2019 zAccessServices_Score2019 zHousing_Score2019 zCommunitySafety_Score2019 zPhysicalEnvironment_Score2019 if sample_LASSO_WALES==1 & country=="Wales", rseed(1234) selection(adaptive)

cvplot

estimates store adaptiveWALES


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_WALES_ADAPTIVE.gph",replace


lassocoef cvWALES minBICWALES adaptiveWALES, sort(coef, standardized) nofvlabel
 

**----------------------------------------------------------------------
**                               |  cvWALES  minBICWALES  adaptiveWALES 
**-------------------------------+--------------------------------------
**     zAccessServices_Score2019 |     x          x      
**         zEmployment_Score2019 |     x          x             x       
**             zHealth_Score2019 |     x          x             x       
**    zCommunitySafety_Score2019 |     x          x      
**zPhysicalEnvironment_Score2019 |     x          x      
**             zIncome_Score2019 |     x          x      
**            zHousing_Score2019 |     x          x      
**          zEducation_Score2019 |                              x       
**                         _cons |     x          x             x       
**----------------------------------------------------------------------

 
lassogof cvWALES minBICWALES adaptiveWALES, over(sample_LASSO_WALES) postselection


**Postselection coefficients
**-------------------------------------------------------------
**Name         sample_L~S |         MSE    R-squared        Obs
**------------------------+------------------------------------
**cvWALES                 |
**                      1 |    .5824592       0.6332        955
**                      2 |    .6391186       0.6136        954
**------------------------+------------------------------------
**minBICWALES             |
**                      1 |    .5824592       0.6332        955
**                      2 |    .6391186       0.6136        954
**------------------------+------------------------------------
**adaptiveWALES           |
**                      1 |    .0680789       0.9349        955
**                      2 |    .0756418       0.9205        954
**-------------------------------------------------------------
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace




capture log close 

capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3B.smcl",replace


use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear


********************************HIERARCHICAL ADAPTIVE LASSO MODELS WITH INTERACTIONS******************************

su zincomescorerate zemploymentscorerate zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zcrimescore zbarrierstohousingandservscore zlivingenvironmentscore



capture drop lncrimescore
gen lncrimescore=ln(crimescore+10) if nation=="England" 

histogram lncrimescore

capture drop zlncrimescore
egen zlncrimescore=std(lncrimescore) if nation=="England"



capture drop zinc_sc
gen zinc_sc=zincomescorerate

capture drop zemp_sc
gen zemp_sc=zemploymentscorerate

capture drop zed_sc
gen zed_sc=zeducationskillsandtrainingscore

capture drop zhlth_sc
gen zhlth_sc=zhealthdeprivationanddisabsc

capture drop zcrm_sc
gen zcrm_sc=zlncrimescore

capture drop zhsserv_sc
gen zhsserv_sc=zbarrierstohousingandservscore

capture drop zlivenv_sc
gen zlivenv_sc=zlivingenvironmentscore

global mainvars zinc_sc zemp_sc zed_sc zhlth_sc zcrm_sc  zhsserv_sc zlivenv_sc

capture drop i_*

foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}


lasso linear zfd $mainvars i_* if sample_LASSO_ENGLAND==1 & nation=="England", rseed(1234) selection(adaptive)



**Lasso linear model                         No. of obs         =     16,422
**                                           No. of covariates  =         28
**Selection: Adaptive                        No. of lasso steps =          2
**
**Final adaptive step results
**--------------------------------------------------------------------------
**         |                                No. of      Out-of-      CV mean
**         |                               nonzero       sample   prediction
**      ID |     Description      lambda     coef.    R-squared        error
**---------+----------------------------------------------------------------
**      90 |    first lambda    3.452653         0       0.0007     .8171704
**     188 |   lambda before    .0003789        28       0.5050     .4048221
**   * 189 | selected lambda    .0003453        28       0.5050     .4048013
**--------------------------------------------------------------------------


lassoselect id=189


cvplot

estimates store adaptiveENGLANDINT


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_ENGLAND_ADAPTIVE_INTERACTIONS.gph",replace


lassocoef, sort(coef, standardized)



**--------------------------------------------
**                        | adaptiveENGLANDINT
**------------------------+-------------------
**                zemp_sc |         x         
**     i_zhlth_sc_zinc_sc |         x         
**                zinc_sc |         x         
*                zcrm_sc |         x         
**   i_zhsserv_sc_zinc_sc |         x         
**      i_zed_sc_zhlth_sc |         x         
**      i_zemp_sc_zinc_sc |         x         
**                 zed_sc |         x         
**             zhsserv_sc |         x         
**      i_zcrm_sc_zinc_sc |         x         
**   i_zemp_sc_zhsserv_sc |         x         
**       i_zed_sc_zemp_sc |         x         
**   i_zcrm_sc_zlivenv_sc |         x         
**  i_zhlth_sc_zhsserv_sc |         x         
**      i_zcrm_sc_zemp_sc |         x         
**i_zhsserv_sc_zlivenv_sc |         x         
**     i_zemp_sc_zhlth_sc |         x         
**   i_zinc_sc_zlivenv_sc |         x         
**               zhlth_sc |         x         
**   i_zemp_sc_zlivenv_sc |         x         
**    i_zed_sc_zlivenv_sc |         x         
**       i_zed_sc_zinc_sc |         x         
**   i_zcrm_sc_zhsserv_sc |         x         
**             zlivenv_sc |         x         
**     i_zcrm_sc_zhlth_sc |         x         
**  i_zhlth_sc_zlivenv_sc |         x         
**    i_zed_sc_zhsserv_sc |         x         
**       i_zcrm_sc_zed_sc |         x         
**                  _cons |         x         
**--------------------------------------------
**Legend:
**  b - base level
**  e - empty cell
**  o - omitted
**  x - estimated



lassogof adaptiveENGLANDINT, over(sample_LASSO_ENGLAND) postselection


**Postselection coefficients
**-------------------------------------------------------------
**Name         sample_L~D |         MSE    R-squared        Obs
**------------------------+------------------------------------
**adaptiveENGLANDINT      |
**                      1 |     .402898       0.5073     16,422
**                      2 |     .404754       0.5153     16,422
**-------------------------------------------------------------



reg zfd zemp_sc i_zhlth_sc_zinc_sc zinc_sc zcrm_sc i_zhsserv_sc_zinc_sc i_zed_sc_zhlth_sc i_zemp_sc_zinc_sc zed_sc zhsserv_sc i_zcrm_sc_zinc_sc i_zemp_sc_zhsserv_sc i_zed_sc_zemp_sc  i_zcrm_sc_zlivenv_sc i_zhlth_sc_zhsserv_sc i_zcrm_sc_zemp_sc i_zhsserv_sc_zlivenv_sc i_zemp_sc_zhlth_sc i_zinc_sc_zlivenv_sc zhlth_sc i_zemp_sc_zlivenv_sc i_zed_sc_zlivenv_sc i_zed_sc_zinc_sc i_zcrm_sc_zhsserv_sc zlivenv_sc i_zcrm_sc_zhlth_sc i_zhlth_sc_zlivenv_sc i_zed_sc_zhsserv_sc i_zcrm_sc_zed_sc if nation=="England"


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear

****************************************SCOTLAND**********************

su zsimd2020v2_income_domain_rank  zsimd2020_employment_domain_rank zsimd2020_education_domain_rank  zsimd2020_access_domain_rank  zsimd2020_crime_domain_rank  zsimd2020_housing_domain_rank if nation=="Scotland"

capture drop zsinc_rank
gen zsinc_rank=zsimd2020v2_income_domain_rank

capture drop zsempl_rank
gen zsempl_rank=zsimd2020_employment_domain_rank

capture drop zsed_rank
gen zsed_rank=zsimd2020_education_domain_rank

capture drop zsacc_rank
gen zsacc_rank=zsimd2020_access_domain_rank

capture drop zscrime_rank
gen zscrime_rank=zsimd2020_crime_domain_rank

capture drop zshouse_rank
gen zshouse_rank=zsimd2020_housing_domain_rank

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace



global mainvars zsinc_rank zsempl_rank zsed_rank zsacc_rank zscrime_rank  zshouse_rank 


capture drop i_*
foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}


lasso linear zfd $mainvars i_* if sample_LASSO_SCOTLAND==1 & nation=="Scotland", rseed(1234) selection(adaptive)


**Lasso linear model                         No. of obs         =      3,488
**                                           No. of covariates  =         21
**Selection: Adaptive                        No. of lasso steps =          2
**
**Final adaptive step results
**--------------------------------------------------------------------------
**         |                                No. of      Out-of-      CV mean
**         |                               nonzero       sample   prediction
**      ID |     Description      lambda     coef.    R-squared        error
**---------+----------------------------------------------------------------
**      81 |    first lambda     2.50637         0       0.0012     1.546238
**     174 |   lambda before     .000438        20       0.4272     .8868124
**   * 175 | selected lambda    .0003991        20       0.4272     .8868122
**     176 |    lambda after    .0003636        20       0.4272     .8868128
**     177 |     last lambda    .0003313        20       0.4272     .8868139
**--------------------------------------------------------------------------


lassoselect id=175


cvplot

estimates store adaptiveSCOTLANDINT


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_SCOTLAND_ADAPTIVE_INTERACTIONS.gph",replace


lassocoef, sort(coef, standardized)


**-------------------------------------------------
**                            | adaptiveSCOTLANDINT
**----------------------------+--------------------
**                 zsinc_rank |          x         
**                 zsacc_rank |          x         
**                zsempl_rank |          x         
**  i_zsacc_rank_zscrime_rank |          x         
**     i_zsed_rank_zsinc_rank |          x         
**               zscrime_rank |          x         
** i_zscrime_rank_zsempl_rank |          x         
**    i_zsed_rank_zsempl_rank |          x         
**               zshouse_rank |          x         
**   i_zsacc_rank_zsempl_rank |          x         
**  i_zscrime_rank_zsinc_rank |          x         
**   i_zsed_rank_zshouse_rank |          x         
**   i_zsempl_rank_zsinc_rank |          x         
**  i_zsacc_rank_zshouse_rank |          x         
**   i_zscrime_rank_zsed_rank |          x         
** i_zscrime_rank_zshouse_rank |          x         
**    i_zsacc_rank_zsinc_rank |          x         
**     i_zsacc_rank_zsed_rank |          x         
**                  zsed_rank |          x         
**                      _cons |          x         
**  i_zshouse_rank_zsinc_rank |          x         
**-------------------------------------------------
**Legend:
**  b - base level
**  e - empty cell
**  o - omitted
**  x - estimated



lassogof adaptiveSCOTLANDINT, over(sample_LASSO_SCOTLAND) postselection


**Postselection coefficients
**-------------------------------------------------------------
**Name         samp~TLAND |         MSE    R-squared        Obs
**------------------------+------------------------------------
**adaptiveSCOTLANDINT     |
**                      1 |    .8775342       0.4332      3,488
**                      2 |    .8981556       0.4287      3,488
**-------------------------------------------------------------


reg zfd zsinc_rank zsacc_rank zsempl_rank i_zsacc_rank_zscrime_rank i_zsed_rank_zsinc_rank zscrime_rank i_zscrime_rank_zsempl_rank i_zsed_rank_zsempl_rank zshouse_rank i_zsacc_rank_zsempl_rank i_zscrime_rank_zsinc_rank i_zsed_rank_zshouse_rank i_zsempl_rank_zsinc_rank i_zsacc_rank_zshouse_rank i_zscrime_rank_zsed_rank i_zscrime_rank_zshouse_rank i_zsacc_rank_zsinc_rank i_zsacc_rank_zsed_rank zsed_rank i_zshouse_rank_zsinc_rank if nation=="Scotland"


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace



*********************WALES********************************

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


su zIncome_Score2019 zEmployment_Score2019 zHealth_Score2019 zEducation_Score2019 zAccessServices_Score2019 zHousing_Score2019 zCommunitySafety_Score2019 PhysicalEnvironment_Score2019


capture drop zwinc_sc
gen zwinc_sc=zIncome_Score2019

capture drop zwemp_sc
gen zwemp_sc=zEmployment_Score2019

capture drop zwhlth_sc
gen zwhlth_sc=zHealth_Score2019

capture drop zwed_sc
gen zwed_sc=zEducation_Score2019

capture drop zwacc_sc
gen zwacc_sc=zAccessServices_Score2019

capture drop zwhouse_sc
gen zwhouse_sc=zHousing_Score2019

capture drop zwcrime_sc
gen zwcrime_sc=zCommunitySafety_Score2019

capture drop zwenv_sc
gen zwenv_sc=PhysicalEnvironment_Score2019

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace



global mainvars zwinc_sc  zwemp_sc zwhlth_sc zwed_sc zwacc_sc zwhouse_sc zwcrime_sc zwenv_sc


capture drop i_*
foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}


lasso linear zfd $mainvars i_* if sample_LASSO_WALES==1 & nation=="Wales", rseed(1234) selection(adaptive)

**Final adaptive step results
**--------------------------------------------------------------------------
**         |                                No. of      Out-of-      CV mean
**         |                               nonzero       sample   prediction
**      ID |     Description      lambda     coef.    R-squared        error
**---------+----------------------------------------------------------------
**      89 |    first lambda    25.04128         0       0.0022     1.584617
**     156 |   lambda before     .049157        16       0.6644     .5329807
**   * 157 | selected lambda    .0447901        17       0.6645     .5328701
**     158 |    lambda after     .040811        17       0.6644     .5329529
**     188 |     last lambda    .0025041        30       0.6606     .5390539
**--------------------------------------------------------------------------

lassoselect id=157


cvplot

estimates store adaptiveWALESINT


graph save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LASSO_WALES_ADAPTIVE_INTERACTIONS.gph",replace


lassocoef, sort(coef, standardized)


**------------------------------------------
**                        | adaptiveWALESINT
**------------------------+-----------------
**               zwemp_sc |        x        
**               zwacc_sc |        x        
**    i_zwemp_sc_zwinc_sc |        x        
**              zwhlth_sc |        x        
**   i_zwemp_sc_zwhlth_sc |        x        
**               zwinc_sc |        x        
**             zwcrime_sc |        x        
**               zwenv_sc |        x        
** i_zwcrime_sc_zwhlth_sc |        x        
**   i_zwacc_sc_zwhlth_sc |        x        
**  i_zwacc_sc_zwcrime_sc |        x        
**             zwhouse_sc |        x        
**  i_zwenv_sc_zwhouse_sc |        x        
**i_zwcrime_sc_zwhouse_sc |        x        
**  i_zwhouse_sc_zwinc_sc |        x        
**    i_zwed_sc_zwhlth_sc |        x        
**  i_zwcrime_sc_zwenv_sc |        x        
**                  _cons |        x        
**------------------------------------------



lassogof adaptiveWALESINT, over(sample_LASSO_WALES) postselection

**Postselection coefficients
**-------------------------------------------------------------
**Name         sample_L~S |         MSE    R-squared        Obs
**------------------------+------------------------------------
**adaptiveWALESINT        |
**                      1 |     .511293       0.6780        955
**                      2 |    .5873504       0.6449        954
**-------------------------------------------------------------



reg zfd zwemp_sc zwacc_sc i_zwemp_sc_zwinc_sc zwhlth_sc i_zwemp_sc_zwhlth_sc zwinc_sc zwcrime_sc zwenv_sc i_zwcrime_sc_zwhlth_sc i_zwacc_sc_zwhlth_sc i_zwacc_sc_zwcrime_sc zwhouse_sc i_zwenv_sc_zwhouse_sc i_zwcrime_sc_zwhouse_sc i_zwhouse_sc_zwinc_sc i_zwed_sc_zwhlth_sc i_zwcrime_sc_zwenv_sc if nation=="Wales"


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace

capture log close


***************************************************************************************************************************

capture log close 

capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3A.smcl",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear


capture drop lncrimescore
gen lncrimescore=ln(crimescore+10) if nation=="England" 

histogram lncrimescore

capture drop zlncrimescore
egen zlncrimescore=std(lncrimescore)


reg zfd zemploymentscorerate zlncrimescore  zbarrierstohousingandservscore zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zlivingenvironmentscore if nation=="England" 

reg zfd  zsimd2020_employment_domain_rank zsimd2020_housing_domain_rank zsimd2020_education_domain_rank   if nation=="Scotland"

reg zfd  zEmployment_Score2019   zHealth_Score2019 zEducation_Score2019 if nation=="Wales"


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3A_MIXED.smcl",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear

capture drop la_codenum
encode laname, gen(la_codenum)

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace

mixed zfd zemploymentscorerate zlncrimescore  zbarrierstohousingandservscore zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zlivingenvironmentscore if nation=="England" || la_codenum:

mixed zfd  zsimd2020_employment_domain_rank zsimd2020_housing_domain_rank zsimd2020_education_domain_rank   if nation=="Scotland" || la_codenum:

mixed zfd  zEmployment_Score2019   zHealth_Score2019 zEducation_Score2019 if nation=="Wales" || la_codenum:

capture log close




capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3B_MIXED.smcl",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear


global mainvars zinc_sc zemp_sc zed_sc zhlth_sc zcrm_sc  zhsserv_sc zlivenv_sc

capture drop i_*

foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}



mixed zfd zemp_sc i_zhlth_sc_zinc_sc zinc_sc zcrm_sc i_zhsserv_sc_zinc_sc i_zed_sc_zhlth_sc i_zemp_sc_zinc_sc zed_sc zhsserv_sc i_zcrm_sc_zinc_sc i_zemp_sc_zhsserv_sc i_zed_sc_zemp_sc  i_zcrm_sc_zlivenv_sc i_zhlth_sc_zhsserv_sc i_zcrm_sc_zemp_sc i_zhsserv_sc_zlivenv_sc i_zemp_sc_zhlth_sc i_zinc_sc_zlivenv_sc zhlth_sc i_zemp_sc_zlivenv_sc i_zed_sc_zlivenv_sc i_zed_sc_zinc_sc i_zcrm_sc_zhsserv_sc zlivenv_sc i_zcrm_sc_zhlth_sc i_zhlth_sc_zlivenv_sc i_zed_sc_zhsserv_sc i_zcrm_sc_zed_sc if nation=="England" || la_codenum:

global mainvars zsinc_rank zsempl_rank zsed_rank zsacc_rank zscrime_rank  zshouse_rank 


capture drop i_*
foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}



mixed zfd  zsinc_rank zsacc_rank zsempl_rank i_zsacc_rank_zscrime_rank i_zsed_rank_zsinc_rank zscrime_rank i_zscrime_rank_zsempl_rank i_zsed_rank_zsempl_rank zshouse_rank i_zsacc_rank_zsempl_rank i_zscrime_rank_zsinc_rank i_zsed_rank_zshouse_rank i_zsempl_rank_zsinc_rank i_zsacc_rank_zshouse_rank i_zscrime_rank_zsed_rank i_zscrime_rank_zshouse_rank i_zsacc_rank_zsinc_rank i_zsacc_rank_zsed_rank zsed_rank i_zshouse_rank_zsinc_rank   if nation=="Scotland" || la_codenum:



global mainvars zwinc_sc  zwemp_sc zwhlth_sc zwed_sc zwacc_sc zwhouse_sc zwcrime_sc zwenv_sc


capture drop i_*
foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}




mixed zfd  zwemp_sc zwacc_sc i_zwemp_sc_zwinc_sc zwhlth_sc i_zwemp_sc_zwhlth_sc zwinc_sc zwcrime_sc zwenv_sc i_zwcrime_sc_zwhlth_sc i_zwacc_sc_zwhlth_sc i_zwacc_sc_zwcrime_sc zwhouse_sc i_zwenv_sc_zwhouse_sc i_zwcrime_sc_zwhouse_sc i_zwhouse_sc_zwinc_sc i_zwed_sc_zwhlth_sc i_zwcrime_sc_zwenv_sc if nation=="Wales" || la_codenum:

save  "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


capture log close

*****************************MIXED MODELS WITH SAFEGARDED DATA LSOA-level groupings: England and Wales*****************************
capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3A_MIXED.smcl",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear

capture drop la_codenum
encode laname, gen(la_codenum)

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace

mixed zfd zemploymentscorerate zlncrimescore  zbarrierstohousingandservscore zeducationskillsandtrainingscore zhealthdeprivationanddisabsc zlivingenvironmentscore if nation=="England" || la_codenum:

mixed zfd  zsimd2020_employment_domain_rank zsimd2020_housing_domain_rank zsimd2020_education_domain_rank   if nation=="Scotland" || la_codenum:

mixed zfd  zEmployment_Score2019   zHealth_Score2019 zEducation_Score2019 if nation=="Wales" || la_codenum:

capture log close



capture log close
capture log using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\OUTPUT\TABLE3B_MIXED_SAFEGUARDED.smcl",replace

********************ENGLAND AND WALES, MIXED-MODELS WITH SAFEGARDED DATA GROUPING*****************

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear
capture drop _merge
sort lsoa
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",replace


clear
import delimited "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LSOA_DZ_SDZ_Lookup.csv"
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LSOA_DZ_SDZ_Lookup.dta",replace
capture rename geography_code lsoa
sort lsoa
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LSOA_DZ_SDZ_Lookup.dta",replace

clear
use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalized.dta",clear
merge lsoa using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\LSOA_DZ_SDZ_Lookup.dta"
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalizedfin.dta",replace



global mainvars zinc_sc zemp_sc zed_sc zhlth_sc zcrm_sc  zhsserv_sc zlivenv_sc

capture drop i_*

foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}



mixed zfd zemp_sc i_zhlth_sc_zinc_sc zinc_sc zcrm_sc i_zhsserv_sc_zinc_sc i_zed_sc_zhlth_sc i_zemp_sc_zinc_sc zed_sc zhsserv_sc i_zcrm_sc_zinc_sc i_zemp_sc_zhsserv_sc i_zed_sc_zemp_sc  i_zcrm_sc_zlivenv_sc i_zhlth_sc_zhsserv_sc i_zcrm_sc_zemp_sc i_zhsserv_sc_zlivenv_sc i_zemp_sc_zhlth_sc i_zinc_sc_zlivenv_sc zhlth_sc i_zemp_sc_zlivenv_sc i_zed_sc_zlivenv_sc i_zed_sc_zinc_sc i_zcrm_sc_zhsserv_sc zlivenv_sc i_zcrm_sc_zhlth_sc i_zhlth_sc_zlivenv_sc i_zed_sc_zhsserv_sc i_zcrm_sc_zed_sc i.supergroup  if nation=="England" || la_codenum:






global mainvars zwinc_sc  zwemp_sc zwhlth_sc zwed_sc zwacc_sc zwhouse_sc zwcrime_sc zwenv_sc


capture drop i_*
foreach v1 in $mainvars {
    foreach v2 in $mainvars {
        if "`v1'" < "`v2'" {
            gen i_`v1'_`v2' = `v1' * `v2'
        }
    }
}




mixed zfd  zwemp_sc zwacc_sc i_zwemp_sc_zwinc_sc zwhlth_sc i_zwemp_sc_zwhlth_sc zwinc_sc zwcrime_sc zwenv_sc i_zwcrime_sc_zwhlth_sc i_zwacc_sc_zwhlth_sc i_zwacc_sc_zwcrime_sc zwhouse_sc i_zwenv_sc_zwhouse_sc i_zwcrime_sc_zwhouse_sc i_zwhouse_sc_zwinc_sc i_zwed_sc_zwhlth_sc i_zwcrime_sc_zwenv_sc i.supergroup if nation=="Wales" || la_codenum:

save  "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalizedfin.dta",replace


capture log close
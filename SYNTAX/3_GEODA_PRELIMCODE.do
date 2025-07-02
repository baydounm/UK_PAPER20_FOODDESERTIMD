
**GEODA PRELIMINARY CODE**

clear
use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\UK_fooddesertdatafinalizedIMD_mergedfinalizedfin.dta" 

capture drop ID_GEODA
encode lsoa, gen(ID_GEODA)

su ID_GEODA

sort ID_GEODA
capture drop ID_FINAL
gen ID_FINAL=_n

su ID_FINAL

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20.dta",replace

export delimited using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20.csv",replace

**Make a copy of the dbf file and call it infuse_lsoa_lyr_2011_clipped FD_IMD.dbf**
clear
import dbase using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped FD_IMD.dbf"

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped FD_IMD.dta",replace
sort geo_code
save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped FD_IMD.dta",replace

**Get a smaller version of the FD_IMD dataset that has only the LSOA code, zFD and zIMD variables and merge with above dataset by LSOA code**

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20.dta",clear

keep lsoa zfd fd zsoa* soa* ID_GEODA ID_FINAL longitude latitude 

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20_SMALL.dta",replace
sort lsoa
capture drop geo_code
gen geo_code=lsoa
sort geo_code

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20_SMALL.dta",replace

use "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\GEODA_DATA_PAPER20_SMALL.dta",clear
merge geo_code using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped FD_IMD.dta"

save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped FD_IMDfin.dta",replace

tab _merge

capture drop _merge


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped_FD_IMDfin.dta",replace

capture rename zsoa_pctinv zIMD 
capture rename zsoa_pctinv_cluster IMDCLUSTER
capture rename zfd zFD

keep lsoa zIMD zFD geo_code geo_label longitude latitude ID_FINAL
order geo_code geo_label zIMD zFD
format geo_code %200s
format lsoa %200s
format geo_label %200s

capture drop geo_code_short
gen geo_code_short = substr(geo_code, 1, 254)
capture drop geo_code
capture drop lsoa
rename geo_code_short geo_code



save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped_FD_IMDfin2.dta",replace

export dbase using "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\DATA\GEODA_DATA\infuse_lsoa_lyr_2011_clipped_FD_IMDfin2.dbf",replace

**Rename all shapefile variable similarly: infuse_lsoa_lyr_2011_clipped_FD_IMDfin2**



**Follow th steps below:



Step 1: Load Your Shapefile
Open GeoDa.

Click File > Open and load your LSOA-level shapefile (e.g., infuse_lsoa_lyr_2011_clipped_FD_IMDfin2.shp).

Ensure that the attribute table includes your variables of interest:

FD: Food Desert score

IMD: Index of Multiple Deprivation

🧩 Step 2: Create a Spatial Weights Matrix
You must first define neighborhood structure via spatial weights:

Option A: Contiguity-Based Weights
Go to Tools > Weights > Create.

Name your weights file (e.g., lsoa_weights.gal).

Choose:

Contiguity-based (e.g., Queen or Rook):

Queen: common borders or vertices

Rook: common borders only

Click Create.

Option B: Distance-Based Weights (if LSOAs are not contiguous)
Choose Distance-based instead.

Use k-nearest neighbors (e.g., k = 4 or 6) or threshold distance.

If using knn, tick standardized to row-standardize weights.

For UK LSOAs, contiguity-based weights are typically appropriate since they are well-defined polygons.

**Used knn=4, distance-based.  

Step 3: Run Univariate Moran's I
Once the weights file is created:

Go to Space > Univariate Moran's I.

Select:

Variable: FD (run first)

Weights file: lsoa_weights.gal

Click Run.

GeoDa will:

Display a Moran scatterplot

Report Moran's I, expected I, z-score, and p-value

Repeat for IMD by selecting it in the Variable dropdown.


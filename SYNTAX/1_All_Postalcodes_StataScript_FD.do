
**E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\UKPOSTALCODES_LOCATIONS\multi_csv

**WORK ON EACH FILE IN THE MULT_CSV FOLDER AND SELECT ONLY THE 2001, 2011 AND 2021 CENSUS POSTCODES ALONG WITH THEIR LAT AND LONG VAR AND THE ACTUAL NAME OF THE POSTCODE**



**USE THE 2011 POSTCODES WHICH MATCH WITH THE 2018 DATA OF FOOD DESERTS. NO NEED TO REMOVE ONE CODE AT THE END. THOSE ARE EXACT MATCHES**

cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FOOD_DESERT_DATA\PCD_LSOA_LONGLAT"

**POSTCODES: APPEND THE 2011 CENSUS POSTCODES FROM THE _SMALL DATASETS ALONG WITH LONGITUDE AND LATITUDE. COLLAPSE TO OBTAIN A SINGLE POSTCODE AND AVERAGE CENTROID LOCATIONS**
clear
local filelist "ONSPD_NOV_2024_UK_SMALL_FILES_LIST.txt"  // Path to the file list

// Ensure any previous file handle named "listfile" is closed
capture file close listfile

// Open the file list
file open listfile using `filelist', read
file read listfile line

// Initialize a counter
local counter = 1

// Loop over the lines in the file list
while r(eof) == 0 {
    // Extract the CSV file name from the line
    local csvfile = subinstr("`line'", "PCD_LSOA_LONGLAT/", "", .)  // Adjust if needed
    local dtafile = subinstr("`csvfile'", ".csv", ".dta", .)

    // Import the CSV file
    display "Importing file (`counter'): `csvfile'"
    import delimited `csvfile', clear

    // Rename variables (adjust as needed)
    capture rename v* longitude
    capture rename lat latitude

    // Save the file as .dta
    display "Saving file (`counter'): `dtafile'"
    save `dtafile', replace

    // Read the next line in the file list
    file read listfile line

    // Increment the counter
    local counter = `counter' + 1
}

// Close the file list
file close listfile


*****************APPEND PART OF THE SCRIPT*****************
clear
local filelist "ONSPD_NOV_2024_UK_SMALL_FILES_LIST.txt"  // Path to the file list

// Ensure any previous file handle named "listfile" is closed
capture file close listfile

// Open the file list
file open listfile using `filelist', read
file read listfile line

// Initialize a counter
local counter = 1

// Loop over the lines in the file list
while r(eof) == 0 {
    // Extract the DTA file name from the line
    local csvfile = subinstr("`line'", "PCD_LSOA_LONGLAT/", "", .)  // Adjust if needed
    local dtafile = subinstr("`csvfile'", ".csv", ".dta", .)

    // Check if it is the first file
    if `counter' == 1 {
        // Load the first file
        display "Loading file (`counter'): `dtafile'"
        use `dtafile', clear
    }
    else {
        // Append subsequent files
        display "Appending file (`counter'): `dtafile'"
        append using `dtafile', force
    }

    // Read the next line in the file list
    file read listfile line

    // Increment the counter
    local counter = `counter' + 1
}

// Close the file list
file close listfile

// Save the appended dataset
save combined_ONSPD.dta, replace

display "All files have been appended and saved as combined_ONSPD.dta."


capture drop longitude_final
gen longitude_final=longitude if longitude~=.
su longitude_final

capture rename longitude_final longitude

capture drop pcd2 oac11 msoa21 oa21

order pcd lsoa11 lsoa21 lsoa01 latitude longitude

save combined_ONSPD.dta, replace


save "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FOOD_DESERT_DATA\combined_ONSPD.dta",replace



**FINAL FILE: "POSTCODES_2011CENSUS_UK.dta"

**DATA MANAGEMENT FOR ENGLAND/WALES AND SCOTTLAND: MERGE WITH LOSA OF 2011 TO OBTAIN LONGITUDE AND LATITUDE**


*****************ENGLAND AND WALES*******************
clear

cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FOOD_DESERT_DATA\

import delimited efdi_england.csv,clear

save "England_fooddesertdata.dta",replace

use "England_fooddesertdata.dta",clear
capture drop Country
gen Country="England and Wales"

save "England_fooddesertdata.dta",replace

capture drop lsoa11
gen lsoa11=lsoaordz

capture drop lsoaordz

sort lsoa11

save "England_fooddesertdatafin.dta",replace

******************SCOTTLAND*************************
clear

cd "E:\16GBBACKUPUSB\BACKUP_USB_SEPTEMBER2014\May Baydoun_folder\UK_BIOBANK_PROJECT\UKB_PAPER20_ADI_FOODESERT\FOOD_DESERT_DATA\

import delimited efdi_scotland.csv,clear

save "Scottland_fooddesertdata.dta",replace

use "Scottland_fooddesertdata.dta",clear
capture drop Country
gen Country="Scottland"

save "Scottland_fooddesertdata.dta",replace

capture drop lsoa11
gen lsoa11=lsoaordz


capture drop lsoaordz

sort lsoa11

save "Scottland_fooddesertdatafin.dta",replace

*******************APPEND ENGLAND/WALES ON SCOTTLAND**********************************

use "England_fooddesertdatafin.dta",clear
append using "Scottland_fooddesertdata.dta"
save "UK_fooddesertdatafin.dta",replace
sort lsoa11
replace lsoa11=lsoaordz if lsoa11==""

save "UK_fooddesertdatafin.dta",replace




use "combined_ONSPD.dta",clear
sort lsoa11
save "combined_ONSPD.dta",replace


use "UK_fooddesertdatafin.dta",clear
sort lsoa11
save "UK_fooddesertdatafin.dta",replace

merge lsoa11 using "combined_ONSPD.dta"

tab _merge


******************COLLAPSE THE COMBINED PCD FILE BY LSOA11***************************

use combined_ONSPD.dta,clear

collapse (mean) longitude latitude, by(lsoa11)
save "combined_ONSPD_collapse_LSOA11",replace
sort lsoa11
save "combined_ONSPD_collapse_LSOA11",replace

********************MERGE WITH ENGLAND/SCOTTLAND COMBINED BY LSOA11***************************

use "UK_fooddesertdatafin.dta",clear
sort lsoa11
save "UK_fooddesertdatafin.dta",replace

merge lsoa11 using "combined_ONSPD_collapse_LSOA11"
save "UK_fooddesertdatafinalized.dta",replace
tab _merge

capture drop fd
gen fd=score


save "UK_fooddesertdatafinalized.dta",replace




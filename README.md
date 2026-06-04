Data and Code for: Linking plant sex to belowground communities and biogeochemical properties: implications for soil food web structure
Authors
Alexandra Bradley Martin, Christopher C. Oxley, Kerri M. Crawford, Dajanae A. Williams, Tingfa Dong, Tingfu Zhang, Junyan Liu, Mingsen Qin, Ashley D. Keiser, James H. Leebens-Mack, Tom E.X. Miller, Matthew A. McCary
Repository Overview
This repository contains all raw data and R code necessary to reproduce the analyses and figures presented in the manuscript.
Software Requirements
Analyses were performed using:
R version 4.3.3
RStudio version 2023.06.2 
Required R packages are listed in the code file.

Data Description
All_files_Metadata.xlsx
Description: Metadata for all the files. Each tab in this file matches with all files’ names.
Variables:
Log_ID: Unique sample identifier.
Site Study location. Levels:
Wi, AWi or Wildflower = LBJWC 
H, BH or Huntsville = PERL
TC, CTC or Turkey Creek = BTNP
Sex Plant sex. Levels:
M = Male
F= Female
CollectionDate Date of sample collection.
Ilex_Soil_Arthropods.csv
Description: Raw arthropod count data.
Berlese_Soil_Weights.csv
Description: All samples’ weight
Ilex_traits.csv
Description: Ilex vomitoria’s morphological traits and soil moisture measures.
Fungal_OTUs.csv
Description: Raw fungal sequencing data.
Rows: Samples.
Columns: Fungal OTU or “taxon”.
Cell values: Sequence read counts.
Bacterial_OTUs.csv
Description: Raw bacterial sequencing data.
Rows: Samples.
Columns: Bacterial OTU or “taxon”.
Cell values: Sequence read counts.
Fungal_taxon.csv
Description: Taxonomy information for each fungal OTU.
Bacterial_taxon.csv
Description: Taxonomy information for each bacterial OTU.
Microbial_Biomass_Data.csv
Description: soil microbial biomass and nutrient measurements.
soil_biogeochem_properties.csv
Other measured soil variables.
Soil_organic_matter_fractions.csv
Description: soil carbon and nitrogen fraction measurements.

Analysis Workflow
Open code file: Ilex_v_code.R
Script should be run in the following order:
Load all packages needed in R
Imports raw datasets.
Performs data cleaning.
2.	Diversity analyses
Calculate richness and diversity metrics.
Generate all figures used in the manuscript.

3.	Community analyses
Conducts Bray-Curtis analyses.
Performs PCoA analyses.
Performs PERMANOVA analyses.
Generate all figures used in the manuscript.

4.	Structural Equation modeling analyses
Fits statistical models.

5.	Supplementary Data
Generate all figures and tables used in the manuscript.

Reproducing the Manuscript
1. Download or clone the repository.
2. Open RStudio.
3. Set the repository folder as the working directory.
4. Install the required R packages listed at the beginning of the script.
5. Run the script "Ilex_v_code.R".

The script:
- Imports all raw data files from the data/ folder.
- Performs all statistical analyses reported in the manuscript.
- Generates all figures used in the manuscript.
- Saves figures to the figures/ folder.
Notes
All raw data are provided exactly as used in the analyses.
Questions regarding the repository can be directed to the corresponding author.

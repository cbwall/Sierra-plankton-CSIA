## Zooplankton reliance on autochthonous and allochthonous resources across mountain lakes revealed through essential amino acid  δ<sup>13</sup>C analysis.  
CB Wall<sup>1,2Δ</sup>, AC Besser<sup>3,4Δ</sup>, CC Symons<sup>5</sup>, S Newsome<sup>4</sup>, JB Shurin<sup>1</sup>. (2026) _Limnology and Oceanography_.  

<sup>1</sup>University of California, San Diego  
<sup>2</sup>University of Hawai'i at Mānoa  
<sup>3</sup>Arizona State University  
<sup>4</sup>University of New Mexico  
<sup>5</sup>University of California, Irvine  
Δ _authors contributed equally and share first authorship_
  
<p align="center">
  <img align="center" src="https://github.com/cbwall/Sierra-plankton-CSIA/blob/main/output//photos/Fig1.sitemap_2.png" width="65%" height="60%">
</p>

*Overview*:  We used CSIA-AA to ask whether zooplankton in high-elevation alpine Lakes of the Sierra Nevada Mountains obtain their essential amino acids (AA<sub>ESS</sub>) from in-water producers (particulate organic matter dominated by microalgae and benthic algae) or from terrestrial C3 plants (sedges, pines, and broadleaf deciduous trees). We sampled 7 Lakes and 1 pond across an elevation gradient (2500-3200m) and measured environmental conditions (temperature, pH, DOC, TDN, TDP, chlorophyll-*a*). We analyzed the AA<sub>ESS</sub> carbon isotope values of three primary producer groups and zooplankton consumers as both a community size-fraction (>350um) and individual taxa.  

# File Directory
The file directory contains folders and scripts (Rmd) to be run in RStudio. The folders house data, output, and raw + polished figures.  
   - `Sierra-plankton-CSIA.Rproj` = the R project for the analysis, load this first to allow code to run from correct directory in R Studio
   - `Sierra Nevada Lake Zoops.Rmd` = the scripts and annotated code chunks here will walk through analyses and produce outputs
   - `Sierra Nevada Lake Zoops.html` = the knit output of the Rmd file. Open this in any browser.

### Folders
   - `data` = contains raw and processed data files
   - `figures` = contains raw exported figures from code in R
   - `output` = contains output subfolders - including final figures - and image jpegs 

# Metadata  
Important data files to generate maps, figures, and run models can be found in the *data folder*.  

## Data  
 - `SNL.CSIA.fulldata.csv` = this finalized, compiled output that is generated through the R script. It contains all data, including: metadata, environmental data, raw and mean-normalized δ<sup>13</sup>C isotope values, linear discriminant values, and MixSIAR dietary estimates.  
  - `Sample_metadata.csv` = contains lake, sampleID, and sample type information.  
  - `Sierra Nevada Lakes AA d13C Data NMX standard_update.csv` = contains metadata for GC-C runs, including: amino acid δ<sup>13</sup>C values and SD from average of 2-3 injections, along with standard names and run notes.  
  - `SNL.envdata.csv` = contains ecological sampling information on the lakes, as well as nutrient analysis data
  - `Symons_Yos_bulkSIA.csv` = contains bulk isotope (δ<sup>13</sup>C, δ<sup>15</sup>N) data of fish, zooplankton, POM, grasses, and pines from the Sierra Nevada, previously published in Symons et al. (2019).

### Subfolders    
 - `producers lit` contains data from end-member essential amino acid producers pulled from the literature.
   -  `SNL_wBacteriaOnly.csv` = contains data from present study (algae, POM, terrestrial plants, zooplankton) + bacteria from Larsen et al. (2009) Besser et al. (2025).
   -  `Producers_microalgae_lit.csv` = contains data from present study (algae, POM, terrestrial plants, zooplankton) plus cultured microalgae from Arsenault et al. 2022, Thorp & Bowes 2017
 - `Sample Log` contains run log from lab and field analysis
   - `Sierra_Nevada_Lakes_AA_Sample_Log.xls` = contains notes on GC runs, taxonomy, and water sampling from the study

## Figures  
 - contains main and supplemental figures executed from the code without modification. Note Figures S8-S10 are probability distributions exported from MixSIAR and are found in the `output/final figures` folder.

## Output  
  - `Bact.zoop.LDA.df.csv` = output from linear discriminant analysis (LDA) for zooplankton assignment to bacteria, plants, algae, POM.  
  - `EAA.mean.stdev.summary.table.csv` = contains summary mean and SD for essential amino acids in sources and zooplankton
    
### Subfolders  
 - `Final figures/` = contains the final figures, executed in R and cleaned/aligned using Affinity software  
   
 - `LDA_models/` = contains linear discriminant analysis (LDA) products. 
     - output: `LDA_Train_Info.txt`  
     - loadings: `SNL_Producers_LDA_loadings_NMX.csv` 
       
 - `MixSIAR_models/` = MixSIAR model data and output
     - `EAA_SampleID_model/` = subfolder where MixSIAR writes output
         - contains sample-by-sample posterior distributions, exported csv of mixing model estimates, run stats, model operation output  
         - exported and saved .RDS from JAGS model is here for load in without re-running the model (memory intensive)  
           
     - `loaded data/` = data called into MixSIAR  
         - `MixSIAR_EAA_TDFs.csv` = trophic discrimination factors  
         - `MixSIAR.df.EAA.zoop.csv` = zooplankton essential amino acid δ<sup>13</sup>C data input for MixSIAR  
         - `source.agg.df.EAA.csv`= source essential amino acid δ<sup>13</sup>C data for MixSIAR aggregated as mean and SD  
           
 -  `photos/`
     - `Conness_mid_up.jpg` = image up from Conness Lake used in .Rmd 
     - `Fig1.sitemap.jpg` = exported site map. Code won't produce this for users unless they have a user-API with Google   
     - `sunrise1_lake.jpg` = a lake site with more terrestrial input  

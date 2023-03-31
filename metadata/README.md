The directory (Climate_projection_PlosOne2018) contains the annual climate projection output for 18 GCMs, 3 RCPs, and each of those combinations has a separate Atlantic/GoMex and Pacific file. This directory has a README within.
Climate_projections_PlosOne2018/ has the SST and SBT (?) projections for North America, developed from a SODA climatology with CMIP5 climate models by Thomas Froelicher.

The directory (nicheMods_PlosOne2018) has the GAM models, two for each species. A presence-absence and a biomass GAM for each. 

The hauls_catch_Dec2017.RData file has the raw trawl data. I believe it has two dataframes, one with the environmental data for each haul, and one with the catch data. The two are linked by a unique haulID. 

The speciesProjectionList.RData has the region factor (in the regionFreq vector) applied for each species projection (named in the projspp vector). The base region was based on where the most occurrences occurred. In other words, the species projection is in units relevant to the catchability of the region factor that was used (this corresponds to a survey, like NEUS Spring).

CEmodels_proj_PresAbs_May2018 has all the projection output for the presence-absence GAMs alone (no biomass model). I have three RCPs (26,45,85) and 18 GCMs for each species. This is what we're using for the economic analysis, and is one of the niche model types I'm using for my current analysis on uncertainty. In order, the climate models (in columns) are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, HadGEM2-ES, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MIROC5, MPI-ESM-LR, NorESM1-ME

CEmodels_proj_Biomass_BCODMO has all the projection output for the biomass GAMs (pres-abs and biomass together). This is from the files on BCO-DMO, since the files on Data Dryad seem corrupt. In order, the climate models are bcc-csm1-1-m, bcc-csm1-1, CanESM2, CCSM4, CESM1-CAM5, CNRM-CM5, GFDL-CM3, GFDL-ESM2M, GFDL-ESM2G, GISS-E2-R, GISS-E2-H, IPSL-CM5A-LR, IPSL-CM5A-MR, MIROC-ESM, MPI-ESM-LR, NorESM1-ME

"Temp proj files" also has climate projection data. Not sure how this relates to Climate_projections_PlosOne2018/

Code is at https://github.com/pinskylab/project_velocity
Code as used for Morley et al. 2018 is https://github.com/pinskylab/project_velocity/releases/tag/v1.1



README.md written by Jim Morley

Edited by Malin Pinsky

Uploaded by Ryan Snow

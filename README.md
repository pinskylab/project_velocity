README
================

project\_velocity Projecting shifts in habitat for 686 North American
marine species under RCP 2.6 and 8.5.

Publication: Morley JW, Selden RL, Latour RJ, Frolicher TL, Seagraves
RJ, Pinsky ML (2018) Projecting shifts in thermal habitat for 686
species on the North American continental shelf. PLoS ONE 13(5):
e0196127. <https://doi.org/10.1371/journal>. pone.0196127

Contact: James Morley University of North Carolina, Chapel Hill
<jwmorley@email.unc.edu>

Malin Pinsky Rutgers University <malin.pinsky@rutgers.edu>

# Directories

## data:

Contains Rdata files used in the R-scripts. These include the raw
temperature projection data in the sub-directory ‘Temp proj files’ and
the formatted temperature projection files in the subdirectory
‘prediction\_files\_Nov2017’.

## figures:

Contains the figures from the projections project.

# Analysis workflow

Note: There are some older scripts here from initial work on this
manuscript, which were not ultimately used. This is noted where
appropriate.

## 1\_data\_allcombine.r

Combine the trawl survey data together and reformat in preparation for
creating climatologies and climate envelope models

## 1\_data\_trawlData.R

An alternative method to combine all the trawl survey data with the
‘trawlData’ package. This is the script that was ultimately used for
the manuscript.

## 1.5\_merge\_hadl\_noaa\_sst.r

Fill in missing SST values for a portion of the trawl hauls, using
satellite temperature observations from the Hadley NOAA/OI

## 2\_create\_climatologies.r

Only used for an earlier version of the manuscript and not for published
version. Because T. Froelicher completed this independently. This
created the historic climatology within each survey region based on in
situ temperature measurements.

## 3\_calculate\_delta.r

Only used for an earlier version of the manuscript and not for published
version. Because T. Froelicher completed this independently. This
calculate the delta (i.e., change in temperature) values within each
grid cell for each projection year, 2007-2100.

## 4\_climateprojection.r

Only used for an earlier version of the manuscript and not for published
version. Because T. Froelicher completed this independently. This
integrates the historic climatology from script 2 with the delta values
from script 3 to compute the temperature projections.

## 4.5\_summarize\_climprojections.r

Only used for an earlier version of the manuscript and not for published
version. Because T. Froelicher completed this independently. Basic
summary statistics and plots for the temperature projections.

## 5\_select\_species.r

Standardizes species names across surveys based on some files that have
species name corrections. Also corrects some species names. Drops
species/life stages that are not of interest (e.g., larval fish). Picks
the species that are ultimately included in the analysis based on
frequency of occurence.

## Projection\_grid.R

Develops the projection grid based on the resolution of the SODA climate
data. This projection grid was used to conduct the future climate
projection and to integrate with the bathymetry environmental data in
‘Bathymetry.R’.

## Bathymetry.R

Depth data is brought in from GEBCO database and a depth data grid is
developed and rugosity is calculated based on depth. Then sediment data
is brought in from a data set that was put together elsewhere (contact
for that is Becca Selden). Finally, the projection grid was developed
based on the SODA climate resolution, in other words, we determined
where these three data sources all overlap and are shallower than 400 m.
Overall, this script is not cleaned up well.

## 5.4\_master\_haul\_file.R

Merge together a number of files containing the environmental data for
hauls, which creates a master file for the trawl-hauls.

## 6\_model\_fitting\_loop.R

Fits the niche model to each species and calculates some sdm model
evaluation statistics.

## 7\_Prep\_Projections.R

Drop species with poorly fit habitat models and then reformat the region
value for conducting projections. Then reads in and formats all the
climate projection data into a data structure that can be used to
predict with the habitat models. Also produces some figures for the
climate projection and some code to qaqc the climate data.

## 7b\_projection\_clusterB\_PA.R

Conducts projections based on the probability of occurrence GAM model
only (i.e., not what was used in the Plos One paper). The script
integrates the habitat models with the climate projection models and
then saves the aggregated (20 year time bins) output. This script is
activated by a shell script, which was used on a CPU cluster.

## 7b\_projection\_clusterB.R

Conducts projections based on the delta biomass GAM model. The script
integrates the habitat models with the climate projection models and
then saves the aggregated (20 year time bins) output. This script is
activated by a shell script, which was used on a CPU cluster.

## 8\_summarize\_projections.R

Some of Malin’s old code from the first round of projections, which were
done on a regional basis. This script was not ultimately used on the
Plos One (2018) paper.

## 8.1B\_plot\_projection\_summaries.R

Some of Malin’s old code from the first round of projections, which were
done on a regional basis. A bunch of figures are made here with the
projections output. This script was not ultimately used on the Plos One
(2018) paper.

## 8.5\_summarize\_projections\_separate\_models.r

Some of Malin’s old code from the first round of projections, which were
done on a regional basis. This script was not ultimately used on the
Plos One (2018) paper.

## abundance\_model\_testing.r

Tests how much variation/deviance is explained by temperature in
abundance GAMs. We did not report these output in the Plos One paper.

## MS\_figures.R

Summarizes the projection output and makes the figures used in the Plos
One (2018) manuscript.

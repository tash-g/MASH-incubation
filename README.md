# Resource allocation underlies parental decision-making during incubation in the Manx shearwater
Natasha Gillies, Oliver Padget, Martyna Syposz, Sarah Bond, Tim Guilford

## Overview
This folder contains codes and data to reproduce both the main and supplementary results of the paper, published in Ornithology (2022); https://doi.org/10.1093/ornithology/ukac006. 

There are 2 folders: "Data_inputs" and "Scripts". As this is a large and multi-year dataset that may be currently used by other members of the project or collaborators, we would appreciate if you could contact Natasha Gillies (gilliesne@gmail.com) if you would like to make use of the dataset. This will reduce duplication of effort and ensure that we can give you any additional information that may be useful for analysis or interpretation.

## Scripts
- **1_MASH-incubation_main-analysis.R** This script contains all the code necessary to reproduce the main results and figures presented in the manuscript.
- **2_MASH-incubation_functions.R** This script contains additional functions to run the analyses in script 1. 

## Data inputs
There are 6 datasets used in the analysis.

- **MASH2015-2019_daily-active-nests.csv** Contains the number of active nests per day and year, and is used to define the simulated changeovers for the randomisation analysis. The columns are:
  -	_date_ Date where active nests recorded
  -	_nburrows_ Number of active nests
  -	_year_ Year of recording

- **MASH2015-2019_daily-foraging-data.csv** Summarised activity for birds across individual foraging trips. The columns are:
  -	_ring_ Individual ring identity of bird
  -	_stintID_ Unique identity of individual stint
  -	_sex_ Female (F) or male (M)
  -	_burrow_ Identity of nest
  -	_start_mass_ Mass (g) at beginning of foraging trip
  -	_year_ Year of recording
  -	_end_mass_ Mass (g) at end of foraging trip
  -	_foragetrip_ Duration (days) of foraging trip
  -	_forageGain_ Mass (g) gained across duration of foraging trip
  -	_neglect_ Whether the trip followed neglect (YES) or not (NO)
  -	_prop_rest_ Proportion of foraging trip spent resting
  -	_prop_flight_ Proportion of foraging trip spent in flight
  -	_prop_forage_ Proportion of foraging trip spent foraging
  -	_phase_ What phase of the trip the bird was in, commuting (commute) or foraging (forage)
	
- **MASH2015-2019_daily-mass-change.csv** Contains daily mass changes for all birds for which mass data were collected. The columns are:
	-	_year_ Year data was collected
	-	_ring_ Individual identity of bird
	-	_sex_ Female (F) or male (M)
	-	_stintID_ Denotes unique identity of stint (within bird)
	-	_DSL_ ‘Days Since Laying’, indicates time since egg was laid
	-	_stintday_ Day of stint
	-	_mass_chg_ Difference in mass (g) from previous day
	-	_mass_chgPer_ Difference in mass (g) from previous day as a percentage of body mass
	-	_exp_ Pair experience, unknown (UNKNOWN), newly formed (NEW) or experienced (EXP)

- **MASH2015-2019_forage-data-by-trip.csv** Contains summarised data from individual foraging trips. The columns are:
	-	_burrow_ Individual identity of nest
	-	_year_ Year data were collected
	-	_stintID_ Denotes unique identity of stint (within bird)
	-	_ring_ Individual identity of bird
	-	_sex_ Female (F) or male (M)
	-	_start_mass_ Mass of bird (g) at beginning of foraging trip
  -	_DSL_ ‘Days Since Laying’, indicates time since egg was laid
  -	_foragetrip_ Duration of foraging trip (days)
  -	_partnerincoming_ Mass of returning partner (g) upon departure on foraging trip
  -	_exp_ Pair experience, unknown (UNKNOWN), newly formed (NEW) or experienced (EXP)

- **MASH2015-2019_incub-data-by-shift.csv** Contains summarised data from individual incubation shifts. The columns are:
  -	_burrow_ Individual identity of nest
  -	_year_ Year data were collected
  -	_stintID_ Denotes unique identity of stint (within bird)
  -	_start_date_ Date that the incubation shift began
  -	_ring_ Individual identity of bird
  -	_sex_ Female (F) or male (M)
  -	_start_mass_ Mass of bird (g) at beginning of foraging trip
  -	_burrowID_ Denotes identity of nest within a given year
  -	_stintlength_ Duration of incubation shift (days)
  -	_DSL_ ‘Days Since Laying’, indicates time since egg was laid
  - _end_mass_ Mass of bird (g) at end of incubation shift
  -	_endsInNeglect_ Whether bird neglected egg at end of stint (YES) or not (NO)
  -	_followsNeglect_ Whether incubation shift commenced following neglect (YES) or not (NO)
  -	_exp_ Pair experience, unknown (UNKNOWN), newly formed (NEW) or experienced (EXP)
  -	_whoFirst_ Factor denoting whether the female (F) or male (M) took the first incubation shift. If the first incubation shift was not recorded for that nest, recorded as ‘nostart’

- **MASH2015-2019_nest-data.csv** Contains summarised data for entire incubation attempts for each nest used in the study. The columns are:
	-	_year_ Year the data were collected
	-	_burrow_ Denotes individual identity of nest
	-	_outcome_ Whether the incubation attempt was successful (HATCHED) or not (FAILED)
	-	_incDur_ Total duration (days) from egg laying to egg hatching (or failure)
	-	_count.F_ Number of incubation shifts taken by the female
	-	_count.M_ Number of incubation shifts taken by the male
	-	_NO.BIRD_ Number of days neither parent was recorded at the nest
	-	_femper_ Proportion of incubation days taken by female
	-	_malper_ Proportion of incubation days taken by male
  -	_exp_ Pair experience, unknown (UNKNOWN), newly formed (NEW) or experienced (EXP)
  -	_julian_ Lay date recorded as Julian date
  -	_whoFirst_ Factor denoting whether the female (F) or male (M) took the first incubation shift. If the first incubation shift was not recorded for that nest, recorded as ‘nostart’



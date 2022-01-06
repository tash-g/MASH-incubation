
### This script contains all the necessary code to run the main analyses and produce
### the figures presented in the main text. 

# PREAMBLE ----------------------------------------------------------------

### Create outputs folder if one doesn't exist already
out.path <- "./Data_outputs/"
if(dir.exists(out.path) == FALSE){
  dir.create(out.path)
}

### Create figures folder if one doesn't exist already
out.path <- "./Figures/"
if(dir.exists(out.path) == FALSE){
  dir.create(out.path)
}

### Set figure colours
female_col <- "darkorange3" # Blood orange
male_col <- "#23a7bc"   # Turquoise


# ~ LOAD DATASETS ~ -------------------------------------------------------

## Daily foraging activity dataset 
daily <- read.csv("Data_inputs/MASH2015-2019_daily-foraging-data.csv")
daily$year <- as.factor(daily$year)

## Dataset of each incubation stint 
incubdf <- read.csv("Data_inputs/MASH2015-2019_incub-data-by-shift.csv")
incubdf$start_date <- as.Date(incubdf$start_date)
incubdf$sex <- as.factor(incubdf$sex)
incubdf$endsInNeglect <- as.factor(incubdf$endsInNeglect)

## Dataset of daily mass changes during incubation
masschanges <- read.csv("Data_inputs/MASH2015-2019_daily-mass-change.csv")

## Dataset of nest data
nestdf <- read.csv("Data_inputs/MASH2015-2019_nest-data.csv")
nestdf$year <- as.factor(nestdf$year)

## Dataset of daily active nests
activeNests <- read.csv("Data_inputs/MASH2015-2019_daily-active-nests.csv")

## Dataset of each foraging trip
foragedf <- read.csv("Data_inputs/MASH2015-2019_forage-data-by-trip.csv")
foragedf$year <- as.factor(foragedf$year)
foragedf$start_mass <- as.numeric(foragedf$start_mass)
foragedf$partnerincoming <- as.numeric(foragedf$partnerincoming)
foragedf$sex <- as.factor(foragedf$sex)
foragedf$ring <- as.factor(foragedf$ring)
foragedf$burrow <- as.factor(foragedf$burrow)


# ~ LOAD PACKAGES ~ -----------------------------------------------------------

# Packages
packages <- c("ggplot2", "lme4", "lmerTest", "emmeans", "dplyr",  "data.table", "reshape2", 
              "glmmTMB", "sjPlot", "effects", "rmcorr", "mgcv")

# Install packages not yet installed
installed_packages <- packages %in% rownames(installed.packages())

if (any(installed_packages == FALSE)) {
  install.packages(packages[!installed_packages], lib = "C:/Users/ngillies/OneDrive - The University of Liverpool/Documents/R/win-library")
}

# Load packages
invisible(lapply(packages, library, character.only = TRUE))

# Load functions
source("Scripts/2_MASH-incubation_functions.R")



# PATTERNS IN INCUBATION ------------------------

## Get sample sizes
length(unique(incubdf$burrow)) # 74 nests
incubdf %>% dplyr::group_by(year) %>% dplyr::summarise(unique_burrows = n_distinct(burrow))

x <- unique(subset(incubdf, year == "2015" | year == "2016" | year == "2017" | year == "2018")$burrow)
y <- unique(subset(incubdf, year == "2019")$burrow)

# ~ LMM 1 # Predictors of lay date ~ ---------------------------------------

# Remove nests where experience unknown
nest.df_Exp <- subset(nestdf, !is.na(exp) & exp != "UNKNOWN")

# Get sample size
nrow(nest.df_Exp); length(unique(nest.df_Exp$burrow))

# Build linear model to predict lay date
lmmLAY <- lmer(julian ~ exp + year + (1|burrow), data = nest.df_Exp)

drop1(lmmLAY, ddf = "lme4", test = "Chisq")
summary(lmmLAY)
emmeans(lmmLAY, "exp", type = "response")
emmeans(lmmLAY, "year", type = "response")

# Get mean lay date
as.Date("2021-01-01") + (as.numeric(round(summary(lmmLAY)$coefficients[1]))) # 2021-05-15


# ~~ TABLE S1 ~~ ----------------------------------------------------------------

## Get lay dates for each year
as.Date("2021-01-01")+130; as.Date("2021-01-01")+127; as.Date("2021-01-01")+133


# ~~ FIGURE S1 ~~ ---------------------------------------------------------------

## Create boxplot of lay date by year
datebreaks = c(115, 120,125,130,135,140,145,150,155,160)
datelabels = format(as.Date("2021-01-01")+datebreaks, "%d-%m")

png("Figures/FIGURES1_layDate-year.png", height = 6, width = 9, units = 'in',res = 300)
ggplot(nest.df_Exp, aes(x = year, y = julian)) +
  geom_boxplot2() +
  theme_custom() +
  labs(x = "Year", y = "Lay date") +
  scale_y_continuous(breaks = datebreaks, labels = datelabels)
dev.off()
  


# ~ GLMM 2 (poisson) # Incubation duration ~ ----------------------------------------------

## Remove outlier identified as skewing results
nest.df_noOutlier <- subset(nest.df_Exp, incDur < 65) 

## Get sample size
nrow(nest.df_noOutlier); length(unique(nest.df_noOutlier$burrow)) # 73 nests, 58 unique

## Get means, sd, min, and max incubation duration for successful and failed nests
aggregate(nest.df_noOutlier$incDur, list(nest.df_noOutlier$outcome), 
          FUN=function(x) c(mn = mean(x), std = sd(x), mini = min(x), maxi = max(x))) 

#Group.1      x.mn     x.std    x.mini    x.maxi
#1  FAILED 36.214286 11.490277 16.000000 67.000000
#2 HATCHED 50.355556  1.653585 47.000000 54.000000

## Truncated poisson to model incubation duration as a function of outcome, experience and year
glmeINCDUR <- glmer(incDur ~ outcome + exp + year + (1|burrow), data = nest.df_noOutlier,
                      family = "poisson")

drop1(glmeINCDUR, ddf = "lme4", test = "Chisq")
summary(glmeINCDUR)
emmeans(glmeINCDUR, "exp", type = "response")
emmeans(glmeINCDUR, "outcome", type = "response")


# ~ GLMM 3 (beta) # Proportional shares of incubation ~ --------------------------

## Melt dataset to extract proportional shares for males and females in each nest
propVars <- c("femper", "malper", "burrow", "year", "outcome", "exp", "incDur", "whoFirst")
nestdf_Prop <- nest.df_Exp[propVars]
nestdf_Prop <- melt(nestdf_Prop, id = c("burrow", "year", "outcome", "exp", "incDur", "whoFirst"))
nestdf_Prop <- subset(nestdf_Prop, whoFirst != "nostart") # Remove nests where incubation data incomplete
nestdf_Prop$value <- as.numeric(nestdf_Prop$value)
colnames(nestdf_Prop)[c(7,8)] <- c("sex", "propShare")

# Get sample sizes
nrow(nestdf_Prop); length(unique(nestdf_Prop$burrow)) # 124 nests, 52 unique

## Construct a beta mixed model (as using proportional data which cannot be 0 or 1) to analyse whether
## incubation shares vary with sex, experience, or year
glmmSHARE <- glmmTMB(propShare ~ sex + exp + year + (1|burrow),
                     data = nestdf_Prop,
                     beta_family())

drop1(glmmSHARE, test = "Chisq")
summary(glmmSHARE)
emmeans(glmmSHARE, "sex", type = "response")


# ~ GLMM 4 (poisson) # Number of incubation shifts ~ ----------------------------------

## Melt dataset to extract number of incubation shifts for males and females in each nest
countVars <- c("count.F", "count.M", "burrow", "year", "outcome", "exp", "incDur", "whoFirst")
nestdf_Count <- nest.df_Exp[countVars]
nestdf_Count <- melt(nestdf_Count, id = c("burrow", "year", "outcome", "exp", "incDur", "whoFirst"))
nestdf_Count <- subset(nestdf_Count, whoFirst != "nostart") # Remove nests where incubation data incomplete
colnames(nestdf_Count)[c(7,8)] <- c("sex", "stintCount")

## Get sample size
nrow(nestdf_Count); length(unique(nestdf_Count$burrow)) # 124 nests, 52 unique

## Construct a Poisson GLMM to analyse whether number of incubation shifts varies between sexes, with
## pair experience or year. Include incubation duration to control for differences in time spent
## incubating
glmmSTINTNO <- glmer(stintCount ~ sex + exp + incDur + year + (1|burrow), data = nestdf_Count,
                      family = "poisson")

drop1(glmmSTINTNO, ddf = "lme4", test = "Chisq")
summary(glmmSTINTNO)
emmeans(glmmSTINTNO, "sex", transform = "response")
emmeans(glmmSTINTNO, "exp", transform = "response")


# ~ GLMM 5 (poisson) # Duration of first incubation shift ~ --------------------------

## Isolate first incubation shift for each nest
incubdf_First <- subset(incubdf, DSL == 1)
incubdf_First <- subset(incubdf_First, !is.na(exp) & exp != "UNKNOWN")

## Get sample sizes
nrow(incubdf_First); length(unique(incubdf_First$burrow)) # 68 nests, 52 unique

## How many first stints taken by males versus females
length(unique(incubdf$burrowID[incubdf$whoFirst == "M"])) ## 44
length(unique(incubdf$burrowID[incubdf$whoFirst == "F"])) ## 34

# Compare female number to null expectation of 0.5
binom.test(34, 78, p = 0.5, alternative = "two.sided")


## Construct poisson GLMM to analyse whether first stint duration differs between
## sexes, with experience, and with year
glmmSEXDUR_first <- glmer(stintlength ~ sex + exp + year + (1|burrow/ring), data = incubdf_First,
                            family = "poisson")

summary(glmmSEXDUR_first)
drop1(glmmSEXDUR_first, ddf = "lme4", test = "Chisq")
summary(glmmSEXDUR_first)
emmeans(glmmSEXDUR_first, "sex", type = "response")
emmeans(glmmSEXDUR_first, "exp", type = "response")
emmeans(glmmSEXDUR_first, "year", type = "response")


# GLMM 6 (poisson) # Duration of all incubation shifts ---------------------------

durVars <- c("stintlength", "sex", "DSL", "exp", "year", "burrow", "ring")
incubdf_Dur <- incubdf[durVars]
incubdf_Dur <- na.omit(incubdf_Dur)
incubdf_Dur <- subset(incubdf_Dur, !is.na(exp) & exp != "UNKNOWN") # Remove nests with unknown experience
incubdf_Dur <- subset(incubdf_Dur, DSL > 1) # Isolate complete incubation periods

# Get samle size
nrow(incubdf_Dur); length(unique(incubdf_Dur$burrow)) # 459 shifts, 58 unique nests

## Construct poisson GLMM to examine whether stint duration varies between the sexes, 
## with egg age (DSL), pair experience, or year
glmmINCDUR <- glmer(stintlength ~ sex + DSL + exp + year + (1|burrow/ring), 
                      data = incubdf_Dur, family = "poisson")

drop1(glmmINCDUR, ddf = "lme4", test = "Chisq")
summary(glmmINCDUR)
emmeans(glmmINCDUR, "sex", type = "response")
emmeans(glmmINCDUR, "exp", type = "response")


# ~ Q # How correlated are trips? ~ -----------------------------------------

corVars <- c("burrow", "year", "stintID", "stintlength")
incubdf_Cor <- incubdf[corVars]

# Construct forward correlation to determine how correlated is stint with the next stint 
incubdf_Cor <- data.frame(   incubdf_Cor %>% 
                         dplyr::group_by(burrow, year, stintID) %>%
                         dplyr::mutate(stintA = stintlength[1],
                                       stintB = stintlength[2]) %>%
                         filter(!duplicated(stintID)))

incubdf_Cor <- subset(incubdf_Cor, !is.na(stintB))
incubdf_Cor$burrowID <- as.factor(paste(incubdf_Cor$burrow, incubdf_Cor$year, sep = "_"))

## Obtain Pearson's correlation coefficient
rmcSHIFT <- rmcorr(burrowID, stintA, stintB, incubdf_Cor)
print(rmcSHIFT)

# ~ Q # Is colony nest attendance random? ~ -----------------------

## Construct dataset where each row indicates a changeover, using incubdf dataset

# Remove first observation - not changeover
incubdf.Chgeovr <- data.frame ( incubdf %>% 
                                  dplyr::group_by(burrowID) %>% 
                                  slice(2:n()) )   

# Remove any stints ending in or following neglect (not changeover)
incubdf.Chgeovr <- subset(incubdf.Chgeovr, endsInNeglect == "NO" & followsNeglect == "NO")

## Get distribution of changeover frequencies from on active nests dataset
activeNests_Dist <- data.frame (  activeNests %>%
                          dplyr::group_by(nburrows, year) %>%
                          dplyr::summarise(Freq = n())  )

## Extract frequencies 
freq.true <- rep(activeNests_Dist$nburrows, activeNests_Dist$Freq)
years.true <- rep(activeNests_Dist$year, activeNests_Dist$Freq)

# Distribution of changeover frequencies across years
freq.years <- data.frame(cbind(freq.true, years.true))


### Randomisation to compare colony changeovers ###

stintLengths <- incubdf.Chgeovr$stintlength
endPoints <- nestdf$incDur
burrowSet <- unique(incubdf$burrowID)
pList <- list()
histList <- list()

for(y in 1:length(unique(incubdf$year))) {
  
  # Get active burrows for target year
  burrows <- burrowSet[grepl(unique(incubdf$year)[y], burrowSet) == TRUE]
  
  for(k in 1:10000) {
    
    ### For each active burrow, simulate a set of incubation stints
    simStints <- list()
    
    for(i in 1:length(burrows)) {
      
      # Randomly select start point and end point
      startpoint <- subset(incubdf, burrowID == burrows[i])$start_date[1]
      endpoint <- sample(endPoints, 1)
      stints <- list()
      
      cumval <- 0
      j <- 1
      
      while (cumval < endpoint) {
        
        ## Sample stints to fill randomised incubation period
        stints[j] <- sample(stintLengths, 1)
        cumval <- sum(unlist(stints))
        j <- j + 1
        
      }
      
      startpoint <- rep(startpoint, length(stints))
      simDF <- as.data.frame(cbind(as.character(unlist(startpoint)), unlist(stints)))
      colnames(simDF) <- c("start_date", "simStints")
      simDF$start_date <- as.Date(as.character(simDF$start_date))
      simDF$simStints <- as.numeric(simDF$simStints)
      
      simDF$cumstints <- cumsum(simDF$simStints)
      simDF$start_date <- ifelse(!is.na(lag(simDF$cumstints)), 
                                 as.character(simDF$start_date+lag(simDF$cumstints)),
                                 as.character(simDF$start_date))
      
      simDF$cumstints <- NULL
      
      simStints[[i]] <- simDF
      
    }
    
    simStints <- data.table::rbindlist(simStints)
    
    ### Calculate simulated distribution
    simChanges <- data.frame(table(simStints$start_date))
    simChanges$Var1 <- as.Date(simChanges$Var1)
    
    # Add in dates with no changeovers
    fullDates.sim <- data.frame(simChanges %>% 
                                  dplyr::mutate(d = list(seq(as.Date(min(Var1)), as.Date(max(Var1)), by = "1 day"))) %>%
                                  tidyr::unnest() ) [,3]
    
    fullDates.sim <- unique(fullDates.sim)
    fullDates.sim <- data.frame(fullDates.sim)
    colnames(fullDates.sim) <- "Var1"
    fullDates.sim$Var1 <- as.Date(fullDates.sim$Var1)
    
    # Merge to changes dataframe
    simChanges <- merge(fullDates.sim, simChanges, by = "Var1", all.x = T)
    simChanges$Freq <- ifelse(is.na(simChanges$Freq), 0, as.numeric(simChanges$Freq))
    colnames(simChanges)[2] <- "nburrows"
    
    # Get distribution
    distDF.sim <- data.frame (  simChanges %>%
                                  dplyr::group_by(nburrows) %>%
                                  dplyr::summarise(Freq = n())  )
    
    
    ### Get mean of distribution 
    
    ## Get frequencies
    freq.sim <- rep(distDF.sim$nburrows, distDF.sim$Freq) 
    
    ## Mean
    #pList[[k]] <- mean(freq.sim)
    pList[[k]] <- length(freq.sim[freq.sim==0])/length(freq.sim)
    print(k)
    
  }
  
  trueVal <- nrow(freq.years[freq.true == 0 & freq.years$years.true == unique(incubdf$year)[y],])/nrow(freq.years[freq.years$years.true == unique(incubdf$year[y]),])
  ranVal <- mean(unlist(pList))
  diffr <- abs(trueVal - ranVal)
  
  pright <- length(unlist(pList[pList < ranVal-diffr]))
  pleft  <- length(unlist(pList[pList > ranVal+diffr])) 
  
  hist(unlist(pList), xlim = c(0,5))
  abline(v = trueVal, col = "red")  
  
  histList[[y]] <- c(pleft, pright, unique(incubdf$year)[y])
  
}

## Extract output of simulated changeover frequencies
sim.output <- data.frame(do.call("rbind",histList))
colnames(sim.output) <- c("pleft", "pright", "year")
sim.output[1:2] <- lapply(sim.output[1:2], as.character)
sim.output[1:2] <- lapply(sim.output[1:2], as.numeric)

## Get p value
1-sum(sim.output$pleft, sim.output$pright)/50000



# MASS CHANGES DURING INCUBATION -------------------------------------------

# ~ LMM 7 # Starting mass ~ -----------------------------------------------

massVars <- c("start_mass", "DSL", "sex", "exp", "year", "ring", "burrow")
incubdf_Mass <- incubdf[massVars]
incubdf_Mass <- na.omit(incubdf_Mass)
incubdf_Mass <- subset(incubdf_Mass, !is.na(exp) & exp != "UNKNOWN")
incubdf_Mass <- incubdf_Mass[order(incubdf_Mass$ring, incubdf_Mass$year, incubdf_Mass$DSL),]

## Get sample sizes
nrow(incubdf_Mass); length(unique(incubdf_Mass$ring)); length(unique(incubdf_Mass$burrow)) # 471 stints, 103 individuals, 52 nests

## Construct LMM to examine effect of egg age (DSL), sex, experience, and year on mass at the 
## start of the incubation stint
lmmSTINTM <- lmer(start_mass ~ DSL + sex + exp + year + (1|ring), data = incubdf_Mass)

drop1(lmmSTINTM, ddf = "lme4", test = "Chisq")
summary(lmmSTINTM)
confint(lmmSTINTM)
emmeans(lmmSTINTM, "sex")
emmeans(lmmSTINTM, "exp")
emmeans(lmmSTINTM, "year")

# Get difference in starting mass between last and first incubation stint
shiftSum <- incubdf_Mass %>% group_by(ring, year) %>%
  mutate(first = first(start_mass), last = last(start_mass), change = last - first,
         ringYear = paste(ring, year)) %>%
  distinct(ringYear, .keep_all = T) %>%
  data.frame() 

mean(shiftSum$change)


# ~~ FIGURE S2 ~~ ---------------------------------------------------------

lmmSTINTM_pdat <- plot_model(lmmSTINTM, type = "pred", terms = c("DSL", "sex"))
lmmSTINTM_pdat <- data.frame(lmmSTINTM_pdat$data)

png("Figures/FIGURES2_startMass-DSL.png", height = 6, width = 9, units = 'in',res = 300)
ggplot() + 
  geom_point(aes(x = DSL, y = start_mass, col = sex), data = incubdf_Mass) + 
  geom_ribbon(data = lmmSTINTM_pdat, aes(x = x, ymin = conf.low, ymax = conf.high, group = group_col),
              alpha = 0.5, fill = "grey") +
  geom_line(data = lmmSTINTM_pdat, aes(x = x, y = predicted, col = group_col), size = 1) +
  theme_custom() + theme(legend.position = c(0.90, 0.85)) + labs(col = "Sex") +
  labs(x = "Days since laying", y = "Mass at start of incubation shift (g)") +
  scale_color_manual(values = c(female_col, male_col), labels = c("Female", "Male")) +
  guides(color=guide_legend(override.aes=list(fill=NA)))
dev.off()


# ~ LMM 8 # Daily mass decline ~ ------------------------------------------

## Get mean mass changes in grams and as a percentage of body mass
mean(masschanges$mass_chg, na.rm = T); sd(masschanges$mass_chg, na.rm = T) # -10.02 +/- 7.42g
mean(masschanges$mass_chgPer, na.rm = T); sd(masschanges$mass_chgPer, na.rm = T) # - 0.02 +/- 0.02%

## Extract relevant variables for model
chgVars <- c("mass_chgPer", "stintday", "sex", "exp", "DSL", "stintID", "year", "ring")
maschanges_Chg <- masschanges[chgVars]
maschanges_Chg <- na.omit(maschanges_Chg)

# Make mass change positive for transformation and remove positives (not real)
maschanges_Chg <- subset(maschanges_Chg, mass_chgPer <= 0)
maschanges_Chg$mass_chgPer <- maschanges_Chg$mass_chgPer*-1
maschanges_Chg <- subset(maschanges_Chg, exp != "UNKNOWN")

## Get sample sizes
nrow(maschanges_Chg); length(unique(maschanges_Chg$ring)) # 2186 days, 102 individuals

## Construct beta mixed GLMM to analyse effect of experience, sex, year, egg age, and day of
## stint on daily mass declines
glmmCHANGES <- glmmTMB(mass_chgPer ~  exp + sex + year + DSL + stintday + 
                     (1|ring), data = maschanges_Chg, ziformula = ~1, family = beta_family())

drop1(glmmCHANGES, test = "Chisq")
summary(glmmCHANGES)
emmeans(glmmCHANGES, "sex", type = "response")
confint(glmmCHANGES, method = "Wald")


# ~ GLMM 9 # Probability shift ends in neglect ~ ------------------------------

# Extract relevant data for model
negVars <- c("start_mass", "stintlength", "DSL", "year", "sex", "exp", "endsInNeglect", "ring")
incubdf_Neg <- incubdf[negVars]
incubdf_Neg$neglect <- ifelse(incubdf_Neg$endsInNeglect == "NO", 0, 1)
incubdf_Neg$neglect <- as.factor(incubdf_Neg$neglect)

incubdf_Neg <- na.omit(incubdf_Neg)
incubdf_Neg <- subset(incubdf_Neg, !is.na(exp) & exp != "UNKNOWN")

## Get sample size
nrow(incubdf_Neg); length(unique(incubdf_Neg$ring)) # 471 stints for 103 individuals

## Construct binomial GLMM to examine effect of start mass, stint length, egg age (DSL), sex
## experience, and year on the probability that a stint ends in neglect
glmeNEG <- glmer(neglect ~ start_mass + stintlength + DSL + sex + exp +  year + (1|ring), 
                 family = binomial("logit"), glmerControl(optimizer = "bobyqa", 
                                                          optCtrl = list(maxfun = 100000)), 
                 data = incubdf_Neg)

drop1(glmeNEG, test = "Chisq", ddf = "lme4")
summary(glmeNEG)
emmeans(glmeNEG, "sex", type = "response")
emmeans(glmeNEG, "exp", type = "response")
glmeNEG_confint <- confint(glmeNEG, method = "Wald")

exp(summary(glmeNEG)$coefficients)
exp(glmeNEG_confint)


# ~~ FIGURE 1 ~~ ----------------------------------------------------------

## Plot probability of neglect as a function of starting mass

glmeNEG_pdat <- plot_model(glmeNEG, type = "pred", vars = c("start_mass"), 
              show.ci = FALSE, show.scatter=FALSE)
glmeNEG_pdat <- glmeNEG_pdat$start_mass$data # gives predicted valus and CIs for each value of mass

incubdf_Neg$neglect <- as.numeric(as.character(incubdf_Neg$neglect))

png("Figures/FIGURE1_neglect_mass.png", height = 6, width = 9, units = 'in',res = 300)
ggplot(aes(x = x, y = predicted), data = glmeNEG_pdat) +
  geom_line(size = 1, col = "blue") +
  geom_count(aes(x = start_mass, y = neglect), alpha = 0.4, data = incubdf_Neg) +
  scale_size(name = expression(paste(italic("n"), " nests")), breaks = c(1,5,10,20)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), data= glmeNEG_pdat, alpha = 0.2) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
        axis.title.y = element_text(colour = "black", size = 17, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
        axis.title.x = element_text(colour = "black", size = 17,margin = margin(t = 12, r = 0, b =0, l = 0)),
        legend.position = c(0.93,0.85)) +
  labs(y = "Probability of neglect", x = "Starting mass (g)")
dev.off()


# ~ GLMM 9B # Time to neglect ~ ---------------------------------------------------------

## Extract stints that ended in neglect
incubdf_Neg.cut <- subset(incubdf_Neg, neglect == 1)

## Get sample size
nrow(incubdf_Neg.cut); length(unique(incubdf_Neg.cut$ring)) # 43 stints for 21 individuals

## Construct Poisson GLMM to analyse effect of start mass, egg age (DSL) and sex on stint
## duration for those ending in neglect
glmeShiftNeg <- glmer(stintlength ~ start_mass + DSL + sex + (1|ring), 
                    data = incubdf_Neg.cut, family = "poisson")

summary(glmeShiftNeg)
drop1(glmeShiftNeg, test = "Chisq")
glmeShiftNeg_confint <- confint(glmeShiftNeg, method = "Wald")
exp(glmeShiftNeg_confint)


# ~~ FIGURE 2 ~~ ----------------------------------------------------------

glmeShiftNeg_pDat <- plot_model(glmeShiftNeg, type = "pred", terms = c("start_mass", "sex"))
glmeShiftNeg_pDat <- data.frame(glmeShiftNeg_pDat$data)

# Plot the relationship 
png("Figures/FIGURE2_neglectShift_mass.png", height = 6, width = 9, units = 'in',res = 300)
ggplot(aes(x = x, y = predicted), data = glmeShiftNeg_pDat) +
  geom_line(size = 1, aes(col = group)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high, group = group), data = glmeShiftNeg_pDat, alpha = 0.1) + 
  geom_point(aes(x = start_mass, y = stintlength, colour = sex), data = incubdf_Neg.cut) +
  theme_bw(base_size =15) +
  theme_custom() + theme(legend.position = c(0.1,0.9)) +
  scale_colour_manual(values = c(female_col, male_col), labels = c("Female", "Male"), name = "Sex") +
  labs(x = "Mass at start of incubation shift (g)", y = "Time to neglect (days)") +
  scale_y_continuous(breaks = c(2,4,6,8,10,12,14))
dev.off()



# ~ GLMM 10 # Probability of hatching ~ ---------------------------------------

hatchVars <- c("outcome", "NO.BIRD", "exp", "year", "burrow")
nestdf_Hatch <- nestdf[hatchVars]
nestdf_Hatch$outcome <- ifelse(nestdf_Hatch$outcome == "HATCHED", 1, 0)
nestdf_Hatch <- subset(nestdf_Hatch, exp != "UNKNOWN")

## Get sample size
nrow(nestdf_Hatch); length(unique(nestdf_Hatch$burrow)) # 73 incubation attempts for 58 nests

## Construct binomial GLMM to analyse effect of neglect, experience, and year on the
## probability of hatching
glmmNEGFAIL <- glmer(outcome ~ NO.BIRD + exp + year + (1|burrow),
                     data = nestdf_Hatch,
                     family = "binomial",
                     glmerControl(optimizer = "bobyqa", 
                                  optCtrl = list(maxfun = 100000)))

drop1(glmmNEGFAIL, test = "Chisq", ddf = "lme4")
summary(glmmNEGFAIL)
glmmNEGFAIL_confint <- confint(glmmNEGFAIL, method = "Wald")


# ~~ FIGURE 3 ~~ ----------------------------------------------------------

## Plot probability of hatching as a function of days of neglect 
glmmNEGFAIL_pDat <- plot_model(glmmNEGFAIL, type = "pred")
glmmNEGFAIL_pDat <- data.frame(glmmNEGFAIL_pDat$NO.BIRD$data)

png("Figures/FIGURE3_hatching_neglect.png", height = 6, width = 9, units = 'in',res = 300)
ggplot(aes(x = x, y = predicted), data = glmmNEGFAIL_pDat) +
  geom_line(size = 1, col = "blue") +
  geom_count(aes(x = NO.BIRD, y = outcome), alpha = 0.4, data = nestdf_Hatch) +
  scale_size(name = expression(paste(italic("n"), " nests")), breaks = c(1,5,10,20,30)) +
  geom_ribbon(aes(ymin = conf.low, ymax = conf.high), data= glmmNEGFAIL_pDat, alpha = 0.2) +
  theme_bw(base_size = 15) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        axis.text.x = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
        axis.title.y = element_text(colour = "black", size = 17, margin = margin(t = 0, r = 12, b = 0, l = 0)),
        axis.text.y = element_text(colour = "black", angle = 0, hjust = 0.5, size=17),
        axis.title.x = element_text(colour = "black", size = 17,margin = margin(t = 12, r = 0, b =0, l = 0)),
        legend.position = c(0.93,0.82)) +
  labs(y = "Probability of hatching", x = "Total days of neglect")
dev.off()


# FORAGING BEHAVIOUR ------------------------------------------------------

# ~ GAM 11 (poisson) # Foraging trip duration ~ -----------------------------------------

## Extract relevant variables for model
forageVars <- c("foragetrip", "start_mass", "partnerincoming", "sex", "DSL", "year", 
                "ring", "burrow")
foragedf_Dur <- foragedf[forageVars]
foragedf_Dur <- na.omit(foragedf_Dur)

## Get sample size
nrow(foragedf_Dur); length(unique(foragedf_Dur$burrow)); length(unique(foragedf_Dur$ring)) # 597 trips for 74 nests and 154 individuals

## Construct Poisson GAMs to test effect of interacting start mass and partner mass, sex and year on foraging trip duration
gamFOR <- gam(foragetrip ~ te(start_mass, partnerincoming) + sex + year + s(ring, bs = "re"), family = poisson, data = foragedf_Dur, select = T)
gamFOR_noYear <- gam(foragetrip ~ te(start_mass, partnerincoming) + sex + s(ring, bs = "re"), family = poisson, data = foragedf_Dur)
gamFOR_noSex <- gam(foragetrip ~ te(start_mass, partnerincoming) + year + s(ring, bs = "re"), family = poisson, data = foragedf_Dur)

summary(gamFOR)
anova(gamFOR, gamFOR_noYear, test = "Chisq")
anova(gamFOR, gamFOR_noSex, test = "Chisq")

gam_summary <- summary(gamFOR)
coefs <- gam_summary$p.coeff
exp(gam_summary$p.coeff)

glmmFOR_cc <- confint(glmmFOR, method = "Wald")


# ~~ FIGURE 4 ~~ ----------------------------------------------------------

## Perspective plots of effect of start mass and partner mass on foraging trip duration

png("Figures/FIGURE4_foragetrip_startmass-partnermass.png", height = 6, width = 9, units = 'in',res = 300)
myvis.gam(x = gamFOR, view = c("start_mass", "partnerincoming"), plot.type = "persp", theta = 35, phi = 15, color = "jet",
          xlab = "Mass at start of foraging trip (g)", ylab = "Partner mass (g)", 
          ticktype = "detailed", type = "response")
dev.off()

png("Figures/FIGURE4_CIs_foragetrip_startmass-partnermass.png", height = 6, width = 9, units = 'in',res = 300)
vis.gam(x = gamFOR, view = c("start_mass", "partnerincoming"), plot.type = "persp", theta = 35, phi = 15,
          xlab = "Mass at start of foraging trip (g)", ylab = "Partner mass (g)", type = "response",
          se = 2, ticktype = "detailed", zlab = "Foraging trip duration (days)")
dev.off() 

# ~ ## LMM 12 ## Percentage daily mass gain ~ ---------------------------------

## Summarise daily foraging activity for each trip
daily_Trip <- data.frame (daily %>%
                         dplyr::group_by(ring, burrow, year, stintID, phase) %>%
                         dplyr::summarise(forage = median(prop_forage, na.rm =T), 
                                          flight = median(prop_flight, na.rm =T), 
                                          rest = median(prop_rest, na.rm =T),
                                          start_mass = start_mass[1],
                                          end_mass = end_mass[1],
                                          foragetrip = foragetrip[1],
                                          forageGain = forageGain[1],
                                          neglect = neglect[1],
                                          sex = sex[1])   )


## Isolate foraging section of trip
daily_Trip <- subset(daily_Trip, phase == "forage")
daily_Trip$phase <- NULL
daily_Trip$gainPer <- (daily_Trip$forageGain/daily_Trip$start_mass)*100
daily_Trip$year <- as.factor(daily_Trip$year)
daily_Trip$sex <- as.factor(daily_Trip$sex)

## Get sample sizes
nrow(daily_Trip); length(unique(daily_Trip$ring)); length(unique(daily_Trip$burrow)) # 161 trips for 54 individuals and 40 nests

## Mean daily foraging gains as a percentage and in grams
mean(daily_Trip$gainPer, na.rm = T); sd(daily_Trip$gainPer, na.rm = T) # 12.6 +/- 9.24 %
mean(daily_Trip$forageGain, na.rm = T); sd(daily_Trip$forageGain, na.rm = T) # 49.0 +/- 33.9%

## Construct model of effect of time spent foraging, sex, trip duration, start mass and year
## on percentage mass gains
gainLMM <- lmer(gainPer ~ forage + sex + foragetrip + 
                  start_mass + year + (1|ring),
                data = daily_Trip)

drop1(gainLMM, ddf = "lme4", test = "Chisq")
summary(gainLMM)
gainLMM_confint <- confint(gainLMM)

emmeans(gainLMM, "sex")
emmeans(gainLMM, "year")


# FIGURE S3 ---------------------------------------------------------------

## Plot effect of start mass and sex on percentage mass gains
gainLMM_pdat <- plot_model(gainLMM, type = "pred", terms = c("start_mass", "sex"))
gainLMM_pdat <- data.frame(gainLMM_pdat$data)

png("Figures/FIGURES3_massGain_mass.png", height = 6, width = 9, units = 'in',res = 300)
ggplot() +
  geom_point(data = daily_Trip, aes(x = start_mass, y = gainPer, col = sex)) +
  geom_ribbon(data = gainLMM_pdat, aes(x = x, ymin = conf.low, ymax = conf.high, group = group),
              alpha = 0.4, fill = "grey") +
  geom_line(data = gainLMM_pdat, aes(x = x, y = predicted, col = group), size = 1) +
  theme_custom() +
  theme(legend.position = c(0.9, 0.85)) + labs(col = "Sex") +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  labs(x = "Mass at start of foraging trip (g)", y = "Total foraging gain (% body mass)") +
  scale_colour_manual(values = c(female_col, male_col), labels = c("Female", "Male")) +
  geom_hline(yintercept = 0, col = "grey", lty = 2) 
dev.off()


# ~ ## GLMM 13/14/15 ## Proportional behaviours ~ -----------------------------

## Extract rows for which we have starting mass
daily_Prop <- subset(daily, !is.na(start_mass))

# Instead of phase, inc dummy variable as prop commuting time 
daily_Prop$propCom <- 2/daily_Prop$foragetrip

# Remove 1 day trips (not realistic) and those with values of 0 or 1 for behaviour
daily_Prop <- subset(daily_Prop, foragetrip > 1)
daily_Prop <- subset(daily_Prop, prop_forage > 0.01 & prop_forage < 1 &
                  prop_rest > 0.01 & prop_rest < 1 &
                  prop_flight > 0.01 & prop_flight < 0.9)

## Get sample size
nrow(daily_Prop); length(unique(daily_Prop$ring)) # 995 days for 54 rings

## Constuct beta mixed GLMM modls to examine effect of start mass, sex, trip duration,
## proportion commuting time, and year on proportion of time spent in each behaviour

# FORAGING # 
glmmFORBEH <- glmmTMB(prop_forage ~ start_mass + sex + foragetrip +
                        propCom + year + (1|ring), 
                      data = daily_Prop, beta_family())

drop1(glmmFORBEH, test = "Chisq")
summary(glmmFORBEH)
glmmFORBEH_conf <- confint(glmmFORBEH)

emmeans(glmmFORBEH, "sex", type = "response")
emmeans(glmmFORBEH, "year", type = "response")

forcoefs <- summary(glmmFORBEH)$coefficients

logit2prob(forcoefs[1]$cond[1]+forcoefs[1]$cond[2]*400)-logit2prob(forcoefs[1]$cond[1]+forcoefs[1]$cond[2]*300)


# RESTING #
glmmRESBEH <- glmmTMB(prop_rest ~ start_mass + sex + foragetrip +
                        propCom + year + (1|ring), 
                      data = daily_Prop, beta_family())

drop1(glmmRESBEH, test = "Chisq")
summary(glmmRESBEH)
glmRESBEH_conf <- confint(glmmRESBEH)

emmeans(glmmRESBEH, "year", type = "response")

restcoefs <- summary(glmmRESBEH)$coefficients

logit2prob(restcoefs[1]$cond[1]+restcoefs[1]$cond[2]*400)-logit2prob(restcoefs[1]$cond[1]+restcoefs[1]$cond[2]*300)

# FLYING #
glmmFLIBEH <- glmmTMB(prop_flight ~ start_mass + sex + foragetrip +
                        propCom + year + (1|ring), 
                      data = modDF, beta_family())

drop1(glmmFLIBEH, test = "Chisq")
summary(glmmFLIBEH)
glmmFLIBEH_conf <- confint(glmmFLIBEH)

emmeans(glmmFLIBEH, "year", type = "response")


# ~~ FIGURE 5 ~~ ----------------------------------------------------------

## Plot proportion time spent in each behaviour as a function of start mass

# FORAGE #
effects_glmmFORBEH <- as.data.frame(effect(term = c("start_mass"), mod = glmmFORBEH))

foragePlot <- ggplot() + 
  geom_point(aes(x = start_mass, y = prop_forage), data = daily_Prop) +
  geom_ribbon(data = effects_glmmFORBEH, aes(x = start_mass, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey") +
  geom_line(data = effects_glmmFORBEH, aes(x = start_mass, y = fit), size = 1, col = "orange") +
  theme_custom() +
  labs(x = "Mass at start of trip (g)", y = "Proportion of time spent in behaviour") +
  ggtitle("Foraging") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,1))  


# REST # 
effects_glmmRESBEH <- as.data.frame(effect(term = c("start_mass"), mod = glmmRESBEH))

restPlot <- ggplot() + 
  geom_point(aes(x = start_mass, y = prop_rest), data = daily_Prop) +
  geom_ribbon(data = effects_glmmRESBEH, aes(x = start_mass, ymin = lower, ymax = upper),
              alpha = 0.5, fill = "grey") +
  geom_line(data = effects_glmmRESBEH, aes(x = start_mass, y = fit), size = 1, col = "orange") +
  theme_custom() +
  labs(x = "Mass at start of trip (g)", y = "") + 
  ggtitle("Resting") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,1))  

# FLIGHT #
flightPlot <- ggplot(aes(x = start_mass, y = prop_flight), data = modDF) + 
  geom_point() +
  theme_tash() +
  labs(x = "Mass at start of trip (g)", y = "") +
  ggtitle("Flight") + theme(plot.title = element_text(hjust = 0.5)) +
  scale_y_continuous(breaks = c(0, 0.2, 0.4, 0.6, 0.8, 1), limits = c(0,1))  

png("Figures/FIGURE5_propBehviour_startMass.png", height = 6, width = 15, units = 'in',res = 300)
ggarrange(foragePlot, restPlot, flightPlot, nrow = 1, ncol = 3)
dev.off()








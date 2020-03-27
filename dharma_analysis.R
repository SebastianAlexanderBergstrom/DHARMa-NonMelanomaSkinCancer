library(DHARMa)
library(curl)
library(MASS)
library(glmmTMB)
library(splines)
source("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/nonMelanomaSkinCancerSweden/master/dataFormating.R")

# Read data ---------------------------------------------------------------

halland_population <- read.csv(curl("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/DHARMa-NonMelanomaSkinCancer/master/halland_population.csv"),
                               header=F,
                               sep=",",
                               encoding="UTF-8")
counts <- formatDF(read.csv(curl("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/nonMelanomaSkinCancerSweden/master/counts.csv"),
                            header=F,sep=";",skip=2,encoding="UTF-8"))
halland_data<-read.csv(curl("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/DHARMa-NonMelanomaSkinCancer/master/Halland-alder.csv"),sep=";")

halland_pop <- read.csv(curl("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/DHARMa-NonMelanomaSkinCancer/master/Hallandfolk.csv"),sep=";")

swedish_population <- read.csv(curl("https://raw.githubusercontent.com/SebastianAlexanderBergstrom/nonMelanomaSkinCancerSweden/master/swedishPopulation.csv"),
                            header=F,sep=";",skip=2,encoding="UTF-8")


# Define functions --------------------------------------------------------

format_population_data<-function(population_Vector,data_frame){
  # input: "populationVector" is a vector of the population in Sweden during the years 1970-2015, dataFrame
  # is the data.frame for which we want to add the population as a variable.
  
  # output: Returns a vector containing the population of Sweden formatted so that it can be added to
  # our data.frame object of interest.
  
  # n = number of times we should replicate every element in populationVector so all the data points for a given
  # year are matched with the correct element in populationVector.
  n <- nrow(data_frame[which(data_frame$years == data_frame$years[1]),])
  m <- nrow(data_frame)
  res_vec <- c()
  
  for(i in 1:46){ # fix hard coding
    res_vec <- append(res_vec,rep(population_Vector[i],n/2))
  }
  res_vec <- rep(res_vec,2)
  return(res_vec)
}

testAgeGrouping <- function(dataFrame,ageGrouping){
  # input: dataFrame is the data.frame object we want to change, ageGrouping is a vector containing
  # the new age grouping we want to try out.
  
  # output: Returns a data frame with the new age grouping that was given as a parameter.
  
  age <- rep(as.factor(ageGrouping),length(dataFrame$age)/18)
  df <- data.frame(dataFrame$cases,dataFrame$years,dataFrame$sexInd,dataFrame$counties,age,dataFrame$pop)
  colnames(df) <- c("cases","years","sexInd","counties","age","pop")
  return(df)
}

testCountyGrouping <- function(dataFrame,countyGrouping){
  # input: "countyGrouping" is list containing 2 vectors of the same length, one containing the names
  # of counties we wish to collapse as one factor level and the other one containing the name for the new
  # level. The new name is matched with the old one. "dataFrame" is the data.frame object whose factor
  # levels for the variable "counties" we want to change.
  
  # output: Returns a new data.frame object with changed factor levels for the variable "counties"
  
  oldNames <- countyGrouping[1]
  newNames <- countyGrouping[2]
  
  df <- data.frame(dataFrame$cases,dataFrame$sexInd,dataFrame$years,dataFrame$age,dataFrame$counties,dataFrame$pop)
  colnames(df) <- c("cases","sexInd","years","age","counties","pop")
  n <- length(unlist(list(oldNames)))
  
  for(i in 1:n){
    levels(df$counties)[levels(df$counties) == unlist(list(oldNames))[i]] <- unlist(list(newNames))[i]
  }
  return(df)
}


# Halland -----------------------------------------------------------------

numerical_values <- halland_population[4:length(halland_population)]
numerical_values <- numerical_values[2:nrow(numerical_values),2:ncol(numerical_values)]
numerical_values <- data.frame(colSums(numerical_values))
colnames(numerical_values) <- c("Antal")
population <- numerical_values$Antal
population <- rev(population)

counts_halland <- counts[which(counts$counties == "Halland"),]
counts_halland$pop <- format_population_data(population,counts_halland)

poisson_model_halland <- glm(cases ~ years + age + sexInd + offset(log(pop)),family=poisson ,data=counts_halland)
negbin_model_halland <- glm.nb(cases ~ years + sexInd + age + offset(log(pop)), data=counts_halland)

poisson_halland_simulation <- simulateResiduals(fittedModel = poisson_model_halland)
negbin_halland_simulation <- simulateResiduals(fittedModel = negbin_model_halland)

plot(poisson_halland_simulation)
plot(negbin_halland_simulation)
# poisson model looks pretty bad, looks like zero-inflation in Hartig's example
testZeroInflation(poisson_halland_simulation) # p-value > 0.05, don't reject H_0
# might also be due to overdispersion, check that
testDispersion(poisson_halland_simulation)
# reject H_0, seems like we have overdispersion

testUniformity(simulationOutput = negbin_halland_simulation)

# ignore poisson, focus on negbin a bit

testZeroInflation(negbin_halland_simulation)
testDispersion(negbin_halland_simulation)
# doesn't look very bad

colnames(halland_data)<-c("Ar","Grupp","M","K")

nr_of_years <- nrow(halland_data)
age_coding <- as.factor(c(2,7,12,17,22,27,32,37,42,47,52,57,62,67,72,77,82,85))
years <- rep(halland_data$Ar,2)
sex <- c(rep(1,nrow(halland_data)),rep(0,nrow(halland_data)))
group <- rep(age_coding,nrow(halland_data)/45)
halland_counts <- c(halland_data$M,halland_data$K)


halland_pop <- halland_pop[5:ncol(halland_pop)]
halland_pop <- colSums(halland_pop)
halland_pop <- rep(log(halland_pop),18)
halland_pop <- rep(halland_pop,2)

age_factors <- as.factor(c(rep(12,5),rep(32,3),42,47,52,57,62,67,72,rep(79.5,2),85))
age_grouping <- factor(rep(age_factors,nrow(halland_data)/45))
halland_cancer_data <- data.frame(halland_counts,age_grouping,years,sex,halland_pop)

negative_binomial_model_thesis <-glm.nb(halland_counts~age_grouping + ns(years,df=9)+ sex + age_grouping:sex + years:sex + years:age_grouping+offset(halland_pop),data=halland_cancer_data)

simulated_residuals <- simulateResiduals(fittedModel = negative_binomial_model_thesis)
plot(simulated_residuals)

# no significant deviation for outlier nor Kolmogorov-Smirnov test. Residual plots look decent,
# so you can't find evidence of any serious wrong-doing in your thesis work.

# All counties ------------------------------------------------------------

lst<-c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49",
       "50-54","55-59","60-64","65-69","70-74","75-59","80-84","85+")
age <- rep(factor(lst,levels=lst),nrow(counts)/length(counts$age))
counts$age <- age

# Change the "sex" variable to an indicator variable with the value 1 for men and 0 for women
counts$sexInd <- gsub("Men",1,counts$sexInd)
counts$sexInd <- gsub("Women",0,counts$sexInd)
counts$sexInd <- as.numeric(counts$sexInd)

# "Sweden" is removed as a factor, "Stockholm" is chosen as refernence level for "counties"
counts <- subset(counts,counties != "Sweden")
counts$counties <- factor(counts$counties)
counts$counties <- relevel(counts$counties,ref="Stockholm")

# We use the youngest age group as reference level.
counts$age <- relevel(counts$age,ref="0-4")

# We center the variable "years"
counts$years <- counts$years - rep(mean(1970:2015),nrow(counts))


numerical_values <- swedish_population[4:length(swedish_population)]
numerical_values <- numerical_values[2:nrow(numerical_values),2:ncol(numerical_values)]
numerical_values <- data.frame(colSums(numerical_values))
colnames(numerical_values) <- c("Antal")
population <- numerical_values$Antal
population <- rev(population)

counts$pop <- format_population_data(population,counts) # WRONG, YOU'RE USING THE HALLAND-DATA FOR POPULATION HERE. FIX THIS

poisson_model <- glm(cases~ years + age + sexInd + counties + offset(log(pop)),family=poisson,data=counts)
negbin_model <- glm.nb(cases ~ years + age + sexInd + counties + offset(log(pop)),data=counts)

poisson_simulation <- simulateResiduals(fittedModel = poisson_model)
negbin_simulation <- simulateResiduals(fittedModel = negbin_model)

plot(poisson_simulation)
plot(negbin_simulation)
# significant deviatitons in both, negbin looks slightly better

# investigate negbin more closely
testUniformity(negbin_simulation)
# can't see anything helpful from here

# can be zero-inflation, have very many zeros in the data
nrow(counts[which(counts$cases == 0),])/nrow(counts)
testZeroInflation(negbin_simulation)
# indicates we have zero inflation

# try a zero-inflated models

zero_inflated_poisson <- glmmTMB(cases ~ years + age + sexInd + counties + offset(log(pop)), ziformula = ~1 ,
                                 family = "poisson",
                                 data = counts,
                                 control = glmmTMBControl(optCtrl = list(iter.max=800,eval.max=800)))
poisson_zero_sim <- simulateResiduals(fittedModel = zero_inflated_poisson)
plot(poisson_zero_sim) # significant deviation
testZeroInflation(poisson_zero_sim) # still significant deviation

# try negbin 
zero_inflated_negbin_no_offset <- glmmTMB(cases ~ years + age + sexInd + counties, ziformula = ~1 ,
                                          family = nbinom2,
                                          data = counts,
                                          control = glmmTMBControl(optCtrl = list(iter.max=500,eval.max=500)))
negbin_zero_sim_no_offset <- simulateResiduals(fittedModel = zero_inflated_negbin_no_offset)
plot(negbin_zero_sim_no_offset)
testZeroInflation(negbin_zero_sim_no_offset)

# the zero inflation might be caused by the younger age gropus, see what happens if these aren't included. THis is mostly from a curiosity standpoints

counts_elderly <- counts[which(!counts$age %in% c("0-4","5-9","10-14","15-19","20-24","25-29","30-34","35-39","40-44","45-49")),]

poisson_elderly <- glm(cases ~ years + age + sexInd + counties + offset(log(pop)), family = "poisson", data = counts_elderly)
poisson_elderly_sim <- simulateResiduals(fittedModel = poisson_elderly)
plot(poisson_elderly_sim) # significant deviations

negbin_elderly <- glm.nb(cases~ years + age +sexInd+counties + offset(log(pop)), data=counts_elderly)
negbin_elderly_sim <- simulateResiduals(fittedModel = negbin_elderly)
plot(negbin_elderly_sim) # no significant deviations

zero_inflated_poisson_elderly <- glmmTMB(cases ~ years + age + sexInd + counties + offset(log(pop)), ziformula = ~1 ,
                                         family = "poisson",
                                         data = counts_elderly,
                                         control = glmmTMBControl(optCtrl = list(iter.max=800,eval.max=800)))
poisson_zero_sim_elderly <- simulateResiduals(fittedModel = zero_inflated_poisson_elderly)
plot(poisson_zero_sim_elderly) # significant deviation
testZeroInflation(poisson_zero_sim_elderly)

zero_inflated_negbin_elderly <- glmmTMB(cases ~ years + age + sexInd + counties + offset(log(pop)), ziformula = ~1 ,
                                        family = nbinom2,
                                        data = counts_elderly,
                                        control = glmmTMBControl(optCtrl = list(iter.max=800,eval.max=800)))

negbin_zero_sim_elderly <- simulateResiduals(fittedModel = zero_inflated_negbin_elderly)
plot(negbin_zero_sim_elderly) # not significant deviation
testZeroInflation(negbin_zero_sim_elderly)

zero_inflated_negbin_elderly <- glmmTMB(cases ~ years + age + sexInd + counties + offset(log(pop)), ziformula = ~1 ,
                                        family = nbinom2,
                                        data = counts_elderly,
                                        control = glmmTMBControl(optCtrl = list(iter.max=800,eval.max=800)))

negbin_zero_sim_elderly <- simulateResiduals(fittedModel = zero_inflated_negbin_elderly)
plot(negbin_zero_sim_elderly) # not significant deviation
testZeroInflation(negbin_zero_sim_elderly)

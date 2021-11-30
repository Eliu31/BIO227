library(tidyverse)
library(cowplot)
library(foreach)

theme_set(theme_cowplot())

########BASIC PARAMETERS##########
num.years <- 249 #Simulation Years
num.patch <- 10 #Number of Patches
J <- 100 # Patch carrying capacity
initial.patch <- 1 #Number of Initially Occupied Patches
Np <- 0.5 #Fraction of Species 1 in initial condition
m <- 0.01 #Dispersal Probability

########DISTURBANCE PARAMETERS##########
Dep.ratio <- 0.95 #Percentage of each patch depleted in disturbance event
Patch.Dep.ratio <- .9 #Percentage of patches depleted in disturbance event
Distubance_events <- 100 #Total Number of Disturbance Events
Disturbance_interval <- 25 #Time between each disturbance event
Disturb_year<- seq(Disturbance_interval,Disturbance_interval*Distubance_events,Disturbance_interval) # Disturbance time Vector

Density_Dependent_Disturbance  <- F #Enable density-dependent disturbance that waits until all patches are saturated Enabling this option disables the two interval-based options above.
Density_Dependent_Disturbance_Thresh <- .99 #Saturation threshold for density-dependent disturbance

########FITNESS PARAMETERS##########
fit.ratio.m <- 1/10 #Dispersal Fitness Ratio
fit.ratio.avg <- vector(length=num.patch)
fit.ratio.avg[] <- c(1.2,1.2) #Reproductive Fitness ratio
freq.dep <- vector(length=num.patch)
freq.dep[] <- 0 # Frequency Dependence

data <- foreach(interval=c(10,25,50,75,100,125,250), .combine="rbind") %do% {

  foreach(simulation=seq(1,3), .combine="rbind") %do% {
    Disturbance_interval <- interval #Time between each disturbance event
    Disturb_year<- seq(Disturbance_interval,Disturbance_interval*Distubance_events,Disturbance_interval) # Disturbance time Vector


    P <- c(rep(J, initial.patch), rep(0, num.patch-initial.patch)) # Defines Patch Population Vector
    init.1 <- Np*J
    COM <- matrix(nrow=J, ncol=num.patch)
    COM[1:init.1,] <- 1; COM[(init.1+1):J,] <- 2
    COM[,P==0]<-0
    year <- 2

    # Set up a matrix to record the frequencies over time.
    freq.1.mat <- matrix(nrow = num.years, ncol = num.patch)
    freq.1.mat[1,] <- (init.1*P)/J^2 #Changed initial frequency calculation due to changes in J/init.1
    # Set up a matrix to record the population sizes over time.
    pop.sizes <- matrix(nrow = num.years, ncol=num.patch)
    pop.sizes[1,] <- P

    ## run simulation
    for (i in 1:(J*num.patch*(num.years-1))) {

      ## choose a patch where a death even will occur
      patch <- sample(1:num.patch,1)
      occupied_patch<-sum(COM!=0)/(J) #Calculates the number of all patches that are currently occupied.

      ## calculate Pr.1 if dispersal occurs
      if (runif(1) < m) {
        freq.1.meta <- sum(COM==1)/(J*occupied_patch) #Frequency calculation divides by the number of living organisms rather than the total number of spots.
        Pr.1 <- fit.ratio.m*freq.1.meta/(fit.ratio.m*freq.1.meta + (1-freq.1.meta))
        Bool<-1  # Added a Boolean marker to signify dispersal case
      } else {

        ## calculate Pr.1 if local reproduction (not dispersal)
        freq.1 <- sum(COM[,patch]==1)/(sum(COM[,patch]!=0)+10^-20); freq.2 <- 1 - freq.1 #Same as above; frequency calculation reflects patch occupancy rather than raw patch number.
        fit.ratio <- exp(freq.dep[patch]*(freq.1-0.5) + log(fit.ratio.avg[patch]))
        Pr.1 <-  fit.ratio*freq.1/(fit.ratio*freq.1 + freq.2)
        Bool<-0 # See above
      }
      ##This line scales local reproduction rate by the total patch population. This prevents empty patches from reproducing.
      ## Newly-colonized patches with few individuals also have a chance of skipping reproduction proportional to % missing population.
      if(Bool==0&&sum(COM[,patch]!=0)<(runif(1)*J)){}else{
        COM[ceiling(J*runif(1)),patch] <- sample(c(1,2), 1, prob=c(Pr.1,1-round(Pr.1, digits=5)))}

      ## record data


      if (i %% (J*num.patch) == 0) {
        freq.1.mat[year,] <- colSums(COM==1)/J
        for (j in 1:num.patch) {
          pop.sizes[year-1,j] <- sum(COM[,j]!=0)
        }
        year <- year + 1
        if(year%in%Disturb_year&&Density_Dependent_Disturbance==F){
          print("Disturb")
          patch_sample<-sample(c(1:num.patch), replace=F, size=(round(num.patch*Patch.Dep.ratio)))
          COM[sample(c(1:J), replace=F, size=(round(J*Dep.ratio))),patch_sample]<-0
        }
        if(sum(COM!=0)>(J*num.patch*Density_Dependent_Disturbance_Thresh)&&Density_Dependent_Disturbance==T){
          print("DENS_Disturb")
          patch_sample<-sample(c(1:num.patch), replace=F, size=(round(num.patch*Patch.Dep.ratio)))
          COM[sample(c(1:J), replace=F, size=(round(J*Dep.ratio))),patch_sample]<-0
        }
      }
    }
    
  

    #######TIDY AND OUTPUT DATA#########
    # Convert frequency matrix into a dataframe.
    freqs <- as.data.frame(freq.1.mat) %>%
      mutate(time=1:n())
    # Tidy the dataframe.
    freqs <- freqs %>%
      group_by(time) %>%
      gather(patch, freq1, -time)
    # Convert population sizes into a dataframe.
    data.pop.sizes <- as.data.frame(pop.sizes) %>%
      mutate(time=1:n()+1)
    # Tidy the dataframe.
    data.pop.sizes <- data.pop.sizes %>%
      group_by(time) %>%
      gather(patch, populationSize, -time)
    # Combine the dataframes on species frequencies and population sizes.
    data.com <- left_join(freqs, data.pop.sizes, by=c("time","patch"))
    data.com %>% mutate(simulation=simulation, interval=interval)
  }

}


# Plot the overall abundance of species 1 across all patches.
data %>%
  group_by(time, interval, simulation) %>%
  mutate(species1=round(freq1*J),
         species2=populationSize-species1) %>%
  summarize(totalfreq1=sum(species1)/sum(populationSize)) %>%
  ggplot() +
  geom_line(aes(x=time, y=totalfreq1, group=factor(simulation))) +
  facet_wrap(~interval) +
  ylim(0,1)

# Calculate the first year in which species 1 reaches fixation
# across all patches for each simulation.
data %>%
  group_by(time, interval, simulation) %>%
  mutate(species1=round(freq1*J),
         species2=populationSize-species1) %>%
  summarize(totalfreq1=sum(species1)/sum(populationSize)) %>%
  ungroup() %>% group_by(interval, simulation) %>%
  filter(totalfreq1==1) %>% top_n(-1, time) %>% ungroup() %>%
  complete(totalfreq1, interval, simulation, fill=list(time=250)) %>%
  ggplot() +
  geom_point(aes(x=interval, y=time), alpha=0.5) +
  ylim(0,250)
# Plot the % of simulations in which species 1 fixes in 250 years.
data %>%
  group_by(time, interval, simulation) %>%
  mutate(numspecies1=freq1*populationSize) %>%
  summarize(totalfreq1=sum(numspecies1)/sum(populationSize)) %>%
  ungroup() %>% group_by(interval, simulation) %>%
  filter(totalfreq1==1) %>% top_n(-1, time) %>% ungroup() %>%
  complete(totalfreq1, interval, simulation, fill=list(time=250)) %>%
  group_by(interval) %>% summarize(numFixed=sum(time<250)) %>%
  ggplot() +
  geom_bar(aes(x=interval, y=numFixed), stat="identity") +
  ylim(0,max(data$simulation))

## graph the results
plot(1:num.years, freq.1.mat[,1], type="l", xlab="Time",
     ylab="Frequency of species 1", ylim=c(0,1))
for (i in 2:(num.patch)) {
  lines(1:num.years,freq.1.mat[,i], type="l", lty=2, ylim=c(0,1))}


# Convert frequency matrix into a dataframe.
freqs <- as.data.frame(freq.1.mat) %>%
  mutate(time=1:n())
# Tidy the dataframe.
freqs <- freqs %>%
  group_by(time) %>%
  gather(patch, freq1, -time)
# Plot species 1 frequencies over time.
freqs %>%
  ggplot() +
  geom_line(aes(x=time, y=freq1, group=factor(patch), color=factor(patch)))

# Convert population sizes into a dataframe.
data.pop.sizes <- as.data.frame(pop.sizes) %>%
  mutate(time=1:n()+1)
# Tidy the dataframe.
data.pop.sizes <- data.pop.sizes %>%
  group_by(time) %>%
  gather(patch, populationSize, -time)
# Plot population sizes over time.
data.pop.sizes %>%
  ggplot() +
  geom_line(aes(x=time, y=populationSize, group=factor(patch), color=factor(patch)))

# Plot the proportion of the population that consists of species 1 and 2.
# Combine the dataframes on species frequencies and population sizes.
data.com <- left_join(freqs, data.pop.sizes, by=c("time","patch"))
data.com %>%
  mutate(species1=round(freq1*J),
         species2=populationSize-species1) %>%
  dplyr::select(-freq1, -populationSize) %>%
  gather(species, populationSize, -patch, -time) %>%
  ggplot() +
  geom_area(aes(x=time, y=populationSize, fill=factor(species))) +
  facet_wrap(~patch)
# Plot the overall frequency of species 1 and 2 across all patches.
data.com %>%
  group_by(time) %>%
  mutate(numspecies1=freq1*populationSize) %>%
  summarize(totalfreq1=sum(numspecies1)/sum(populationSize)) %>%
  ggplot() +
  geom_line(aes(x=time, y=totalfreq1)) +
  ylim(0,1)
data.com %>%
  group_by(time) %>%
  mutate(species1=round(freq1*J),
         species2=populationSize-species1) %>%
  summarize(totalfreq1=sum(species1)/sum(populationSize)) %>%
  ggplot() +
  geom_line(aes(x=time, y=totalfreq1)) +
  ylim(0,1)

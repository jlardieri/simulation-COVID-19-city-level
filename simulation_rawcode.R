#load libraries
  
library(EpiModel)
library(extrafont)
library(animation)
library(tidyverse)
library(magrittr)
library(lubridate)
library(stringr)
library(tibble)
library(broom)
library(ggplot2)
library(knitr)
library(devtools)
library(DiagrammeR)
library(parallel)
library(foreach)
library(tictoc)
suppressMessages(library(EpiModel))
library(incidence)
library(earlyR)
library(future)

#monkey patch the EpiModel extension

source_files <- c("_icm.mod.init.seiqhrf.R", "_icm.mod.status.seiqhrf.R", 
                  "_icm.mod.vital.seiqhrf.R", "_icm.control.seiqhrf.R", "_icm.utils.seiqhrf.R", 
                  "_icm.saveout.seiqhrf.R", "_icm.icm.seiqhrf.R")

src_path <- paste0("./_posts/2020-03-18-modelling-the-effects-of-public-health-", 
                   "interventions-on-covid-19-transmission-part-2/")

gist_url <- "https://gist.github.com/timchurches/92073d0ea75cfbd387f91f7c6e624bd7"

local_source <- FALSE

for (source_file in source_files) {
  if (local_source) {
    source(paste(src_path, source_file, sep = ""))
  } else {
    source_gist(gist_url, filename = source_file)
  }
}

#read in data
setwd("")
#read in and remove note on first line of csv
df <- read.csv("HRR_Population.csv",header = TRUE,colClasses = "character")[-1,]
attach(df)
#include population columns only, will do own projection calculations
df <- df[1:9]

#clean hospital data
#create city, state fields
#replace any space in City name with an underscore
hospital <- df%>%separate(HRR,c("City","State"),sep=",\\s*")

#change character values to numeric
hospital[,3:10] <- lapply(hospital[,3:10],function(hospital){as.numeric(gsub(",", "", hospital))})

#total beds available
hospital <- hospital %>% mutate(Available_Beds_Total = (Available.Hospital.Beds + Available.ICU.Beds))
hospital <- hospital %>% mutate(Potential_Available_Beds_Total = Potentially.Available.Hospital.Beds. + Potentially.Available.ICU.Beds.)

# The EpiModel extension is named the Susceptible-Exposed-Infectious-Quarantined-Hospitalized-Recovered-Fatal model (`SEIQHRF`)
### City Simulation #1 : Boston, MA

boston <- hospital %>% filter(City=="Boston")
city1 <- boston
#scale
city1[3:12] <- city1[3:12]/100


# function to set-up and run the baseline simulations
# i kept many of the defaults here, but made some tweaks
simulate_bos <- function( #control
  type = "SEIQHRF", 
  nsteps = 366, #number of sim days, day 1 is initialization so 366 for full year
  nsims = 5, #take average of n number of sims
  ncores = 4, # cores to use
  prog.rand = FALSE, #false for random draws from weibull distribution, random binomial draws were computationally difficult
  rec.rand = FALSE, # random draws from weibull distribtuion at prog.rate for the progression from exposed + asymptomatic to infected + symptomatic
  fat.rand = TRUE, #random draws from weibull distribution at fat.rate.base for case fatality progression from hospitalization to fatal case
  quar.rand = FALSE, #random sample with sample faction at quar.rate for self-iso progression from infected + symptomatic to quarantine
  hosp.rand = FALSE, #random sample with a sample fraction at hosp.rate for progression from infected + symptomatic or quarantine to hospitalization
  disch.rand = TRUE, #random sample with sample fraction at disch.rate for progression from hospitalization to discharge
  infection.FUN = infection.seiqhrf.icm, #default function for infection process
  recovery.FUN = progress.seiqhrf.icm, #default function for recovery process
  departures.FUN = departures.seiqhrf.icm, #default function for non-COVID deaths
  arrivals.FUN = arrivals.icm, #default function for arrivals (births)
  get_prev.FUN = get_prev.seiqhrf.icm, #default utility function to store values
  
  # initial 
  s.num=city1$Adult.Population, #initial population of susceptible persons
  e.num=0, #initial number of exposed + asymptomatic
  i.num = 1, #initial number of infected + symptomatic
  q.num=0, # initial number of quarantined population
  h.num=0, #initial number of COVID needing hospitalization
  r.num = 0, #intial number of COVID recovered
  f.num = 0, #initial number of COVID deaths
  
  # params
  inf.prob.e = 0.025, #prob. of an exposed + asymp person infecting a susceptible person 
  act.rate.e = 10, # number of actions between exposed + asymp and a sus. person per day
  inf.prob.i = 0.05, # prob. of infected person infecting a sus. person 
  act.rate.i = 10, #number of actions between infected + symp person and a sus. person per day
  inf.prob.q = 0.02, #prob. of infected quarantine person infecting a sus. person
  act.rate.q = 2.5, # number of actions between quarantine person and a sus. person per day                    
  quar.rate = 1/30, # rate per day an infected person enters quarantine, very low rate where most people are not self-quarantine
  hosp.rate = 1/100, #rate per day an infected or quarantine person needs to be hospitalized 1% per day
  disch.rate = 1/15, #rate per day of discharges from hospital
  prog.rate = 1/10, #rate per day an exposed person becomes an infected person
  prog.dist.scale = 5,
  prog.dist.shape = 1.5,
  rec.rate = 1/20,# rate per day an infected person recovers
  rec.dist.scale = 35, #
  rec.dist.shape = 1.5,
  fat.rate.base = 1/50, #1/50 baseline COVID mortality rate per day for those needing hospitalization
  hosp.cap = city1$Available_Beds_Total, #number of available hospital beds
  fat.rate.overcap = 1/25, #mortality rate per day for people needing beds but who can't get one
  fat.tcoeff = 0.5, #time coeff. for increasing mortality rate as a hospitalized person spends more time in the hospital
  vital = TRUE, # enables arrivals and departures 
  a.rate = 0, # arrival rate into the population
  a.prop.e = 0.01,
  a.prop.i = 0.001,
  a.prop.q = 0.01,
  ds.rate = (7/365)/1000, #non-COVID deaths
  de.rate = (7/365)/1000, 
  di.rate = (7/365)/1000, #non-COVID exit rate for infected
  dq.rate = (7/365)/1000,
  dh.rate = (20/365)/1000,
  dr.rate = (7/365)/1000, #non-COVID exit rate for recovered 
  out="mean" #summary function for simulations - take average
) {
  control <- control.icm(type = type, 
                         nsteps = nsteps, 
                         nsims = nsims,
                         ncores = ncores,
                         prog.rand = prog.rand,
                         rec.rand = rec.rand,
                         infection.FUN = infection.FUN,
                         recovery.FUN = recovery.FUN,
                         arrivals.FUN = arrivals.FUN,
                         departures.FUN = departures.FUN,
                         get_prev.FUN = get_prev.FUN)
  
  init <- init.icm(s.num = s.num,
                   e.num = e.num,
                   i.num = i.num,
                   q.num = q.num,
                   h.num = h.num,
                   r.num = r.num,
                   f.num = f.num)
  
  param <-  param.icm(inf.prob.e = inf.prob.e, 
                      act.rate.e = act.rate.e,
                      inf.prob.i = inf.prob.i, 
                      act.rate.i = act.rate.i,
                      inf.prob.q = inf.prob.q, 
                      act.rate.q = act.rate.q,                    
                      quar.rate = quar.rate,
                      hosp.rate = hosp.rate,
                      disch.rate = disch.rate,
                      prog.rate = prog.rate,
                      prog.dist.scale = prog.dist.scale,
                      prog.dist.shape = prog.dist.shape,
                      rec.rate = rec.rate,
                      rec.dist.scale = rec.dist.scale,
                      rec.dist.shape = rec.dist.shape,
                      fat.rate.base = fat.rate.base,
                      hosp.cap = hosp.cap,
                      fat.rate.overcap = fat.rate.overcap,
                      fat.tcoeff = fat.tcoeff,
                      vital = vital,
                      a.rate = a.rate, 
                      a.prop.e = a.prop.e,
                      a.prop.i = a.prop.i,
                      a.prop.q = a.prop.q,
                      ds.rate = ds.rate, 
                      de.rate = de.rate, 
                      di.rate = di.rate,
                      dq.rate = dq.rate,
                      dh.rate = dh.rate,
                      dr.rate = dr.rate)
  
  sim <- icm.seiqhrf(param, init, control)
  sim_df <- as.data.frame(sim, out=out)
  
  return(list(sim=sim, df=sim_df))
}
set.seed(123)
baseline_sim <- simulate_bos(ncores=4)


#Check the duration distributions from the baseline simulation
# function from developer that allows user to extract times and look at duration periods 
get_times <- function(simulate_results) {
  
  sim <- simulate_results$sim
  
  for (s in 1:sim$control$nsims) {
    if (s == 1) {
      times <- sim$times[[paste("sim", s, sep = "")]]
      times <- times %>% mutate(s = s)
    } else {
      times <- times %>% bind_rows(sim$times[[paste("sim", 
                                                    s, sep = "")]] %>% mutate(s = s))
    }
  }
  
  times <- times %>% mutate(infTime = ifelse(infTime < 0, -5, 
                                             infTime), expTime = ifelse(expTime < 0, -5, expTime)) %>% 
    mutate(incubation_period = infTime - expTime, illness_duration = recovTime - 
             expTime, illness_duration_hosp = dischTime - expTime, 
           hosp_los = dischTime - hospTime, quarantine_delay = quarTime - 
             infTime, survival_time = fatTime - infTime) %>% 
    select(s, incubation_period, quarantine_delay, illness_duration, 
           illness_duration_hosp, hosp_los, survival_time) %>% 
    pivot_longer(-s, names_to = "period_type", values_to = "duration") %>% 
    mutate(period_type = factor(period_type, levels = c("incubation_period", 
                                                        "quarantine_delay", "illness_duration", "illness_duration_hosp", 
                                                        "hosp_los", "survival_time"), labels = c("Incubation period", 
                                                                                                 "Delay entering isolation", "Illness duration", "Illness (hosp)", 
                                                                                                 "Care required (hosp)", "Survival time of case fatalities"), 
                                ordered = TRUE))
  return(times)
}
times <- get_times(baseline_sim)

times %>% filter(duration <= 30) %>% ggplot(aes(x = duration)) +
  geom_bar(fill="darkblue") + facet_grid(period_type ~ ., scales = "free_y") +theme(strip.text.y = element_text(size=16))+
  labs(title = "Boston Baseline Simulation Distributions") + theme(axis.text = element_text(size=16))


baseline_plot_df <- baseline_sim$df %>% # use only the person columns
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  #first 120 days
  filter(time <= 120) %>% pivot_longer(-c(time), names_to = "Compartment", #individual compartmental model (ICM) 
                                       values_to = "count")

# colors to represent compartments
compcols <- c(s.num = "lightblue", e.num = "goldenrod1", i.num = "orangered1", 
              q.num = "cyan", h.num = "purple", r.num = "lightgreen", 
              f.num = "darkgray")
complabels <- c(s.num = "Susceptible", e.num = "Exposed", 
                i.num = "Infected", q.num = "Quarantine", h.num = "Hospitalized", 
                r.num = "Recovered", f.num = "Fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = `Compartment`)) + 
  geom_line(size = 1.2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                           labels = complabels) + theme(plot.background = element_rect(fill="gray")) + labs(title = "Boston Baseline Simulations (low rate of self-quarantine, social distance)", 
                                                                                                                                            x = "Days Since Start of COVID-19", y = "Persons (scaled 100)")


baseline_plot_df %>% filter(Compartment %in% c("e.num", "i.num", 
                                               "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, y = count, 
                                                                                          colour = Compartment)) + geom_line(size = 1.2, alpha = 0.7) + 
  scale_colour_manual(values = compcols, labels = complabels) + 
  theme(plot.background = element_rect(fill="gray")) + labs(title = "BOS Baseline Simulations (low rate of self-quarantine, social distance)", x = "Days Since Start of COVID-19", 
                                                            y = "Persons (scaled 100)")
# incidence rate
case_counts <- baseline_sim$df %>% select(time, se.flow)
# uncount them
case_dates <- case_counts %>% uncount(se.flow) %>% 
  pull(time)
# convert to an incidence object
case_all <- case_dates %>% incidence(.)
# find the peak of the epidemic curve
peak_of_epidemic_curve <- find_peak(case_all)
peak_of_epidemic_curve
# plot the case curve

plot(case_all,xlab="Time",ylab = "Daily Case Count",color = "darkblue")

### Intervention #1 :
## Majority of `Infected` self-isolate. We can gradually increase the self-quarantine rate from 1/30 to 7/10 over two weeks. 

quar_rate_inc <- function(t) {
ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
(0.7 - 0.0333)/15, 0.70))
}
inc_quar_rate_sim <- simulate_bos(quar.rate = quar_rate_inc(1:366))

inc_quar_rate_sim_plot_df <- inc_quar_rate_sim$df %>%
# use only the persons columns
select(time, s.num, e.num, i.num, q.num, 
h.num, r.num, f.num) %>%
# examine only the first 120 days
filter(time <= 120) %>%
pivot_longer(-c(time),
names_to="Compartment",
values_to="count")
inc_quar_rate_sim_plot_df %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) 

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_quar_rate_sim_plot_df %>%
mutate(experiment="Gradual inc. self-isolation")) %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
facet_grid(experiment ~ .) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray")) +
labs(title="BOS Simulation of baseline vs. gradual increase in self-isolation",
x="Days Since Start of COVID-19",
y="Persons (scaled 100)")

inc_quar_rate_sim_plot_df %>% filter(Compartment %in% c("e.num", "i.num"
)) %>% ggplot(aes(x = time, y = count, 
colour = Compartment)) + geom_line(size = 1.2, alpha = 0.7) + 
scale_colour_manual(values = compcols, labels = complabels) + 
theme(plot.background = element_rect(fill="gray")) + labs(title = "BOS Simulation of baseline vs. gradual increase in self-isolation)", x = "Days Since Start of COVID-19", 
y = "Persons (scaled 100)")
inc_quar_rate_sim_plot_df


baseline_plot_df %>%
  mutate(experiment="Baseline") %>%
  bind_rows(inc_quar_rate_sim_plot_df %>%
              mutate(experiment="Gradual inc. self-isolation")) %>%
  filter(Compartment %in% c("h.num",
                            "f.num")) %>%
  ggplot(aes(x=time, y=count, colour=Compartment)) +
  geom_line(size=1.2, alpha=0.7) +
  #change hospital count per use case
  geom_hline(yintercept = city1$Available_Beds_Total, colour="red", alpha=1) + 
  facet_grid(experiment ~ .) +
  scale_colour_manual(values = compcols, labels=complabels) +
  theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs gradual self-isolation simulations",
     x="Days Since Start of COVID-19",
     y="Persons")

### Intervention #2 : 
## The *HRR dataset* includes potentially available bed capacity. The metric assumes that a potential bed becomes available if 50% of used beds could free up. Intervention #2 will take into account the gradual increase in self-isolation plus a gradual increase in bed capacity. 
quar_rate_inc <- function(t) {
  ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
                                  (0.7 - 0.0333)/15, 0.7))
}
quar_rate_inc_sim <- simulate_bos(quar.rate = quar_rate_inc(1:366)) 

hosp_cap_inc <- function(t) {
  ifelse(t < 15, city1$Available_Beds_Total, ifelse(t <= 45, city1$Available_Beds_Total + (t - 15) * (city1$Potential_Available_Beds_Total - 
                                                                                                        city1$Available_Beds_Total)/30, city1$Potential_Available_Beds_Total))
}

inc_hosp_cap_sim_plot_df <- inc_hosp_cap_sim$df %>%
  # use only the persons columns
  select(time, s.num, e.num, i.num, q.num, 
         h.num, r.num, f.num) %>%
  # examine only the first 120 days
  filter(time <= 120) %>%
  pivot_longer(-c(time),
               names_to="Compartment",
               values_to="count")
inc_hosp_cap_sim_plot_df %>%
  filter(Compartment %in% c("e.num","i.num",
                            "q.num","h.num",
                            "f.num")) 


baseline_plot_df %>%
  mutate(experiment="Baseline") %>%
  bind_rows(inc_hosp_cap_sim_plot_df %>%
              mutate(experiment="Gradual inc. self-quarantine, hospital capacity")) %>%
  filter(Compartment %in% c("e.num","i.num",
                            "q.num","h.num",
                            "f.num")) %>%
  ggplot(aes(x=time, y=count, colour=Compartment)) +
  geom_line(size=1.2, alpha=0.7) +
  facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
  scale_colour_manual(values = compcols, labels=complabels) +
  theme(plot.background = element_rect(fill="gray")) +
  labs(title="BOS Baseline vs. gradual inc. self-quarantine, hospital capacity",
       x="Days Since Start of COVID-19",
       y="Persons (scaled 100)")

baseline_plot_df %>%
  mutate(experiment="Baseline") %>%
  bind_rows(inc_hosp_cap_sim_plot_df %>%
              mutate(experiment="Gradual inc. self-quarantine, hospital capacity")) %>%
  filter(Compartment %in% c("h.num",
                            "f.num")) %>%
  ggplot(aes(x=time, y=count, colour=Compartment)) +
  geom_line(size=1.2, alpha=0.7) +
  #change hospital count per use case
  geom_hline(yintercept = city1$Potential_Available_Beds_Total, colour="red", alpha=1) + 
  facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
  scale_colour_manual(values = compcols, labels=complabels) +
  theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs Gradual inc. self-quarantine, hospital capacity",
     x="Days Since Start of COVID-19",
     y="Persons")

###Intervention #3 :
quar_rate_inc <- function(t) {
ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
(0.7 - 0.0333)/15, 0.7))
}
quar_rate_inc_sim <- simulate_bos(quar.rate = quar_rate_inc(1:366)) 

hosp_cap_inc <- function(t) {
ifelse(t < 15, city1$Available_Beds_Total, ifelse(t <= 45, city1$Available_Beds_Total + (t - 15) * (city1$Potential_Available_Beds_Total - 
city1$Available_Beds_Total)/30, city1$Potential_Available_Beds_Total))
}
social_dist_inc <- function(t) {
ifelse(t < 30, 10, ifelse(t <= 45, 10 - (t - 30) * (10 - 
5)/15, 5))
}

social_dist_inc_sim <- simulate_bos(act.rate.i = social_dist_inc(1:366), 
act.rate.e = social_dist_inc(1:366),quar.rate = quar_rate_inc(1:366),hosp.cap = hosp_cap_inc(1:366))

social_dist_inc_sim_plot_df <- social_dist_inc_sim$df %>%
# use only the persons columns
select(time, s.num, e.num, i.num, q.num, 
h.num, r.num, f.num) %>%
# examine only the first 120 days
filter(time <= 120) %>%
pivot_longer(-c(time),
names_to="Compartment",
values_to="count")
social_dist_inc_sim_plot_df %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) 



baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_hosp_cap_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray")) +
labs(title="BOS Baseline vs. gradual inc. self-quarantine, social dist, hospital capacity",
x="Days Since Start of COVID-19",
y="Persons (scaled 100)")

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(social_dist_inc_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
#change hospital count per use case
geom_hline(yintercept = city1$Potential_Available_Beds_Total, colour="red", alpha=0.5) + 
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs Gradual inc. hospital capacity",
x="Days Since Start of COVID-19",
y="Persons")


plot_df <- baseline_sim$df %>% select(time, s.num, e.num, i.num, 
q.num, h.num, r.num, f.num) %>% mutate(experiment = "#0 Baseline") %>% 
bind_rows(inc_quar_rate_sim$df %>% select(time, s.num, e.num, 
i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#1")) %>% 
bind_rows(inc_hosp_cap_sim$df %>% select(time, s.num, e.num, 
i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#2")) %>% 
bind_rows(social_dist_inc_sim$df %>% select(time, s.num, 
e.num, i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#3")) %>% 
filter(time <= 150) %>% pivot_longer(-c(time, experiment), 
names_to = "compartment", values_to = "count") %>% filter(compartment %in% 
c("e.num", "i.num", "q.num", "h.num", "f.num"))

plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
geom_line(size = 1.2, alpha = 0.7) + facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) + scale_colour_manual(values = compcols, labels = complabels) + 
theme(plot.background = element_rect(fill="gray")) + labs(title = "Baseline vs experiments", x = "Days since start of COVID-19", 
y = "Persons (scaled 100)")

plot_df %>% filter(compartment %in% c("h.num", "f.num")) %>% 
ggplot(aes(x = time, y = count, colour = compartment)) + 
geom_line(size = 1.2, alpha = 0.7) + facet_grid(experiment ~ 
.) + scale_colour_manual(values = compcols, labels = complabels) + 
theme(plot.background = element_rect(fill="gray")) + labs(title = "Baseline vs experiments", x = "Days since start of COVID-19", 
y = "Persons")

##City Simulation #2 : Atlanta, GA
atlanta <- hospital %>% filter(City =="Atlanta")
city2 <- atlanta
#scale
city2[3:12] <- city2[3:12]/100


# function to set-up and run the baseline simulations
# i kept many of the defaults here, but made some tweaks
simulate_atl <- function( #control
type = "SEIQHRF", 
nsteps = 366, #number of sim days, day 1 is initialization so 366 for full year
nsims = 5, #take average of n number of sims
ncores = 4, # cores to use
prog.rand = FALSE, #false for random draws from weibull distribution, random binomial draws were computationally difficult
rec.rand = FALSE, # random draws from weibull distribtuion at prog.rate for the progression from exposed + asymptomatic to infected + symptomatic
fat.rand = TRUE, #random draws from weibull distribution at fat.rate.base for case fatality progression from hospitalization to fatal case
quar.rand = FALSE, #random sample with sample faction at quar.rate for self-iso progression from infected + symptomatic to quarantine
hosp.rand = FALSE, #random sample with a sample fraction at hosp.rate for progression from infected + symptomatic or quarantine to hospitalization
disch.rand = TRUE, #random sample with sample fraction at disch.rate for progression from hospitalization to discharge
infection.FUN = infection.seiqhrf.icm, #default function for infection process
recovery.FUN = progress.seiqhrf.icm, #default function for recovery process
departures.FUN = departures.seiqhrf.icm, #default function for non-COVID deaths
arrivals.FUN = arrivals.icm, #default function for arrivals (births)
get_prev.FUN = get_prev.seiqhrf.icm, #default utility function to store values

# initial 
s.num=city2$Adult.Population, #initial population of susceptible persons
e.num=0, #initial number of exposed + asymptomatic
i.num = 1, #initial number of infected + symptomatic
q.num=0, # initial number of quarantined population
h.num=0, #initial number of COVID needing hospitalization
r.num = 0, #intial number of COVID recovered
f.num = 0, #initial number of COVID deaths

# params
inf.prob.e = 0.025, #prob. of an exposed + asymp person infecting a susceptible person 
act.rate.e = 10, # number of actions between exposed + asymp and a sus. person per day
inf.prob.i = 0.05, # prob. of infected person infecting a sus. person 
act.rate.i = 10, #number of actions between infected + symp person and a sus. person per day
inf.prob.q = 0.02, #prob. of infected quarantine person infecting a sus. person
act.rate.q = 2.5, # number of actions between quarantine person and a sus. person per day                    
quar.rate = 1/30, # rate per day an infected person enters quarantine, very low rate where most people are not self-quarantine
hosp.rate = 1/100, #rate per day an infected or quarantine person needs to be hospitalized 1% per day
disch.rate = 1/15, #rate per day of discharges from hospital
prog.rate = 1/10, #rate per day an exposed person becomes an infected person
prog.dist.scale = 5,
prog.dist.shape = 1.5,
rec.rate = 1/20,# rate per day an infected person recovers
rec.dist.scale = 35, #
rec.dist.shape = 1.5,
fat.rate.base = 1/50, #1/50 baseline COVID mortality rate 
hosp.cap = city2$Available_Beds_Total, #number of available hospital beds
fat.rate.overcap = 1/25, #mortality rate per day for people needing beds but who can't get one
fat.tcoeff = 0.5, #time coeff. for increasing mortality rate as a hospitalized person spends more time in the hospital
vital = TRUE, # enables arrivals and departures 
a.rate = 0, # arrival rate into the population
a.prop.e = 0.01,
a.prop.i = 0.001,
a.prop.q = 0.01,
ds.rate = (7/365)/1000, #non-COVID deaths
de.rate = (7/365)/1000, 
di.rate = (7/365)/1000, #non-COVID exit rate for infected
dq.rate = (7/365)/1000,
dh.rate = (20/365)/1000,
dr.rate = (7/365)/1000, #non-COVID exit rate for recovered 
out="mean" #summary function for simulations - take average
) {
  control <- control.icm(type = type, 
                         nsteps = nsteps, 
                         nsims = nsims,
                         ncores = ncores,
                         prog.rand = prog.rand,
                         rec.rand = rec.rand,
                         infection.FUN = infection.FUN,
                         recovery.FUN = recovery.FUN,
                         arrivals.FUN = arrivals.FUN,
                         departures.FUN = departures.FUN,
                         get_prev.FUN = get_prev.FUN)
  
  init <- init.icm(s.num = s.num,
                   e.num = e.num,
                   i.num = i.num,
                   q.num = q.num,
                   h.num = h.num,
                   r.num = r.num,
                   f.num = f.num)
  
  param <-  param.icm(inf.prob.e = inf.prob.e, 
                      act.rate.e = act.rate.e,
                      inf.prob.i = inf.prob.i, 
                      act.rate.i = act.rate.i,
                      inf.prob.q = inf.prob.q, 
                      act.rate.q = act.rate.q,                    
                      quar.rate = quar.rate,
                      hosp.rate = hosp.rate,
                      disch.rate = disch.rate,
                      prog.rate = prog.rate,
                      prog.dist.scale = prog.dist.scale,
                      prog.dist.shape = prog.dist.shape,
                      rec.rate = rec.rate,
                      rec.dist.scale = rec.dist.scale,
                      rec.dist.shape = rec.dist.shape,
                      fat.rate.base = fat.rate.base,
                      hosp.cap = hosp.cap,
                      fat.rate.overcap = fat.rate.overcap,
                      fat.tcoeff = fat.tcoeff,
                      vital = vital,
                      a.rate = a.rate, 
                      a.prop.e = a.prop.e,
                      a.prop.i = a.prop.i,
                      a.prop.q = a.prop.q,
                      ds.rate = ds.rate, 
                      de.rate = de.rate, 
                      di.rate = di.rate,
                      dq.rate = dq.rate,
                      dh.rate = dh.rate,
                      dr.rate = dr.rate)
  
  sim <- icm.seiqhrf(param, init, control)
  sim_df <- as.data.frame(sim, out=out)
  
  return(list(sim=sim, df=sim_df))
}
set.seed(123)
baseline_sim <- simulate_atl(ncores=4)


#Check the duration distributions from the baseline simulation
# function from developer that allows user to extract times and look at duration periods 
get_times <- function(simulate_results) {
  
  sim <- simulate_results$sim
  
  for (s in 1:sim$control$nsims) {
    if (s == 1) {
      times <- sim$times[[paste("sim", s, sep = "")]]
      times <- times %>% mutate(s = s)
    } else {
      times <- times %>% bind_rows(sim$times[[paste("sim", 
                                                    s, sep = "")]] %>% mutate(s = s))
    }
  }
  
  times <- times %>% mutate(infTime = ifelse(infTime < 0, -5, 
                                             infTime), expTime = ifelse(expTime < 0, -5, expTime)) %>% 
    mutate(incubation_period = infTime - expTime, illness_duration = recovTime - 
             expTime, illness_duration_hosp = dischTime - expTime, 
           hosp_los = dischTime - hospTime, quarantine_delay = quarTime - 
             infTime, survival_time = fatTime - infTime) %>% 
    select(s, incubation_period, quarantine_delay, illness_duration, 
           illness_duration_hosp, hosp_los, survival_time) %>% 
    pivot_longer(-s, names_to = "period_type", values_to = "duration") %>% 
    mutate(period_type = factor(period_type, levels = c("incubation_period", 
                                                        "quarantine_delay", "illness_duration", "illness_duration_hosp", 
                                                        "hosp_los", "survival_time"), labels = c("Incubation period", 
                                                                                                 "Delay entering isolation", "Illness duration", "Illness (hosp)", 
                                                                                                 "Care required (hosp)", "Survival time of case fatalities"), 
                                ordered = TRUE))
  return(times)
}
times <- get_times(baseline_sim)

times %>% filter(duration <= 30) %>% ggplot(aes(x = duration)) +
  geom_bar(fill="darkblue") + facet_grid(period_type ~ ., scales = "free_y") +theme(strip.text.y = element_text(size=16))+
  labs(title = "Atlanta Baseline Simulation Distributions") + theme(axis.text = element_text(size=16))

baseline_plot_df <- baseline_sim$df %>% # use only the person columns
  select(time, s.num, e.num, i.num, q.num, h.num, r.num, f.num) %>% 
  #first 120 days
  filter(time <= 120) %>% pivot_longer(-c(time), names_to = "Compartment", #individual compartmental model (ICM) 
                                       values_to = "count")

# colors to represent compartments
compcols <- c(s.num = "lightblue", e.num = "goldenrod1", i.num = "orangered1", 
              q.num = "cyan", h.num = "purple", r.num = "lightgreen", 
              f.num = "darkgray")
complabels <- c(s.num = "Susceptible", e.num = "Exposed", 
                i.num = "Infected", q.num = "Quarantine", h.num = "Hospitalized", 
                r.num = "Recovered", f.num = "Fatality")

baseline_plot_df %>% ggplot(aes(x = time, y = count, colour = `Compartment`)) + 
  geom_line(size = 1.2, alpha = 0.7) + scale_colour_manual(values = compcols, 
                                                           labels = complabels) + theme(plot.background = element_rect(fill="gray")) + labs(title = "Atlanta Baseline Simulations (low rate of self-quarantine, social distance)", 
                                                                                                                                            x = "Days Since Start of COVID-19", y = "Persons (scaled 100)")


baseline_plot_df %>% filter(Compartment %in% c("e.num", "i.num", 
                                               "q.num", "h.num", "f.num")) %>% ggplot(aes(x = time, y = count, 
                                                                                          colour = Compartment)) + geom_line(size = 1.2, alpha = 0.7) + 
  scale_colour_manual(values = compcols, labels = complabels) + 
  theme(plot.background = element_rect(fill="gray")) + labs(title = "Atlanta Baseline Simulations (low rate of self-quarantine, social distance)", x = "Days Since Start of COVID-19", 
                                                            y = "Persons (scaled 100)")

# incidence rate
case_counts <- baseline_sim$df %>% select(time, se.flow)
# uncount them
case_dates <- case_counts %>% uncount(se.flow) %>% 
  pull(time)
# convert to an incidence object
case_all <- case_dates %>% incidence(.)
# find the peak of the epidemic curve
peak_of_epidemic_curve <- find_peak(case_all)
peak_of_epidemic_curve
# plot the case curve

plot(case_all,xlab="Time",ylab = "Daily Case Count",color = "darkblue")

quar_rate_inc <- function(t) {
  ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
                                  (0.7 - 0.0333)/15, 0.7))
}
inc_quar_rate_sim <- simulate_atl(quar.rate = quar_rate_inc(1:366))

inc_quar_rate_sim_plot_df <- inc_quar_rate_sim$df %>%
  # use only the persons columns
  select(time, s.num, e.num, i.num, q.num, 
         h.num, r.num, f.num) %>%
  # examine only the first 120 days
  filter(time <= 120) %>%
  pivot_longer(-c(time),
               names_to="Compartment",
               values_to="count")
inc_quar_rate_sim_plot_df %>%
  filter(Compartment %in% c("e.num","i.num",
                            "q.num","h.num",
                            "f.num")) 

baseline_plot_df %>%
  mutate(experiment="Baseline") %>%
  bind_rows(inc_quar_rate_sim_plot_df %>%
              mutate(experiment="Gradual inc. self-isolation")) %>%
  filter(Compartment %in% c("e.num","i.num",
                            "q.num","h.num",
                            "f.num")) %>%
  ggplot(aes(x=time, y=count, colour=Compartment)) +
  geom_line(size=1.2, alpha=0.7) +
  facet_grid(experiment ~ .) +
  scale_colour_manual(values = compcols, labels=complabels) +
  theme(plot.background = element_rect(fill="gray")) +
  labs(title="ATL Simulation of baseline vs. gradual increase in self-isolation",
       x="Days Since Start of COVID-19",
       y="Persons (scaled 100)")

inc_quar_rate_sim_plot_df %>% filter(Compartment %in% c("e.num", "i.num"
)) %>% ggplot(aes(x = time, y = count, 
                  colour = Compartment)) + geom_line(size = 1.2, alpha = 0.7) + 
  scale_colour_manual(values = compcols, labels = complabels) + 
  theme(plot.background = element_rect(fill="gray")) + labs(title = "ATL Simulation of baseline vs. gradual increase in self-isolation)", x = "Days Since Start of COVID-19", 
                                                            y = "Persons (scaled 100)")

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_quar_rate_sim_plot_df %>%
mutate(experiment="Gradual inc. self-isolation")) %>%
filter(Compartment %in% c("h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
#change hospital count per use case
geom_hline(yintercept = city2$Available_Beds_Total, colour="red", alpha=1) + 
facet_grid(experiment ~ .) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs gradual self-isolation simulations",
x="Days Since Start of COVID-19",
y="Persons")

### Intervention #2 : 
## The *HRR dataset* includes potentially available bed capacity. The metric assumes that a potential bed becomes available if 50% of used beds could free up. Intervention #2 will take into account the gradual increase in self-isolation plus a gradual increase in bed capacity. 

quar_rate_inc <- function(t) {
ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
(0.7 - 0.0333)/15, 0.7))
}

hosp_cap_inc <- function(t) {
ifelse(t < 15, city2$Available_Beds_Total, ifelse(t <= 45, city2$Available_Beds_Total + (t - 15) * (city2$Potential_Available_Beds_Total - 
city2$Available_Beds_Total)/30, city2$Potential_Available_Beds_Total))
}
inc_hosp_cap_sim <- simulate_atl(quar.rate = quar_rate_inc(1:366),hosp.cap = hosp_cap_inc(1:366)) 

inc_hosp_cap_sim_plot_df <- inc_hosp_cap_sim$df %>%
# use only the persons columns
select(time, s.num, e.num, i.num, q.num, 
h.num, r.num, f.num) %>%
# examine only the first 120 days
filter(time <= 120) %>%
pivot_longer(-c(time),
names_to="Compartment",
values_to="count")
inc_hosp_cap_sim_plot_df %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) 

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_hosp_cap_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, hospital capacity")) %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray")) +
labs(title="ATL Baseline vs. gradual inc. self-quarantine, hospital capacity",
x="Days Since Start of COVID-19",
y="Persons (scaled 100)")

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_hosp_cap_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, hospital capacity")) %>%
filter(Compartment %in% c("h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
#change hospital count per use case
geom_hline(yintercept = city2$Potential_Available_Beds_Total, colour="red", alpha=1) + 
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs Gradual inc. self-quarantine, hospital capacity",
x="Days Since Start of COVID-19",
y="Persons")



###Intervention #3 :

quar_rate_inc <- function(t) {
ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
(0.7 - 0.0333)/15, 0.7))
}

hosp_cap_inc <- function(t) {
ifelse(t < 15, city2$Available_Beds_Total, ifelse(t <= 45, city2$Available_Beds_Total + (t - 15) * (city2$Potential_Available_Beds_Total - 
city2$Available_Beds_Total)/30, city2$Potential_Available_Beds_Total))
}
social_dist_inc <- function(t) {
ifelse(t < 15, 10, ifelse(t <=30, 10 - (t - 15) * (10 - 
7.5)/15, 7.5))
}

social_dist_inc_sim <- simulate_atl(act.rate.i = social_dist_inc(1:366), 
act.rate.e = social_dist_inc(1:366),quar.rate = quar_rate_inc(1:366),hosp.cap = hosp_cap_inc(1:366))

social_dist_inc_sim_plot_df <- social_dist_inc_sim$df %>%
# use only the persons columns
select(time, s.num, e.num, i.num, q.num, 
h.num, r.num, f.num) %>%
# examine only the first 120 days
filter(time <= 120) %>%
pivot_longer(-c(time),
names_to="Compartment",
values_to="count")
social_dist_inc_sim_plot_df %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) 

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_hosp_cap_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray")) +
labs(title="ATL Baseline vs. gradual inc. self-quarantine, social dist, hospital capacity",
x="Days Since Start of COVID-19",
y="Persons (scaled 100)")

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(social_dist_inc_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
#change hospital count per use case
geom_hline(yintercept = city2$Potential_Available_Beds_Total, colour="red", alpha=0.5) + 
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs Gradual inc. hospital capacity",
x="Days Since Start of COVID-19",
y="Persons")



###Intervention #4 :

quar_rate_inc <- function(t) {
ifelse(t < 15, 0.0333, ifelse(t <= 30, 0.0333 + (t - 15) * 
(0.7 - 0.0333)/15, 0.7))
}


hosp_cap_inc <- function(t) {
ifelse(t < 15, city2$Available_Beds_Total, ifelse(t <= 45, city2$Available_Beds_Total + (t - 15) * (city2$Potential_Available_Beds_Total - 
city2$Available_Beds_Total)/30, city2$Potential_Available_Beds_Total))
}

social_dist_inc2 <- function(t) {
ifelse(t < 15, 10, ifelse(t <= 30, 10 - (t - 15) * (10 - 
7.5)/15, ifelse(t<=60, 7.5- (t-30) * (7.5 - 10)/30,10)))
}


social_dist_inc_sim2 <- simulate_atl(act.rate.i = social_dist_inc2(1:366), 
act.rate.e = social_dist_inc2(1:366),quar.rate = quar_rate_inc(1:366),hosp.cap = hosp_cap_inc(1:366))

social_dist_inc_sim_plot_df <- social_dist_inc_sim2$df %>%
# use only the persons columns
select(time, s.num, e.num, i.num, q.num, 
h.num, r.num, f.num) %>%
# examine only the first 120 days
filter(time <= 120) %>%
pivot_longer(-c(time),
names_to="Compartment",
values_to="count")
social_dist_inc_sim_plot_df %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) 


baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(inc_hosp_cap_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("e.num","i.num",
"q.num","h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray")) +
labs(title="Baseline vs. gradual inc. self-quarantine, social dist, hospital capacity",
x="Days Since Start of COVID-19",
y="Persons (scaled 100)")

baseline_plot_df %>%
mutate(experiment="Baseline") %>%
bind_rows(social_dist_inc_sim_plot_df %>%
mutate(experiment="Gradual inc. self-quarantine, social dist, hospital capacity")) %>%
filter(Compartment %in% c("h.num",
"f.num")) %>%
ggplot(aes(x=time, y=count, colour=Compartment)) +
geom_line(size=1.2, alpha=0.7) +
#change hospital count per use case
geom_hline(yintercept = city2$Potential_Available_Beds_Total, colour="red", alpha=0.5) + 
facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) +
scale_colour_manual(values = compcols, labels=complabels) +
theme(plot.background = element_rect(fill="gray"))
labs(title="Baseline vs Gradual inc. hospital capacity",
x="Days Since Start of COVID-19",
y="Persons")


##Compare the Interventions
plot_df <- baseline_sim$df %>% select(time, s.num, e.num, i.num, 
q.num, h.num, r.num, f.num) %>% mutate(experiment = "#0 Baseline") %>% 
bind_rows(inc_quar_rate_sim$df %>% select(time, s.num, e.num, 
i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#1")) %>% 
bind_rows(inc_hosp_cap_sim$df %>% select(time, s.num, e.num, 
i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#2")) %>% 
bind_rows(social_dist_inc_sim$df %>% select(time, s.num, 
e.num, i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#3")) %>% 
bind_rows(social_dist_inc_sim2$df %>% select(time, s.num, 
e.num, i.num, q.num, h.num, r.num, f.num) %>% mutate(experiment = "#4")) %>% 
filter(time <= 150) %>% pivot_longer(-c(time, experiment), 
names_to = "compartment", values_to = "count") %>% filter(compartment %in% 
c("e.num", "i.num", "q.num", "h.num", "f.num"))

plot_df %>% ggplot(aes(x = time, y = count, colour = compartment)) + 
geom_line(size = 1.2, alpha = 0.7) + facet_grid(experiment ~ .,labeller = label_wrap_gen(width = 35,multi_line = TRUE)) + scale_colour_manual(values = compcols, labels = complabels) + 
theme(plot.background = element_rect(fill="gray")) + labs(title = "Baseline vs experiments", x = "Days since start of COVID-19", 
y = "Persons (scaled 100)")

plot_df %>% filter(compartment %in% c("h.num", "f.num")) %>% 
ggplot(aes(x = time, y = count, colour = compartment)) + 
geom_line(size = 1.2, alpha = 0.7) + facet_grid(experiment ~ 
.) + scale_colour_manual(values = compcols, labels = complabels) + 
theme(plot.background = element_rect(fill="gray")) + labs(title = "Baseline vs experiments", x = "Days since start of COVID-19", 
y = "Persons")


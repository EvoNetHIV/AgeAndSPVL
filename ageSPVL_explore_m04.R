
library(evonet)
options(error=browser) # go into debug mode on error
#--------------------------------------------------------------

initial_pop = 1000
mean_sqrt_age_diff = 0.3
meandeg = 0.7

param_list=list(
  model_name = "age_mixing04",
  nw_form_terms = "~edges+absdiff('sqrt_age') + offset(nodematch('role',diff=TRUE, keep=1:2))",
  target_stats = c(initial_pop*meandeg/2, mean_sqrt_age_diff*initial_pop*meandeg/2),
  #mean_sex_acts_day = 0.2,
  min_age = 18,
  max_age = 55,

  #age_dist = seq(50, 10, -10/9)/1110,
  initial_agedata_male = "linear_decrease",
  nw_coef_form = c(-Inf, -Inf),
  prob_sex_by_age	= TRUE,
  prob_sex_age_19	= 0.2,
  max_age_sex	= 55,
  
  nsims = 1,
  initial_pop = initial_pop,
  initial_infected = 100,
  n_steps = 365*20,
  popsumm_frequency=30,
  fast_edgelist=T,
  plot_nw=F
)

evoparams <- do.call(evonet_setup,param_list)
nw <- nw_setup(evoparams)
modules <- c(
  "aging",
  "testing",
  "treatment",
  "viral_update",
  "coital_acts",
  "transmission",
  "deaths",
  "births",
  "summary_module"
)

evomodel <- evorun(modules,evoparams,nw)
ageSPVL_m04 <- evomodel
save(ageSPVL_m04, file="ageSPVL_m04.rda")

#to specify name of output pdf and pathway, use arguments
# name="xyz.pdf",outpath="/path to folder"
evoplot(model=evomodel, name="ageSPVL_m04")

popatts1 <- ageSPVL_m01$pop[[1]]
plot(popatts1$age_infection[popatts1$Time_Inf>0],
     log(popatts1$SetPoint[popatts1$Time_Inf>0],10))
popatts4 <- ageSPVL_m04$pop[[1]]
plot(popatts4$age_infection[popatts4$Time_Inf>0],
     log(popatts4$SetPoint[popatts4$Time_Inf>0],10))

summary(glm(log(popatts1$SetPoint[popatts1$Time_Inf>0],10)~
              popatts1$age_infection[popatts1$Time_Inf>0]+
              popatts1$Time_Inf[popatts1$Time_Inf>0]))
summary(glm(log(popatts4$SetPoint[popatts4$Time_Inf>0],10)~
              popatts4$age_infection[popatts4$Time_Inf>0]+
              popatts4$Time_Inf[popatts4$Time_Inf>0]))
summary(glm(log(popatts1$SetPoint[popatts1$Time_Inf>0],10)~
              popatts1$age_infection[popatts1$Time_Inf>0]))
summary(glm(log(popatts4$SetPoint[popatts4$Time_Inf>0],10)~
              popatts4$age_infection[popatts4$Time_Inf>0]))

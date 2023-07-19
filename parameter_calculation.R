### PACKAGES ### ----
library(rstudioapi) # to set the working directory
library(car) # Anova function (type II anova)
library(tidyr) # for the separate function
library(reshape2) # for data.frame formatting
library(nlme) # for mixed linear models
#library(terra)
library(emdbook) # to use the lambertW function
library(stats) # to use the nls function
library(doBy)
library(cowplot) # to arrange graphs
library(ggrepel) # to use geom_repel (labels)

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
source('figure_settings.R')
path_data="Data/"
path_figure="Figures/"

### GENERAL VARIABLES #### ----
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
fresh_to_dry_litter=0.65 # Joly et la. (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
second_to_day=1/(24*3600)
hour_to_day=1/24
month_to_day=30
year_to_day=365
celsius_to_kelvin=273.15
kB=8.62e-5 # Bolzmann's constant
water_to_soil=1/0.2 # Heathman et al. (2003) volume of soil corresponding to a volume of water
cell_to_fresh=1e-9 # Makarieva et al. (2005) (fresh mass in mg)
molC_to_mgC=12.011e3 # conversion to mg !!!!!
molN_to_mgN=14.0067e3 # conversion to mg !!!!!
volume_to_surface=1/0.1 # 1x1x0.1 m parallelepiped
soil_mass_to_volume=1/(1.4e3) # mg m-3 Bowden et al. (2014) or 1.4 g.cm-3 Heathman et al. (2003)
ha_to_m2=1e4

allometric_parameter<-function(bodymass,a_x,x0,s_x,E_x,kB,temperature){
  return(a_x*x0*bodymass^s_x*exp(-E_x/(kB*temperature)))
}

functional_response<-function(prey_biomass,predator_bodymass,a_a,a0,s_a,E_a,a_h,h0,s_h,E_h,kB,temperature){
  a=allometric_parameter(predator_bodymass,a_a,a0,s_a,E_a,kB,temperature)
  h=allometric_parameter(predator_bodymass,a_h,h0,s_h,E_h,kB,temperature)
  FR=a*prey_biomass/(1+a*h*prey_biomass)
  return(FR)
}

### METABOLISM ### ----
# data from Johnston and Sibly (2018) (10.1038/s41559-018-0648-6)
joule_to_carbon=1/20.1*0.5363
# metabolic rate : J.h-1
# body mass : mg fresh weight
data_metabo<-read.csv(file=paste(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaMetabolicData.csv",sep=""),
                      header=TRUE)
names(data_metabo)<-c("study","trophic_group","taxon","metabolic_rate","bodymass","temperature")
  
ggplot(data=data_metabo[data_metabo$trophic_group=="Macrofauna",])+
  geom_point(aes(bodymass,metabolic_rate/bodymass,colour=taxon))+
  theme+theme(legend.title=element_blank())+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Fresh body mass (mg)"))+
  ylab(expression("Respiration rate (J hr"^-1*")"))

# original values before conversion (same results as Johnston)
data_metabo$temperature_arrhenius<--1/(celsius_to_kelvin+data_metabo$temperature)*1/kB
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~trophic_group+log(bodymass):trophic_group+temperature:trophic_group)
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")

# conversion
data_metabo$metabolic_rate<-data_metabo$metabolic_rate*joule_to_carbon/hour_to_day
data_metabo$bodymass<-data_metabo$bodymass*fresh_to_dry*dry_to_C
data_metabo$metabolic_rate<-data_metabo$metabolic_rate/data_metabo$bodymass # conversion into mass metabolic rate
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~trophic_group+log(bodymass):trophic_group+temperature_arrhenius:trophic_group)
par(mfrow=c(2,2))
plot(model)
params<-summary(model)
params
Anova(model,type=2,test.statisticf="F")
params<-(params$coefficients)[,1]
names(params)[1]="trophic_groupMacrofauna"
params[2]=params[2]+params[1] # total effect for mesofauna
params[3]=params[3]+params[1] # total effect for microbes
params<-as.data.frame(params)
params$names<-row.names(params)
params<-params %>% separate(names,c("trophic_group","parameter"),sep=":",fill="right")
params<-params %>% separate(trophic_group,c(NA,"trophic_group"),sep=13)
params$parameter[is.na(params$parameter)]="R0"
params$parameter<-as.factor(params$parameter)
levels(params$parameter)=c("s_R","R0","E_R")
params<-dcast(params,trophic_group~parameter,value.var="params")
params

data<-data.frame(bodymass_fresh=10^seq(-9,4,1))
data$bodymass=fresh_to_dry*dry_to_C*data$bodymass_fresh
data$Microbe<-allometric_parameter(data$bodymass,1,exp(params$R0[3]),params$s_R[3],params$E_R[3],kB,celsius_to_kelvin+15)
data$Mesofauna<-allometric_parameter(data$bodymass,1,exp(params$R0[2]),params$s_R[2],params$E_R[2],kB,celsius_to_kelvin+15)
data$Macrofauna<-allometric_parameter(data$bodymass,1,exp(params$R0[1]),params$s_R[1],params$E_R[1],kB,celsius_to_kelvin+15)
data<-melt(data, id.vars = c("bodymass_fresh","bodymass"),
           variable.name = "trophic_group", 
           value.name = "respiration")

p1<-ggplot(data=data)+
  geom_point(data=data_metabo,aes(bodymass/fresh_to_dry/dry_to_C,metabolic_rate,colour=trophic_group))+
  geom_line(aes(bodymass_fresh,respiration,colour=trophic_group),size=1.5)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE))+
  theme+theme(legend.title=element_blank())
legend<-get_legend(p1)

p1<-ggplot(data=data)+
      geom_point(data=data_metabo,aes(bodymass/fresh_to_dry/dry_to_C,metabolic_rate,colour=trophic_group))+
      geom_line(aes(bodymass_fresh,respiration,colour=trophic_group),size=1.5)+
      scale_colour_viridis(discrete=T,
                           guide=guide_legend(reverse=TRUE))+
      theme+theme(legend.position="none",
                  legend.title=element_blank())+
      scale_x_log10(breaks = 10^seq(-10,5,2),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      y_axis_log10+
      xlab(expression("Fresh body mass (mg)"))+
      ylab(expression("Respiration rate (mgC mgC"^-1*" day"^-1*")"))

# linear regression on the residuals to graphically disentangle the effects of body mass and température
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~trophic_group+temperature_arrhenius:trophic_group)
residuals_bodymass<-model$residuals
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~trophic_group+log(bodymass):trophic_group)
residuals_temperature<-model$residuals
data<-cbind(data_metabo,residuals_bodymass)
data<-cbind(data,residuals_temperature)

p2<-ggplot(data=data)+
      geom_point(data=data_metabo,aes(bodymass/fresh_to_dry/dry_to_C,residuals_bodymass,colour=trophic_group))+
      geom_smooth(aes(bodymass/fresh_to_dry/dry_to_C,residuals_bodymass,colour=trophic_group),method="lm",se=FALSE)+
      scale_colour_viridis(discrete=T,
                           guide=guide_legend(reverse=TRUE))+
      theme+theme(legend.position="none",
                  legend.title=element_blank())+
      x_axis_log10+
      #y_axis_log10+
      xlab(expression("Fresh body mass (mg)"))+
      ylab("Partial residuals")

p3<-ggplot(data=data)+
      geom_point(data=data_metabo,aes(temperature,residuals_temperature,colour=trophic_group))+
      geom_smooth(aes(temperature,residuals_temperature,colour=trophic_group),method="lm",se=FALSE)+
      scale_colour_viridis(discrete=T,
                           guide=guide_legend(reverse=TRUE))+
      theme+theme(legend.position="none",
                  legend.title=element_blank())+
      #x_axis_log10+
      #y_axis_log10+
      xlab("Temperature (°C)")+
      ylab("Partial residuals")

# mass specific metabolic rate calculated for the average body mass of the main soil taxa
data_bodymass<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv",header=T,sep=",")
data_bodymass$bodymass<-data_bodymass$Individual_bodymass_mg_DM*dry_to_C # bondy mass in C
data_bodymass$bodymass_fresh<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry # fresh body mass
data_bodymass$metabolic_rate<-allometric_parameter(data_bodymass$bodymass,1,exp(params$R0[1]),params$s_R[1],params$E_R[1],kB,celsius_to_kelvin+15)
data_bodymass$metabolic_rate[data_bodymass$Fauna_group=="Microbe"]<-allometric_parameter(data_bodymass$bodymass[data_bodymass$Fauna_group=="Microbe"],1,exp(params$R0[3]),params$s_R[3],params$E_R[3],kB,celsius_to_kelvin+15)
data_bodymass$metabolic_rate[data_bodymass$Fauna_group=="Mesofauna"]<-allometric_parameter(data_bodymass$bodymass[data_bodymass$Fauna_group=="Mesofauna"],1,exp(params$R0[2]),params$s_R[2],params$E_R[2],kB,celsius_to_kelvin+15)

p4<-ggplot(data=data_bodymass)+
      geom_point(data=data_metabo,aes(bodymass/fresh_to_dry/dry_to_C,metabolic_rate,colour=trophic_group),alpha=0.2)+
      geom_point(aes(bodymass_fresh,metabolic_rate,fill=Fauna_group),shape = 21,colour="red",size=4)+
      geom_label_repel(aes(bodymass_fresh,metabolic_rate,label=Soil_biota_group,colour=Fauna_group),
                       point.padding=0.5,max.overlaps = Inf,size=5,
                       fill = "white")+ # allow the labels to go beyond the edges of the panel
      scale_colour_viridis(discrete=T,
                           guide=guide_legend(reverse=TRUE))+
      scale_fill_viridis(discrete=T,
                         guide=guide_legend(reverse=TRUE))+
      theme+theme(legend.position="none",
                  legend.title=element_blank())+
      scale_x_log10(breaks = 10^seq(-10,5,2),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      y_axis_log10+
      xlab(expression("Fresh body mass (mg)"))+
      ylab(expression("Respiration rate (mgC mgC"^-1*" day"^-1*")"))

graph<-ggdraw(xlim = c(0, 2.3), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 1, 1)+
  draw_plot(p4, 1, 1, 1, 1)+
  draw_plot(p2, 0, 0, 1, 1)+
  draw_plot(p3, 1, 0, 1, 1)+
  draw_plot(legend, 2.05, 0.75, 0.15, 0.5)+
  draw_plot_label(c("A","B","C","D"), c(0,1,0,1), c(2,2,1,1), size = 30)
ggsave(paste(path_figure,"supp_metabolism.pdf",sep=""),graph, width = 18, height = 12, device = cairo_pdf)

### MORTALITY ### ----
# data from McCoy and Gillooly (2008) (10.1111/j.1461-0248.2008.01190.x)
data_mort<-read.csv(file=paste(path_data,"data_mortality_McCoy_2008.csv",sep=""),
                      header=TRUE)
data_mort<-na.omit(data_mort)
data_mort<-data_mort[data_mort$mortality>0,]
T0=293.15 # reference temperature 20°C
# dry_mass : g
# temperature : celsius
# mortality : year-1
names(data_mort)<-c("group","species","bodymass","temperature","mortality","reference")
data_metabo$temperature<--1/(celsius_to_kelvin+data_metabo$temperature)*1/kB
data_mort$group<-as.factor(data_mort$group) #levels(data_mort$group)
data_mort$bodymass<-data_mort$bodymass*4*fresh_to_dry*dry_to_C*1e3 # conversion from g (dry) to mgC
data_mort$mortality<-data_mort$mortality/year_to_day
# model with all taxons
model<-lm(data=data_mort,
          formula=log(mortality)~group+log(bodymass):group+temperature:group)
par(mfrow=c(2,2))
plot(model)
params<-summary(model)
params
Anova(model,type=2,test.statisticf="F")

params<-(params$coefficients)[,1]
names(params)[1]="groupBird"
for (i in 2:5){
  params[i]=params[i]+params[1] # total effect for each group for \mu0
}
params<-as.data.frame(params)
params$names<-row.names(params)
params<-params %>% separate(names,c("group","parameter"),sep=":",fill="right")
params<-params %>% separate(group,c(NA,"group"),sep=5)
params$parameter[is.na(params$parameter)]="mu0"
params$parameter<-as.factor(params$parameter)
levels(params$parameter)=c("s_mu","mu0","E_mu")
params<-dcast(params,group~parameter,value.var="params")

# model for invertebrates only
data_mort<-data_mort[data_mort$group=="Invertebrate",]
model<-lm(data=data_mort,
          formula=log(mortality)~log(bodymass)+temperature)
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")

# linear mortality and self-regulation
# data<-expand.grid(predator_bodymass_fresh=seq(-7,4,1),
#                   prey_biomass=10^seq(-2,6,0.1))
# data$predator_bodymass=fresh_to_dry*dry_to_C*10^data$predator_bodymass_fresh
# data$FR<-functional_response(data$prey_biomass,data$predator_bodymass,a0,s_a,E_a,h0,s_h,E_h,kB,293.15)
# #data$a<-allometric_parameter(data$predator_bodymass,a0,s_a,E_a,kB,293.15)
# data$h<-allometric_parameter(data$predator_bodymass,h0,s_h,E_h,kB,293.15)
# 
# p1<-ggplot(data=data)+
#   geom_line(aes(prey_biomass,FR,colour=as.factor(predator_bodymass_fresh)),size=1.5)+
#   scale_colour_viridis(discrete=T,
#                        guide=guide_legend(reverse=TRUE),
#                        name="log10 of predator\nfresh body mass (mg)")+
#   theme+
#   x_axis_log10+
#   y_axis_log10+
#   xlab(expression("Prey biomass density (mgC "*m^-2*")"))+
#   ylab(expression("Biomass eaten (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))
# ggsave(paste(path_figure,"supp_mortality.pdf",sep=""),p1, width = 9, height = 5.5, device = cairo_pdf)

### FUNCTIONAL RESPONSE OF CARNIVORES ### ----
# data from Li et al. 2017 (Oikos)
T0=293.15 # reference temperature 20°C
# time duration and handling time : s
# attack rate : m²/s
# body mass : mg (fresh)
# temperature : °C
# surface : m²
data_FR<-read.csv(file=paste(path_data,"FR_li_2017_oikos.csv",sep=""),header=TRUE)
#data_FR<-data_FR[which(data_FR$ecosystem.type=="terrestrial"),]
#data_FR<-data_FR[which(data_FR$dimensionality=="2D"),]
#data_FR<-data_FR[which(data_FR$predator.ana.group=="invertebrate"),]
data_FR<-data_FR[,which(names(data_FR)%in%c("dimensionality","experimental.duration.seconds","temperature.degree.celcius","predator.mass.mg","prey.mass.mg","Starvation.Y.N","attack.rate","handling.time","original.publication.short","arena.size"))]
names(data_FR)=c("publication","dimensionality","arena_size","duration","temperature","predator_mass","prey_mass","attack_rate","handling_time","starvation")
data_FR$temperature<-celsius_to_kelvin+data_FR$temperature

# original values before conversion (same results as Li)
data_FR$temperature_arrhenius<-(data_FR$temperature-T0)/(kB*data_FR$temperature*T0)
# attack rate
model=lme(data=data_FR,
         log(attack_rate)~dimensionality+log(predator_mass)+temperature_arrhenius+log(duration),
         random = ~ 1+ log(predator_mass)| publication,
         method = "ML",
         control=lmeControl(opt="optim",
                            msMaxIter = 200,
                            msMaxEval = 500,
                            msVerbose = F))
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")
# handling time
model=lme(data=data_FR,
          log(handling_time)~starvation+log(predator_mass)+temperature_arrhenius,
          random = ~ 1| publication/arena_size,
          method = "ML",
          control=lmeControl(opt="optim",
                             msMaxIter = 200,
                             msMaxEval = 500,
                             msVerbose = F))
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")

# conversion
data_FR$handling_time<-second_to_day*data_FR$handling_time*data_FR$predator_mass # conversion from individual based (s Ind_pred mg_prey-1) handling time to biomass based
data_FR$predator_mass<-fresh_to_dry*dry_to_C*data_FR$predator_mass # already in mg
data_FR$duration<-second_to_day*data_FR$duration
data_FR$attack_rate<-data_FR$attack_rate/(second_to_day*data_FR$predator_mass) # conversion from individual based attack rate to biomass based
data_FR$temperature<--1/(data_FR$temperature*kB)
# attack rate
model=lme(data=data_FR,
          log(attack_rate)~dimensionality+log(predator_mass)+temperature, # duration is dropped +log(duration)
          random = ~ 1+ log(predator_mass)| publication,
          method = "ML",
          control=lmeControl(opt="optim",
                             msMaxIter = 200,
                             msMaxEval = 500,
                             msVerbose = F))
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")
params<-summary(model)
params<-(params$coefficients)[[1]]
a0=exp(params[1])
s_a=params[3]
E_a=params[4]

# handling time
model=lme(data=data_FR,
          log(handling_time)~log(predator_mass)+temperature,
          random = ~ 1| publication/arena_size,
          method = "ML",
          control=lmeControl(opt="optim",
                             msMaxIter = 200,
                             msMaxEval = 500,
                             msVerbose = F))
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")
params<-summary(model)
params<-(params$coefficients)[[1]]
h0=exp(params[1])
s_h=params[2]
E_h=params[3]

# functional response
data<-expand.grid(predator_bodymass_fresh=seq(-7,4,1),
                  prey_biomass=10^seq(-2,6,0.1))
data$predator_bodymass=fresh_to_dry*dry_to_C*10^data$predator_bodymass_fresh
data$FR<-functional_response(data$prey_biomass,data$predator_bodymass,1,a0,s_a,E_a,1,h0,s_h,E_h,kB,293.15)
#data$a<-allometric_parameter(data$predator_bodymass,a0,s_a,E_a,kB,293.15)
data$h<-allometric_parameter(data$predator_bodymass,1,h0,s_h,E_h,kB,293.15)

p1<-ggplot(data=data)+
  geom_line(aes(prey_biomass,FR,colour=as.factor(predator_bodymass_fresh)),size=1.5)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE),
                       name="log10 of predator\nfresh body mass (mg)")+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Prey biomass density (mgC "*m^-2*")"))+
  ylab(expression("Biomass eaten (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))
ggsave(paste(path_figure,"supp_FR_carnivores.pdf",sep=""),p1, width = 9, height = 5.5, device = cairo_pdf)

ggplot(data=data)+
  geom_line(aes(10^predator_bodymass_fresh,1/h),size=1.5)+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab("Predator fresh body mass (mg)")+
  ylab(expression("Maximum ingestion rate ("*italic(h)*") (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))

data<-expand.grid(predator_bodymass_fresh=c(-7,4),
                  prey_biomass=10^seq(-2,6,0.1),
                  a_a=10^seq(-6,0),
                  a_h=10^seq(0,6))
data$predator_bodymass=fresh_to_dry*dry_to_C*10^data$predator_bodymass_fresh
data$FR<-functional_response(data$prey_biomass,data$predator_bodymass,data$a_a,a0,s_a,E_a,data$a_h,h0,s_h,E_h,kB,293.15)
data$predator_bodymass_fresh<-as.factor(data$predator_bodymass_fresh)
levels(data$predator_bodymass_fresh)<-c(expression(10^-7*" mg (fw)"),expression(10^4*" mg (fw)"))

p1<-ggplot(data=data[data$a_h==1,])+
  geom_line(aes(prey_biomass,FR,colour=as.factor(a_a)),size=1.5)+
  facet_wrap(~predator_bodymass_fresh,labeller=label_parsed)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE),
                       name=expression(italic(a[a])))+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Prey biomass density (mgC "*m^-2*")"))+
  ylab(expression("Biomass eaten (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))

p2<-ggplot(data=data[data$a_a==1,])+
  geom_line(aes(prey_biomass,FR,colour=as.factor(a_h)),size=1.5)+
  facet_wrap(~predator_bodymass_fresh,labeller=label_parsed)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE),
                       name=expression(italic(a[h])))+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Prey biomass density (mgC "*m^-2*")"))+
  ylab(expression("Biomass eaten (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figure,"supp_FR_carnivores_tuning.pdf",sep=""),graph, width = 18, height = 6, device = cairo_pdf)

### FUNCTIONAL RESPONSE OF DETRITIVORES ### ----
# data from Ott et al. 2012 (PRSB) # ----
data_ott<-read.csv(file=paste(path_data,"data_Ott_2012.csv",sep=""),sep=";",header=TRUE)
data_ott<-data_ott[data_ott$predno=="3",] # removes controls, spontaneous decay is negligible
data_ott<-data_ott[data_ott$litter!="buche",] # removes data with beech tree
data_ott<-data_ott[,which(names(data_ott)%in%c("day","temp","size","litter","replica","id","startmass","endmass","meanisopodweight","meanindividualweight"))]
area=pi*(0.063/2)^2 # surface of a jar (6.3cm diameter)

# REPRODUCTION OF THE RESULTS OF OTT ET ALL (2012) # ----
# fit of the attack rate and the handling time #
# see Bolker (2008) (ISBN: 978-0-691-12522-0) p.355 and Ott et al. (2012) for the detailed method
data_ott$id<-as.factor(data_ott$id)
data_ott$startmass<-data_ott$startmass/area*1e-3
data_ott$endmass<-data_ott$endmass/area*1e-3
data_ott$meanisopodweight<-data_ott$meanisopodweight/area*1e-3
data_ott$meanindividualweight<-data_ott$meanindividualweight*1e-3
data_ott$temp<-data_ott$temp+celsius_to_kelvin
id=levels(data_ott$id)
id=levels(data_ott$id)
n=length(id)
TS<-vector("list",n)
for(i in 1:n){
  TS[[i]]=data_ott[which(data_ott$id%in%id[i]),]
  TS[[i]]$meanisopodweight=TS[[i]]$meanisopodweight[1]
  TS[[i]]$meanindividualweight=TS[[i]]$meanindividualweight[1]
  if(is.na(TS[[i]]$meanisopodweight[1])){
    TS[[i]]$meanisopodweight=TS[[i]]$meanindividualweight[1]*3/area # 3 individuals are in each jar
  }
  data_ott$meanindividualweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanindividualweight[1]
}

# first fitting to initialise a complete fitting with all variables later
# Rogers random-predator equation
rogers_allomtric<-function(litter,startmass, meanindividualweight, day, temp,
                           a0_esche, s_a_esche, E_a_esche, h0_esche, s_h_esche, E_h_esche,
                           a0_hainbu, s_a_hainbu, E_a_hainbu, h0_hainbu, s_h_hainbu, E_h_hainbu){
  a=(litter=="esche")*allometric_rate(a0_esche, s_a_esche, E_a_esche, temp, meanindividualweight)+
    (litter=="hainbu")*allometric_rate(a0_hainbu, s_a_hainbu, E_a_hainbu, temp, meanindividualweight)
  h=(litter=="esche")*allometric_rate(h0_esche, s_h_esche, E_h_esche, temp, meanindividualweight)+
    (litter=="hainbu")*allometric_rate(h0_hainbu, s_h_hainbu, E_h_hainbu, temp, meanindividualweight)
  return(lambertW(a * h * startmass * exp(-a * (3 * day - h * startmass)))/(a * h)) # returns the remaining biomass of prey
}
allometric_rate<-function(x0, s_x, E_x, temp, meanindividualweight){
  kB=8.62e-5
  return(x0*meanindividualweight^s_x*exp(-E_x/(kB*temp)))
}
rogers_ind<-function(startmass, day, a, h){
  return(lambertW(a * h * startmass * exp(-a * (3 * day - h * startmass)))/(a * h)) # returns the remaining biomass of prey
}

init=list(a0_esche=exp(13),
          s_a_esche=-0.5,
          E_a_esche=0.43,
          h0_esche=exp(-14.32),
          s_h_esche=-0.75,
          E_h_esche=-0.3,
          a0_hainbu=exp(13),
          s_a_hainbu=-0.5,
          E_a_hainbu=0.43,
          h0_hainbu=exp(-14.32),
          s_h_hainbu=-0.75,
          E_h_hainbu=-0.3)

data<-data_ott
data$fit<-rogers_allomtric(data$litter, data$startmass, data$meanindividualweight, data$day, data$temp,
                           init$a0_esche, init$s_a_esche, init$E_a_esche, init$h0_esche, init$s_h_esche, init$E_h_esche,
                           init$a0_hainbu, init$s_a_hainbu, init$E_a_hainbu, init$h0_hainbu, init$s_h_hainbu, init$E_h_hainbu)

a=attack_rate(init$a0, init$s_a, init$E_a, data$temp[1], data$meanindividualweight[1])
h=handling_time(init$h0, init$s_h, init$E_h, data$temp[1], data$meanindividualweight[1])
a * h * data$startmass[1] * exp(-a * (data$meanisopodweight[1] * data$day[1] - h * data$startmass[1]))

fit<-nls(endmass~rogers_allomtric(litter,startmass, meanindividualweight, day, temp,
                                  a0_esche, s_a_esche, E_a_esche, h0_esche, s_h_esche, E_h_esche,
                                  a0_hainbu, s_a_hainbu, E_a_hainbu, h0_hainbu, s_h_hainbu, E_h_hainbu),
         data_ott[data_ott$endmass>0,],
         start=init,
         algorithm="plinear",
         control = list(maxiter=50,warnOnly=T))
summary(fit)

fit<-nls(endmass~rogers_ind(startmass, day, a, h),
         data_ott,
         start=list(a=0.01,h=15),
         algorithm="port",
         control = list(maxiter=50,warnOnly=T),
         lower=c(1e-8,1e-4), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e2,1e4))
fit<-nls(endmass~rogers_ind(startmass, day, a, h),
         data_ott,
         start=list(a=coef(fit)[1],h=coef(fit)[2]),
         algorithm="port",
         control = list(maxiter=50,warnOnly=T),
         lower=c(1e-8,1e-4), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e2,1e4))
summary(fit)

# TEST # ----
# see Bolker (2008) (ISBN: 978-0-691-12522-0) p.355 and Ott et al. (2012) for the detailed method
data_ott$id<-as.factor(data_ott$id)
data_ott$startmass<-data_ott$startmass/area*1e-3
data_ott$endmass<-data_ott$endmass/area*1e-3
data_ott$meanisopodweight<-data_ott$meanisopodweight/area*1e-3
data_ott$meanindividualweight<-data_ott$meanindividualweight*1e-3
id=levels(data_ott$id)
n=length(id)
TS<-vector("list",n)
for(i in 1:n){
  TS[[i]]=data_ott[which(data_ott$id%in%id[i]),]
  TS[[i]]$meanisopodweight=TS[[i]]$meanisopodweight[1]
  TS[[i]]$meanindividualweight=TS[[i]]$meanindividualweight[1]
  if(is.na(TS[[i]]$meanisopodweight[1])){
    TS[[i]]$meanisopodweight=TS[[i]]$meanindividualweight[1]*3/area # 3 individuals are in each jar
  }
  data_ott$meanindividualweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanindividualweight[1]
  data_ott$meanisopodweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanisopodweight[1]
}
TS[[1]]

# estimation of initial values
holling<-data_ott[data_ott$day==23,]
ggplot(data=holling)+
  geom_point(aes(startmass,(startmass-endmass),colour=temp))+
  facet_grid(as.factor(temp)~litter,scale="free")

litter<-levels(as.factor(holling$litter))
init<-expand.grid(litter=litter,
                  temp=c(10,15,20))
init$a0=0 # estimated initial values of the attack rate for the fitting
init$h0=0 # estimated initial values of the handling time for the fitting
for(i in 1:2){
  for(j in c(10,15,20)){
    subset<-holling[which(holling$litter==litter[i]&holling$temp==j),]
    index_max<-which(subset$startmass==max(subset$startmass))
    linear<-lm(data=subset,(startmass-endmass)~startmass)
    coef<-summary(linear)$coefficients
    init$a0[which(init$litter==litter[i]&init$temp==j)]=coef[2,1]
    init$h0[which(init$litter==litter[i]&init$temp==j)]=subset$meanisopodweight[index_max]/(subset$startmass[index_max]-subset$endmass[index_max])
  }
}

# fitting
# Rogers random-predator equation
rogers<-function(startmass, meanisopodweight, a, h, day){
  lambertW(a * h * startmass * exp(-a * (meanisopodweight * day - h * startmass)))/(a * h) # returns the remaining biomass of prey
}

data<-data_ott[,which(names(TS[[1]])%in%c("day","endmass","startmass","meanisopodweight","temp"))]
#data<-data[data$startmass-data$endmass>0,]
a0=mean(init$a0)
h0=mean(init$h0)
fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
         data,
         start=list(a=a0,h=h0),
         algorithm="port",
         control = list(maxiter = 50, warnOnly=T),
         lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e4,1e4))
coef<-coef(fit)
fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
         data,
         start=list(a=coef[1],h=coef[2]),
         algorithm="port",
         control = list(warnOnly=T),
         lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e4,1e4))

# FINAL FITTING # ----
# fit of the attack rate and the handling time #
# see Bolker (2008) (ISBN: 978-0-691-12522-0) p.355 and Ott et al. (2012) for the detailed method
data_ott$id<-as.factor(data_ott$id)
data_ott$startmass<-data_ott$startmass/area*fresh_to_dry_litter*dry_to_C*1e-3
data_ott$endmass<-data_ott$endmass/area*fresh_to_dry_litter*dry_to_C*1e-3
data_ott$meanisopodweight<-data_ott$meanisopodweight/area*fresh_to_dry*dry_to_C*1e-3
data_ott$meanindividualweight<-data_ott$meanindividualweight*fresh_to_dry*dry_to_C*1e-3
id=levels(data_ott$id)
n=length(id)
TS<-vector("list",n)
for(i in 1:n){
  TS[[i]]=data_ott[which(data_ott$id%in%id[i]),]
  TS[[i]]$meanisopodweight=TS[[i]]$meanisopodweight[1]
  TS[[i]]$meanindividualweight=TS[[i]]$meanindividualweight[1]
  if(is.na(TS[[i]]$meanisopodweight[1])){
    TS[[i]]$meanisopodweight=TS[[i]]$meanindividualweight[1]*3/area # 3 individuals are in each jar
  }
  data_ott$meanindividualweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanindividualweight[1]
  data_ott$meanisopodweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanisopodweight[1]
}

# estimation of initial values
holling<-data_ott[data_ott$day==27,]
ggplot(data=holling)+
  geom_point(aes(startmass,(startmass-endmass),colour=replica))+
  facet_grid(as.factor(temp)~litter,scale="free")

litter<-levels(as.factor(holling$litter))
size<-levels(as.factor(holling$size))
init<-expand.grid(size=size,
                  temp=c(10,15,20))
init$a0=0 # estimated initial values of the attack rate for the fitting
init$h0=0 # estimated initial values of the handling time for the fitting
for(i in 1:length(size)){
  for(j in c(10,15,20)){
    subset<-holling[which(holling$size==size[i]&holling$temp==j),]
    index_max<-which(subset$startmass==max(subset$startmass))
    linear<-lm(data=subset,(startmass-endmass)~startmass) # a is roughly equivalent to the slope and h to P/max_consumption
    coef<-summary(linear)$coefficients
    init$a0[which(init$size==size[i]&init$temp==j)]=coef[2,1]
    init$h0[which(init$size==size[i]&init$temp==j)]=subset$meanisopodweight[index_max]/(subset$startmass[index_max]-subset$endmass[index_max])
  }
}
a0=mean(init$a0)
h0=mean(init$h0)

# first fitting to initialise a complete fitting with all variables later
# Rogers random-predator equation
rogers<-function(startmass, meanisopodweight, a, h, day){
  lambertW(a * h * startmass * exp(-a * (meanisopodweight * day - h * startmass)))/(a * h) # returns the remaining biomass of prey
}

init$a=0
init$h=0
for(i in 1:dim(init)[1]){
  data<-data_ott[which(data_ott$size==init$size[i]&data_ott$temp==init$temp[i]&data_ott$litter=="esche"),]
  fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
           data,
           start=list(a=init$a0[i],h=init$h0[i]),
           algorithm="port",
           control = list(warnOnly=T,scaleOffset = 1),
           lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
           upper=c(1e4,1e4))
  coef<-coef(fit)
  fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
           data,
           start=list(a=coef[1],h=coef[2]),
           algorithm="port",
           control = list(warnOnly=T),
           lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
           upper=c(1e4,1e4))
  init$a[i]=coef[1]
  init$h[i]=coef[2]
}
results<-as.data.frame(matrix(0,n,2))
names(results)<-c("a","h")
for(i in 1:n){
  fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
           TS[[i]],
           start=list(a=0.01,h=15),
           algorithm="port",
           control = list(warnOnly=T,scaleOffset = 1),
           lower=c(1e-15,1e-15), # having a lower bond is vital, otherwise negative parameters can be estimated
           upper=c(1e15,1e15))
  coef<-coef(fit)
  fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
           TS[[i]],
           start=list(a=coef[1],h=coef[2]),
           algorithm="port",
           control = list(warnOnly=T),
           lower=c(1e-10,1e-10), # having a lower bond is vital, otherwise negative parameters can be estimated
           upper=c(1e8,1e8))
  coef<-coef(fit)
  results$a[i]=coef[1]
  results$h[i]=coef[2]
}
ggplot(data=TS[[2]])+
  geom_point(aes(day,endmass))
#data<-data[data$startmass-data$endmass>0,]
a0=mean(init$a0)
h0=mean(init$h0)
fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
         data,
         start=list(a=a0,h=h0),
         algorithm="port",
         control = list(maxiter = 50, warnOnly=T),
         lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e4,1e4))
coef<-coef(fit)
fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
         data,
         start=list(a=coef[1],h=coef[2]),
         algorithm="port",
         control = list(warnOnly=T),
         lower=c(1e-8,1e-8), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e4,1e4))
# GENERAL FITTING # ----
# fit of the attack rate and the handling time #
# see Bolker (2008) (ISBN: 978-0-691-12522-0) p.355 and Ott et al. (2012) for the detailed method
data_ott$id<-as.factor(data_ott$id)
data_ott$startmass<-data_ott$startmass/area*fresh_to_dry_litter*dry_to_C*1e-3
data_ott$endmass<-data_ott$endmass/area*fresh_to_dry_litter*dry_to_C*1e-3
data_ott$meanisopodweight<-data_ott$meanisopodweight/area*fresh_to_dry*dry_to_C*1e-3
data_ott$meanindividualweight<-data_ott$meanindividualweight*fresh_to_dry*dry_to_C*1e-3
data_ott$temp<-data_ott$temp+celsius_to_kelvin
id=levels(data_ott$id)
id=levels(data_ott$id)
n=length(id)
TS<-vector("list",n)
for(i in 1:n){
  TS[[i]]=data_ott[which(data_ott$id%in%id[i]),]
  TS[[i]]$meanisopodweight=TS[[i]]$meanisopodweight[1]
  TS[[i]]$meanindividualweight=TS[[i]]$meanindividualweight[1]
  if(is.na(TS[[i]]$meanisopodweight[1])){
    TS[[i]]$meanisopodweight=TS[[i]]$meanindividualweight[1]*3/area # 3 individuals are in each jar
  }
  data_ott$meanindividualweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanindividualweight[1]
  data_ott$meanisopodweight[which(data_ott$id%in%id[i])]=TS[[i]]$meanisopodweight[1]
}

# first fitting to initialise a complete fitting with all variables later
# Rogers random-predator equation
rogers_allomtric<-function(startmass, meanisopodweight, meanindividualweight, day, temp, a0, s_a, E_a, h0, s_h, E_h){
  a=attack_rate(a0, s_a, E_a, temp, meanindividualweight)
  h=handling_time(h0, s_h, E_h, temp, meanindividualweight)
  return(lambertW(a * h * startmass * exp(-a * (meanisopodweight * day - h * startmass)))/(a * h)) # returns the remaining biomass of prey
}
attack_rate<-function(a0, s_a, E_a, temp, meanindividualweight){
  kB=8.62e-5
  return(a0*meanindividualweight^s_a*exp(-E_a/(kB*temp)))
}
handling_time<-function(h0, s_h, E_h, temp, meanindividualweight){
  kB=8.62e-5
  return(h0*meanindividualweight^s_h*exp(-E_h/(kB*temp)))
}

init=list(a0=exp(13),
          s_a=-0.5,
          E_a=0.43,
          h0=exp(-14.32),
          s_h=-0.75,
          E_h=-0.3)

attack_rate(init$a0, init$s_a, init$E_a, 288.3, 0.003)
handling_time(init$h0, init$s_h, init$E_h, 288.3, 0.003)
data<-data_ott
data$a<-attack_rate(init$a0, init$s_a, init$E_a, data$temp, data$meanindividualweight)
data$h<-handling_time(init$h0, init$s_h, init$E_h, data$temp, data$meanindividualweight)
data$fit<-rogers_allomtric(data$startmass, data$meanisopodweight, data$meanindividualweight, data$day, data$temp, init$a0, init$s_a, init$E_a, init$h0, init$s_h, init$E_h)

a=attack_rate(init$a0, init$s_a, init$E_a, data$temp[1], data$meanindividualweight[1])
h=handling_time(init$h0, init$s_h, init$E_h, data$temp[1], data$meanindividualweight[1])
a * h * data$startmass[1] * exp(-a * (data$meanisopodweight[1] * data$day[1] - h * data$startmass[1]))

fit<-nls(endmass~rogers_allomtric(startmass, meanisopodweight, meanindividualweight, day, temp, a0, s_a, E_a, h0, s_h, E_h),
         data_ott[data_ott$endmass>0,],
         start=init,
         algorithm="plinear",
         control = list(warnOnly=T))
coef<-coef(fit)

fit<-nls(endmass~rogers(startmass, meanisopodweight, a, h, day),
         data,
         start=list(a=0.01,h=15),
         algorithm="port",
         control = list(warnOnly=T),
         lower=c(1e-5,1e-5), # having a lower bond is vital, otherwise negative parameters can be estimated
         upper=c(1e3,1e3))
summary(fit)

init$a[i]=coef[1]
init$h[i]=coef[2]

# DIRECT CONVERSION # ----
T0=288.15 # 15°C
# attack rate original values
s_a_1=0.201 # scaling for individual attack rate
s_a_2=0.745 # scaling for individual attack rate
s_a=(s_a_1+s_a_2)/2
a0_1=2.74e-4 # m2 d-1 g-s_a
a0_2=1.18e-3 # m2 d-1 g-s_a
a0=(a0_1+a0_2)/2
E_a_1=2.3e-2 # eV
E_a_2=0.274 # eV
E_a=(E_a_1+E_a_2)/2
# conversion
s_a=s_a-1
a0=a0*exp(E_a/(kB*T0))*1e-3

# handling time
s_h=-0.165 # appropriate scaling
h0_1=2.07e2 # day g -s_h
h0_2=5e2 # day g -s_h
h0=(h0_1+h0_2)/2
E_h_1=-1.36 # eV
E_h_2=-0.86 # eV
E_h=(E_h_1+E_h_2)/2
# conversion
s_h=s_h+1
h0=h0*exp(E_h/(kB*T0))

# functional response
data<-expand.grid(predator_bodymass_fresh=seq(-7,4,1),
                  prey_biomass=10^seq(-2,6,0.1))
data$predator_bodymass=fresh_to_dry*dry_to_C*10^data$predator_bodymass_fresh
data$FR<-functional_response(data$prey_biomass,data$predator_bodymass,1,a0,s_a,E_a,1,h0,s_h,E_h,kB,293.15)
data$h<-allometric_parameter(data$predator_bodymass,1,h0,s_h,E_h,kB,293.15)

p1<-ggplot(data=data)+
  geom_line(aes(prey_biomass,FR,colour=as.factor(predator_bodymass_fresh)),size=1.5)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE),
                       name="log10 of detritivore\nfresh body mass (mg)")+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Detritus biomass density (mgC "*m^-2*")"))+
  ylab(expression("Biomass eaten (mgC"["det"]*" mgC"["pred"]^-1*" day"^-1*")"))
ggsave(paste(path_figure,"supp_FR_detritivores.pdf",sep=""),p1, width = 9, height = 5.5, device = cairo_pdf)

ggplot(data=data)+
  geom_line(aes(10^predator_bodymass_fresh,1/h),size=1.5)+
  theme+
  x_axis_log10+
  y_axis_log10+
  xlab("Predator fresh body mass (mg)")+
  ylab(expression("Maximum ingestion rate ("*italic(h)*") (mgC"["prey"]*" mgC"["pred"]^-1*" day"^-1*")"))

### FAECES ALLOMETRY ### ----
data_faeces<-read.csv(file=paste(path_data,"data_faeces_Ganault_2022.csv",sep=""),header=TRUE)
data_faeces$area<-data_faeces$faecal_specific_area*data_faeces$faeces_dry_mass
data_faeces$radius<-sqrt(data_faeces$area)/pi
data_faeces$body_shape="compact"
data_faeces$body_shape[which(data_faeces$species%in%c("Cylindroiulus caerulocinctus","Ommatoiuls sabulosus"))]="elongated"
# conversion of units
data_faeces$average_body_mass<-data_faeces$average_body_mass*1e3*fresh_to_dry*dry_to_C # conversion from fresh g to mgC
data_faeces$faeces_wet_mass<-data_faeces$faeces_wet_mass*1e3 # conversion from g to mg
data_faeces$faeces_dry_mass<-data_faeces$faeces_dry_mass*1e3 # conversion from g to mg
data_faeces$radius<-data_faeces$radius*1e3 # conversion from mm to µm

# faecal pellet radius
p1<-ggplot(data=data_faeces)+
  geom_point(aes(average_body_mass,radius,colour=species))+
  geom_smooth(aes(average_body_mass,radius,linetype=body_shape),method="lm")+
  scale_linetype_discrete(name="body shape")+
  theme+
  xlab("Body mass (mgC)")+
  ylab("Radius (µm)")
ggsave(paste(path_figure,"supp_faeces_regression.pdf",sep=""),p1, width = 9, height = 4, device = cairo_pdf)

model<-lm(data=data_faeces,
          formula=radius~average_body_mass*body_shape)

summary(model)
Anova(model)

params<-summary(model)
params<-(params$coefficients)[,1]
names(params)<-c("r_int_compact","r_slope_compact","r_int_elongated","r_slope_elongated")
params[3]=params[1]+params[3]
params[4]=params[2]+params[4]
params

# faecal pellet "body mass"
ggplot(data=data_faeces)+
  geom_point(aes(average_body_mass,faeces_wet_mass,colour=species))+
  geom_smooth(aes(average_body_mass,faeces_wet_mass,linetype=body_shape),method="lm")+
  theme

model<-lm(data=data_faeces,
          formula=faeces_wet_mass~average_body_mass*body_shape)

summary(model)
Anova(model)

databis<-data_faeces[which(data_faeces$body_shape=="compact"),]
model<-lm(data=databis,
          formula=faeces_wet_mass~average_body_mass)

summary(model)
Anova(model)

# we just do a general ratio to establish the allometric relationship between invertebrate body mass and their faeces mass

ratio=mean(data_faeces$faeces_wet_mass)/mean(data_faeces$average_body_mass)
ratio

### MICROBES ### ----
molSub_to_molC=6 # there are 6C in one molecule of substrate, the 10 other C atoms are from the fluorescent dye
temperature=293.15

# half-saturation K of FOM and SOM converted to mgC m-2
s_K=0.034
K0=24.3 # µmol L-1 : original value from German et al (2012) after correction of unit mistakes
K0_OM=K0*1e-6*molSub_to_molC*molC_to_mgC/(1e-3*water_to_soil*volume_to_surface*exp(celsius_to_kelvin*s_K))

# maximum C consumption rate \phi of FOM and SOM converted to day-1 from Perveen et al (2014)
s_phi=0.063
# FOM
phi=0.0093
phi0_FOM=phi/exp(s_phi*temperature)
# SOM
phi=0.0318
phi0_SOM=phi/exp(s_phi*temperature)

# half saturation K of N converted to mgN m-2
K=1e-4 # mol m-3 : original value from Grover (2003)
K0_N=K*molN_to_mgN/(water_to_soil*volume_to_surface)/exp(s_K*temperature)

# maximum N consumption rate \phi of N converted to day-1
phi=7.7e-15 # mol cell-1 day-1
phi0_N=phi*molN_to_mgN/(cell_to_fresh*fresh_to_dry*dry_to_C)/exp(s_phi*temperature)

# half saturation K of DOC converted to mgC m-2
K=1e-3 # mol m-3 : original value from Grover (2003)
K0_DOC=K*molC_to_mgC/(water_to_soil*volume_to_surface*exp(s_K*temperature))

# maximum C consumption rate \phi of DOC converted to day-1
phi=52e-15 # mol cell-1 day-1
phi0_DOC=phi*molC_to_mgC/(cell_to_fresh*fresh_to_dry*dry_to_C)/exp(s_phi*temperature)

data<-data.frame(a=10^seq(-8,3,0.1))
data$K_OM=data$a*K0_OM*exp(s_K*temperature)
data$K_DOC=data$a*K0_DOC*exp(s_K*temperature)
data$K_N=data$a*K0_N*exp(s_K*temperature)
data$phi_FOM=data$a*phi0_FOM*exp(s_phi*temperature)
data$phi_SOM=data$a*phi0_SOM*exp(s_phi*temperature)
data$phi_DOC=data$a*phi0_DOC*exp(s_phi*temperature)
data$phi_N=data$a*phi0_N*exp(s_phi*temperature)

data<-melt(data, id.vars = c("a"),
           variable.name = "parameter", 
           value.name = "value")
levels(data$parameter)<-c(expression(italic(K["FOM,SOM"])),
                          expression(italic(K["DOC"])),
                          expression(italic(K["N"])),
                          expression(italic("\u03C6"["FOM"]^"max")),
                          expression(italic("\u03C6"["SOM"]^"max")),
                          expression(italic("\u03C6"["DOC"]^"max")),
                          expression(italic("\u03C6"["N"]^"max")))

p1<-ggplot(data=data)+
      geom_line(aes(a,value,colour=parameter,linetype=parameter),size=1.5)+
      geom_vline(xintercept=1,linetype="dashed")+
      geom_hline(yintercept=1,linetype="dashed")+
      scale_colour_manual(values=c("orange3","black","chartreuse4","orange1","orange4","black","chartreuse4"),
                          labels = parse_format())+
      scale_linetype_manual(values=c(rep("22",3),rep("solid",4)),
                            labels = parse_format())+
      theme+theme(legend.title=element_blank())+
      scale_x_log10(breaks = c(10^(-8:3)),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      y_axis_log10_short+
      xlab(expression(italic(a[x])))+
      ylab("Uptake parameter value")
ggsave(paste(path_figure,"supp_uptake_tuning.pdf",sep=""),p1, width = 8, height = 5, device = cairo_pdf)

### NONALLOMETRIC PARAMETERS ### ----
# FOM
I_FOM=239 # litter input gC m-2 year-1 Bowden et al 2014
I_FOM=I_FOM*1e3/year_to_day
I_FOM=c(111,667) # coniferous litter input gC m-2 year-1 Wunderlich et al 2012
I_FOM=I_FOM*1e3/year_to_day
I_FOM=c(206,447) # deciduous litter input gC m-2 year-1 Wunderlich et al 2012
I_FOM=I_FOM*1e3/year_to_day
I_FOM<-as.data.frame(matrix(0,nrow=8,ncol=3))
names(I_FOM)<-c("biome","input","SE")
I_FOM$biome<-c("Tundra and cold steppe","Boreal forest","Temperate forest","Temperate grassland",
               "Mediterranean vegetation","Desert","Tropical grassland","Tropical forest")
I_FOM$input<-c(1702,2032,3221,2997,2974,638,3893,5413) # Mean annual litterfall (kg ha−1 year−1) (dry weight I assume)
I_FOM$SE<-c(706,1094,1394,1053,1480,448,1894,1394) # Standard error (kg ha−1 year−1)
I_FOM$input<-I_FOM$input*1e6/ha_to_m2/year_to_day*dry_to_C # conversion to mgC m-2 day-1
I_FOM$SE<-I_FOM$SE*1e6/1e4/year_to_day*dry_to_C # conversion to mgC m-2 day-1
min(I_FOM$input-I_FOM$SE)
max(I_FOM$input+I_FOM$SE)
# DOC
I_DOC=c(0.1,5) # gC kg_soil-1 month-1
I_DOC=I_DOC*1e3/(1e6*soil_mass_to_volume*volume_to_surface)/month_to_day
# Microbe dormancy
Q=1 # maximum transition rate hour-1
Q=Q/hour_to_day
# Net Primary Production (NPP)
data_NPP<-read.csv(file=paste(path_data,"data_NPP.csv",sep=""),header=TRUE) # g C year-1 m-2
data_NPP$NPP_gC_m2_yr<-1e3*data_NPP$NPP_gC_m2_yr/year_to_day # mg C day-1 m-2
data_NPP$exsudate<-0.1*data_NPP$NPP_gC_m2_yr

### INITIAL DENSITIES ### ----
N=63*month_to_day
DOC=I_DOC*month_to_day
FOM=I_FOM*month_to_day
micro=exp(6) # g m-2 (fresh biomass)
micro=micro*1e3*fresh_to_dry*dry_to_C
meso=exp(4) # g m-2 (fresh biomass)
meso=meso*1e3*fresh_to_dry*dry_to_C
macro=exp(5) # g m-2 (fresh biomass)
macro=macro*1e3*fresh_to_dry*dry_to_C

### CN RATIOS ###----
data_CN<-read.csv(file=paste(path_data,"data_stoichio_Andrieux_2021.csv",sep=""),header=TRUE) # data from Andrieux et al. 2021 (10.1111/geb.13265)
# microbes
data_CN_microbes<-data_CN[data_CN$Group=="Microbe",]
ggplot(data=data_CN_microbes)+
  geom_point(aes(Diet_full,CN_ratio,colour=Diet_full))
data_CN_microbes<-data_CN_microbes[which(is.na(data_CN_microbes$CN_ratio)==F)
                                   & data_CN_microbes$CN_ratio<20,which(names(data_CN_microbes)%in%c("Diet_full","CN_ratio"))]
summaryBy(CN_ratio~Diet_full, data=data_CN_microbes, FUN=c(mean, sd, min, max))
# invertebrates
data_CN_invert<-data_CN[which(data_CN$Group=="Invert"
                        & data_CN$Habitat=="Terrestrial"
                        & data_CN$Diet_full%in%c("C","D")),]
ggplot(data=data_CN_invert)+
  geom_point(aes(Diet_full,CN_ratio,colour=Diet_full))
data_CN_invert<-data_CN_invert[which(is.na(data_CN_invert$CN_ratio)==F),which(names(data_CN_invert)%in%c("Diet_full","CN_ratio"))]
CN_params<-summaryBy(CN_ratio~Diet_full, data=data_CN_invert, FUN=c(mean, sd, min, max))
CN_params$max=CN_params$CN_ratio.mean*2-CN_params$CN_ratio.min
# other data set for bacteria and fungi
data_CN<-read.csv(file=paste(path_data,"data_stoichio_Mouginot_2014.csv",sep=""),header=TRUE) # data from Mouginot et al. 2014 (10.1016/j.soilbio.2014.05.011)
data_CN<-data_CN %>% separate(C.N,c("CN_ratio",NA),sep="±")
data_CN<-data_CN %>% separate(C.P,c("CP_ratio",NA),sep="±")
data_CN<-data_CN %>% separate(N.P,c("NP_ratio",NA),sep="±")
data_CN$CN_ratio<-as.numeric(data_CN$CN_ratio)
data_CN$CP_ratio<-as.numeric(data_CN$CP_ratio)
data_CN$NP_ratio<-as.numeric(data_CN$NP_ratio)
ggplot(data=data_CN)+
  geom_point(aes(Realm,CN_ratio,colour=Realm))
# separation between bacteria and fungi
summary<-summaryBy(CN_ratio~Realm, data=data_CN, FUN=c(mean, sd, min, max))
ggplot(data=data.frame(x=c(min(summary$CN_ratio.min),max(summary$CN_ratio.max))), aes(x)) +
  stat_function(fun = dnorm, n = 500, args = list(mean = summary$CN_ratio.mean[summary$Realm=="Bacteria"], sd = summary$CN_ratio.sd[summary$Realm=="Bacteria"]),colour="purple")+
  stat_function(fun = dnorm, n = 500, args = list(mean = summary$CN_ratio.mean[summary$Realm=="Fungi"], sd = summary$CN_ratio.sd[summary$Realm=="Fungi"]),colour="orange")
# all microbes together
mean_CN=mean(data_CN$CN_ratio)
sd_CN=sd(data_CN$CN_ratio)
min_CN=min(data_CN$CN_ratio)
max_CN=2*mean_CN-min_CN
### PREFRENCE FOR PREY:PREDATOR BODY MASS RATIO ### ----
pref_M<-function(x,theta,sigma){
  return(exp(-((x-theta)/sigma)^2))
}
theta=-2
data<-expand.grid(x=seq(-5,1,0.01),sigma=c(0.1,0.5,1,2,5))
data$pref<-pref_M(data$x,theta,data$sigma)
data$sigma<-as.factor(data$sigma)
text<-data.frame(x=c(-1.7,-1.2,-0.5,-0.5,-0.5),sigma=c(0.1,0.5,1,2,5))
text$pref<-pref_M(text$x,theta,text$sigma)
text$sigma<-as.factor(text$sigma)

p1<-ggplot(data=data)+
  geom_line(aes(x,pref,colour=sigma),linewidth=1.5)+
  geom_vline(xintercept=-2,linetype="dashed")+
  geom_vline(xintercept=0,linetype="dashed")+
  geom_text(data=text,aes(x+0.1,pref+0.1,label=sigma,colour=sigma),size=7,family="serif",show.legend = FALSE)+
  annotate(geom="text",label="italic(\u03B8[i]^M)",x=-1.7,y=1.07,size=8,family="serif",parse=TRUE)+
  annotate(geom="text",label="Predator\nbody mass",x=0.7,y=1.06,size=6,family="serif")+
  scale_colour_viridis(discrete=T,
                       name=expression(italic("\u03C3"[i]^"M")))+
  scale_x_continuous(breaks=seq(min(data$x),max(data$x),1))+
  scale_y_continuous(breaks=seq(0,1,0.25),
                     limits=c(0,1.1))+
  theme+
  xlab("Log of prey to predator body mass ratio")+
  ylab("Preference for the prey")

ggsave(paste(path_figure,"methods_pref.pdf",sep=""),p1, width = 8, height = 4, device = cairo_pdf)

### GROWTH FUNCTIONS ### ----
fun<-function(h){
  return(exp(h)/(exp(h)-1))
}
G_rel=0.5
h1=log(1/(1-G_rel))

data<-data.frame(h=seq(0,3,0.01))
data$fun<-fun(data$h)*G_rel
for(i in 1:dim(data)[1]){
  data$fun[i]=min(data$fun[i],1)
}

p1<-ggplot(data=data)+
      geom_line(aes(h,fun),size=1.5)+
      geom_vline(xintercept=h1,linetype="dashed")+
      geom_hline(yintercept=G_rel,linetype="dashed")+
      annotate(geom="text",x=2.7,y=G_rel-0.2,label="italic(frac(G[k],G[k]+B[k]))",size=8,family="serif",parse=TRUE)+
      theme+
      scale_x_continuous(limits=c(0,3),
                         breaks=c(0,h1,1,2,3),
                         labels=c(0,expression(italic(h[1]),1,2,3)))+
      scale_y_continuous(limits=c(0,1.5),
                         breaks=c(0,1))+
      xlab(expression("Body mass interval width "*italic(h)))+
      ylab(expression("Fraction of biomass entering class "*italic(k)+1))
ggsave(paste(path_figure,"figure_growth_frac.pdf",sep=""),p1, width = 7, height = 6, device = cairo_pdf)

### RESPIRATION - INGESTION BALANCE ### ----
# ingestion detritivore
T0=288.15 # 15°C
s_h=-0.165 # appropriate scaling
h0_1=2.07e2 # day g -s_h
h0_2=5e2 # day g -s_h
h0=(h0_1+h0_2)/2
E_h_1=-1.36 # eV
E_h_2=-0.86 # eV
E_h_detri=(E_h_1+E_h_2)/2
s_h_detri=s_h+1
h0_detri=h0*exp(E_h_detri/(kB*T0))

# ingestion carnivores
T0=293.15 # reference temperature 20°C
data_FR<-read.csv(file=paste(path_data,"FR_li_2017_oikos.csv",sep=""),header=TRUE)
data_FR<-data_FR[,which(names(data_FR)%in%c("dimensionality","experimental.duration.seconds","temperature.degree.celcius","predator.mass.mg","prey.mass.mg","Starvation.Y.N","attack.rate","handling.time","original.publication.short","arena.size"))]
names(data_FR)=c("publication","dimensionality","arena_size","duration","temperature","predator_mass","prey_mass","attack_rate","handling_time","starvation")
data_FR$temperature<-celsius_to_kelvin+data_FR$temperature
data_FR$handling_time<-second_to_day*data_FR$handling_time*data_FR$predator_mass # conversion from individual based (s Ind_pred mg_prey-1) handling time to biomass based
data_FR$predator_mass<-fresh_to_dry*dry_to_C*data_FR$predator_mass # already in mg
data_FR$duration<-second_to_day*data_FR$duration
data_FR$temperature<--1/(data_FR$temperature*kB)

model=lme(data=data_FR,
          log(handling_time)~log(predator_mass)+temperature,
          random = ~ 1| publication/arena_size,
          method = "ML",
          control=lmeControl(opt="optim",
                             msMaxIter = 200,
                             msMaxEval = 500,
                             msVerbose = F))
params<-summary(model)
params<-(params$coefficients)[[1]]
h0_pred=exp(params[1])
s_h_pred=params[2]
E_h_pred=params[3]

# respiration
joule_to_carbon=1/20.1*0.5363
data_metabo<-read.csv(file=paste(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaMetabolicData.csv",sep=""),
                      header=TRUE)
names(data_metabo)<-c("study","trophic_group","taxon","metabolic_rate","bodymass","temperature")
data_metabo$temperature<--1/(celsius_to_kelvin+data_metabo$temperature)*1/kB
data_metabo$metabolic_rate<-data_metabo$metabolic_rate*joule_to_carbon/hour_to_day
data_metabo$bodymass<-data_metabo$bodymass*fresh_to_dry*dry_to_C
data_metabo$metabolic_rate<-data_metabo$metabolic_rate/data_metabo$bodymass # conversion into mass metabolic rate
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~trophic_group+log(bodymass):trophic_group+temperature:trophic_group)
params<-summary(model)
params<-(params$coefficients)[,1]
names(params)[1]="trophic_groupMacrofauna"
params[2]=params[2]+params[1] # total effect for mesofauna
params[3]=params[3]+params[1] # total effect for microbes
params<-as.data.frame(params)
params$names<-row.names(params)
params<-params %>% separate(names,c("trophic_group","parameter"),sep=":",fill="right")
params<-params %>% separate(trophic_group,c(NA,"trophic_group"),sep=13)
params$parameter[is.na(params$parameter)]="R0"
params$parameter<-as.factor(params$parameter)
levels(params$parameter)=c("s_R","R0","E_R")
params<-dcast(params,trophic_group~parameter,value.var="params")

# balance
data<-data.frame(bodymass_fresh=10^seq(-9,4,1))
data$bodymass=fresh_to_dry*dry_to_C*data$bodymass_fresh
data$Microbe<-allometric_parameter(data$bodymass,1,exp(params$R0[3]),params$s_R[3],params$E_R[3],kB,293.15)
data$Mesofauna<-allometric_parameter(data$bodymass,1,exp(params$R0[2]),params$s_R[2],params$E_R[2],kB,293.15)
data$Macrofauna<-allometric_parameter(data$bodymass,1,exp(params$R0[1]),params$s_R[1],params$E_R[1],kB,293.15)
data<-melt(data, id.vars = c("bodymass_fresh","bodymass"),
           variable.name = "trophic_group", 
           value.name = "respiration")
data$h_detri<-1/allometric_parameter(data$bodymass,1,h0_detri,s_h_detri,E_h_detri,kB,293.15)
data$h_pred<-1/allometric_parameter(data$bodymass,1,h0_pred,s_h_pred,E_h_pred,kB,293.15)
data$trophic_group<-factor(data$trophic_group,levels=c("Macrofauna","Mesofauna","Microbe"))

p1<-ggplot(data=data)+
  geom_line(aes(bodymass_fresh,respiration,colour=trophic_group),size=1.5)+
  geom_line(aes(bodymass_fresh,h_detri),linetype="dashed")+
  geom_line(aes(bodymass_fresh,h_pred),linetype="dotted")+
  annotate(geom="text",x=1e-3,y=1e4,label="Detritivore maximum\ningestion rate",family="serif",size=6)+
  annotate(geom="text",x=1e2,y=1e1,label="Carnivore maximum\ningestion rate",family="serif",size=6)+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE))+
  theme+theme(legend.title=element_blank())+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Fresh body mass (mg)"))+
  ylab(expression("Respiration rate (mgC mgC"^-1*" day"^-1*")"))
ggsave(paste(path_figure,"supp_balance.pdf",sep=""),p1, width = 9, height = 5.5, device = cairo_pdf)

### ASSIMILATION EFFICIENCY ### ----
data_lang<-read.csv(file=paste(path_data,"Assimilation efficiency Lang et al 2017 Oikos/Lang_et_al_2017_data.csv",sep=""),sep=",",header=TRUE)
data_lang<-data_lang[,c("consumer.type","assimilation.efficiency")]
data_lang<-data_lang[data_lang$assimilation.efficiency>-999,]
data_lang<-aggregate(data=data_lang,assimilation.efficiency~consumer.type, FUN = function(x) c(mean = mean(x), se = sd(x)/length(x)))

                     
### DORMANCT FUNCTION FIGURE ### ----
activation<-function(G,R){
  return(G/(G+R))
}
G=seq(0,1,length.out=30)
R=0.5

data<-data.frame(G=G,
                 R=R,
                 activation=activation(G,R))
data$dormancy<-1-data$activation
data<-melt(data, id.vars = c("G","R"),
           variable.name = "fun", 
           value.name = "value")
levels(data$fun)<-c(expression("activation "*italic(q[act])),expression("dormancy "*italic(q[dorm])))

p1<-ggplot(data=data)+
      geom_line(aes(G,value,colour=fun),linewidth=2)+
      geom_vline(xintercept=0.5,linetype="dashed")+
      annotate(geom="text",label="dormant",x=0.3,y=0.85,family="serif",size=7)+
      annotate(geom="text",label="active",x=0.7,y=0.85,family="serif",size=7)+
      theme+theme(legend.title=element_blank())+
      scale_colour_manual(values=c("blue","red"),
                          labels = parse_format())+
      guides(colour = guide_legend(label.hjust = 0))+
      scale_x_continuous(breaks=c(0,0.5),
                         labels=c(0,"Loss (metabolism)"),
                         name="Gain (resource consumption)")+
      ylab("Sensitivity")
ggsave(paste(path_figure,"methods_function_dormancy.pdf",sep=""),p1, width = 7, height = 4, device = cairo_pdf)

### BIOMASS VARIATION FUNCTION FIGURE ### ----
x=seq(1,7,length.out=50)
data<-data.frame(x=x,
                 y=(1-exp(-x)))

p1<-ggplot(data=data)+
  geom_line(aes(x,y),linewidth=2)+
  geom_hline(yintercept=1,linetype="dashed")+
  #annotate(geom="text",label=as.character(expression(italic(B[i]*"(t)"))),x=2,y=1.02,parse=T,family="serif",size=10)+
  theme+theme(axis.ticks=element_blank(),
              axis.text=element_blank(),
              panel.grid=element_blank(),
              plot.title=element_text(hjust = 0,size=30),
              text = element_text(family="serif",size=30))+
  xlab(expression(italic("\u0394"*"t")))+
  ylab(expression(italic("\u0394"*"X"[i])))+
  ggtitle(expression(italic(B[i]*"(t)")))
ggsave(paste(path_figure,"methods_integration.pdf",sep=""),p1, width = 4, height = 3, device = cairo_pdf)

### FAECES DECOMPITION AND SPECIFIC AREA ### ----
data_faeces<-read.csv(file=paste0(path_data,"Faeces Joly 2018/Phys_Chem_characteristics.csv"),header=TRUE)
data_faeces<-data_faeces[,c("Tree.species","Form","Replicate","Litter.and.faecal.pellet.specific.area","Faeces.particles.specific.area")]
names(data_faeces)<-c("tree","form","replicate","litter_specific_area","particle_specific_area")
data_faeces$specific_area<-data_faeces$litter_specific_area
data_faeces$specific_area[data_faeces$form=="faeces"]<-data_faeces$particle_specific_area[data_faeces$form=="faeces"]

temp<-read.csv(file=paste0(path_data,"Faeces Joly 2018/CN_dynamics.csv"),header=TRUE)
temp<-temp[,c("Tree.species","Form","Replicate","C.loss....")]
names(temp)<-c("tree","form","replicate","C_loss")
data_faeces<-merge(data_faeces,temp,by=c("tree","form","replicate"),all = TRUE)

# ratios
data<-data_faeces[data_faeces$form=="faeces",]
data$specific_area<-data$specific_area/data_faeces$specific_area[data_faeces$form=="litter"]
data$C_loss<-data$C_loss/data_faeces$C_loss[data_faeces$form=="litter"]
data<-data[is.na(data$specific_area)==FALSE,]

ggplot(data)+
  geom_point(aes(specific_area,C_loss,colour=tree),size=2)+
  theme

data_faeces<-read.csv(file=paste0(path_data,"Faeces Joly 2020/faeces_litter_quality.csv"),header=TRUE)
temp<-read.csv(file=paste0(path_data,"Faeces Joly 2020/faeces_litter_decomposition.csv"),header=TRUE)
data_faeces<-merge(data_faeces,temp,by=c("litter","animal","type","rep"))

data<-data_faeces[data_faeces$type=="faeces",]
data$specific_area_particles<-data$specific_area_particles/data_faeces$specific_area_particles[data_faeces$type=="litter"]
data$C_loss<-data$C_loss/data_faeces$C_loss[data_faeces$type=="litter"]

ggplot(data)+
  geom_smooth(aes(specific_area_particles,C_loss),method = "lm")+
  geom_point(aes(specific_area_particles,C_loss,colour=litter),size=2)+
  theme


### CONTRIBUTION OF INVERTEBRATES TO LITTER DECOMPOSITION ### ----
data_wall<-read.csv(file=paste(path_data,"data_wall_2008.csv",sep=""),sep=",",header=TRUE)
data<-data_wall[data_wall$climat%in%c("temperate","wet_tropical"),]
data<-dcast(data,climate+site~treatment,value.var="decomposition_rate_per_day")
data$ratio<-(1-data$naph/data$no_naph)*100
mean(data$ratio)
sd(data$ratio)

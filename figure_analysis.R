##################
### LIBRARIRES ### ----
##################

library(rstudioapi) # to set the working directory
library(reshape2) # to format data frames
library(ggpubr)
library(ggrepel) # to use geom_repel (labels)
library(cowplot) # to arrange graphs
library(tidyr) # for the function separate
library(grid) # for rectGrob
#library(gridExtra) # for grid.arrange
library(car) # for Anova

source("figure_settings.R")

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
path_potapov<-"Data/Potapov et al 2021 Ecology/DataS1/"
path_hedenec="Data/Hedenec 2022/"
path_figures<-"Figures/"
path_data<-"Data/"

############# ----
# MAIN TEXT # ----
############# ----
# load data # ----
path_results<-"results_final_main/"#"results_23_02_22/"#"results_23_02_03/"#"results_23_01_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

data_size_structured<-Rload_data(path_results,3)
palette_size_structured<-Rmake_palettes(data_size_structured$species_traits,
                                        data_size_structured$n_tot_species,
                                        data_size_structured$n_faeces)

data_detritivore_spectrum<-Rload_data(path_results,2)
palette_detritivore_spectrum<-Rmake_palettes(data_detritivore_spectrum$species_traits,
                                             data_detritivore_spectrum$n_tot_species,
                                             data_detritivore_spectrum$n_faeces)
# species data
data_species<-rbind(data_multi_channel$data_species,
                    data_size_structured$data_species)
data_species$trophic_group_file<-as.factor(data_species$trophic_group_file)
levels(data_species$trophic_group_file)=c(label_multichannel,label_size_structured)
# detritus data
data_detritus<-rbind(data_multi_channel$data_detritus,
                     data_size_structured$data_detritus)
data_detritus$trophic_group_file<-as.factor(data_detritus$trophic_group_file)
levels(data_detritus$trophic_group_file)=c(label_multichannel,label_size_structured)
# detritus production data
production<-rbind(data_multi_channel$data_detritus_production,
                  data_size_structured$data_detritus_production)
production$trophic_group_file<-as.factor(production$trophic_group_file)
levels(production$trophic_group_file)=c(label_multichannel,label_size_structured)

# temp<-Rload_data(path_results,trophic_group_file)
# list2env(temp,.GlobalEnv)
# temp<-Rmake_palettes(species_traits,n_faeces)
# list2env(temp,.GlobalEnv)

# FOOD WEB STRUCTURE # ----
# Trophic levels and diet # ----
# trophic levels
data<-data_species
data$TL[data$TL==0]=NA

data_potapov<-read.table(paste(path_potapov,"RawData_MassDensityTrophicLevel.csv",sep=""),sep=",",header=T)
data_potapov<-data_potapov[,c("SizeClass","TL")]
data_potapov$TL<-data_potapov$TL-1
data_potapov<-data_potapov %>% separate(SizeClass, c("bodymass",NA), sep = "-")
data_potapov$bodymass[data_potapov$bodymass==">1000"]=1000
data_potapov$bodymass<-as.numeric(data_potapov$bodymass)
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")

p1<-ggplot(data=data)+
      geom_point(aes(bodymass,TL,colour=trophic_group),size=4)+
      palette_multi_channel$scale_colour_line_trophic_groups+
      theme+theme(legend.position="bottom",
                  legend.key.size=unit(30,'pt'))+
      guides(colour = guide_legend(reverse=FALSE))
legend<-get_legend(p1)

p1<-ggplot(data=data)+
      geom_smooth(aes(bodymass,TL), colour="black",fill="lightgrey", linetype="blank")+
      geom_point(aes(bodymass,TL,colour=trophic_group),size=4)+
      geom_smooth(data=data_potapov, aes(bodymass,TL), se=FALSE, colour="black", linetype="dashed")+
      geom_smooth(aes(bodymass,TL), colour="black", se=FALSE, linetype="solid")+
      facet_wrap(~trophic_group_file,nrow=2)+
      theme+theme(legend.position="none")+
      palette_multi_channel$scale_colour_line_trophic_groups+
      scale_x_log10(minor_breaks = 10^seq(-10,5,1),
                    breaks = 10^(seq(-8,2,2)),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      xlab(label_bodymass)+
      ylab(label_TL)

# diet
data<-data_species[data_species$trophic_group_file==label_multichannel,]
diet<-read.table(paste(path_results,"diet_",data$simu_ID[1],".txt",sep=""),sep=",",header=T)
diet$pred<-data_multi_channel$names_species # add the names of the predators
diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
data<-data[order(data$names),]
diet[data$TL==0,c(1:(ncol(diet)-1))]=NA # removes the diet of extinct carnivores
diet<-diet[which(data$trophic_type%in%c("microbivores","carnivores")),c(1:data_multi_channel$n_tot_species,ncol(diet))] # keeps only consumers
diet<-Rmake_table_for_plot(diet,"pred","prey","flow") # makes diet a long table
levels(diet$prey)<-data_multi_channel$names_species
survivor<-data$names[data$TL>0 & data$trophic_type%in%c("microbivores","carnivores")]
diet<-diet[diet$pred%in%survivor,]
# sum of prey biomass by trophic group
temp<-data_species[,c("names","trophic_group")]
names(temp)[1]="prey"
diet<-merge(diet,temp,by="prey")
diet$prey<-NULL
names(diet)[3]="prey"
names(temp)[1]="pred"
diet<-merge(diet,temp,by="pred")
diet<-aggregate(flow~pred+prey+trophic_group, data=diet, sum)

n_micro<-length(which(data$names[data$trophic_group%in%c("micro-food web microbivores","micro-food web carnivores")]%in%survivor))
n_macro<-length(which(data$names[data$trophic_group%in%c("macro-food web carnivores")]%in%survivor))

p2<-ggplot(data=diet)+
      geom_bar(aes(pred,flow,fill=prey),stat="identity",position="fill")+
      geom_point(aes(pred,-0.05,colour=trophic_group),size=4)+
      annotate("errorbarh",y=1.03,xmin=0.7,xmax=n_micro+0.3,height = 0.05)+
      annotate("errorbarh",y=1.03,xmin=n_micro+0.7,xmax=n_micro+n_macro+0.45,height = 0.05)+
      annotate("text",x=n_micro/2,y=1.08,label="micro-food web",family="serif",size=6,angle=-90)+
      annotate("text",x=n_micro+n_macro/2,y=1.08,label="macro-food web",family="serif",size=6,angle=-90)+
      theme+theme(legend.position="none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.line.y = element_line(arrow = arrow(angle = 15, length = unit(.3,"inches"),type = "closed")),
                  #axis.line.y=element_blank(),
                  panel.grid=element_blank())+
      palette_multi_channel$scale_colour_line_trophic_groups_consumers+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      coord_flip()+
      xlab("Increasing body\nsize of predators")+
      ylab(label_diet)

legend_height=0.15
graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 1+legend_height)) +
  draw_plot(p1, 0, legend_height, 1, 1)+
  draw_plot(p2, 1, legend_height, 1, 1)+
  draw_plot(legend, 0, 0, 2, legend_height)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1)+legend_height, size = 30)
ggsave(paste(path_figures,"figure_TL_diet.pdf",sep=""), graph, width = 12, height = 8, device=cairo_pdf)

# ECOSYSTEM FUNCTIONING # ----
# Decomposition and respiration # ----
# detritus stocks
data<-aggregate(biomass~type+input_FOM+input_DOC+trophic_group_file, data=data_detritus, sum)
p1<-ggplot(data)+
  #geom_col(aes(type,biomass,fill=type))+
  geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(biomass)*0.5,ymax=biomass,fill=type))+
  facet_wrap(~trophic_group_file)+
  theme+theme(legend.position="none",
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank())+
  palette_multi_channel$scale_colour_bar_detritus_classes+
  scale_x_continuous(breaks=1:length(levels(data$type)),
                     labels=levels(data$type))+
  y_axis_log10+
  ylab(label_stock)

data<-rbind(data_multi_channel$data_detritus,
            data_size_structured$data_detritus)
data$trophic_group_file<-as.factor(data$trophic_group_file)
levels(data$trophic_group_file)=c(label_multichannel,label_size_structured)
data<-aggregate(biomass~type+trophic_group_file, data=data, sum)
data<-data[data$type=="FOM",]
data$biomass[data$trophic_group_file==label_multichannel]/data$biomass[data$trophic_group_file==label_size_structured]

# decomposition by each decomposer trophic group
consumption=NULL
# multichannel
temp<-read.table(paste(path_results,"consumption_",1,".txt",sep=""),sep=",",header=T)
names(temp)<-data_multi_channel$names_tot
temp<-temp[,data_multi_channel$names_detritus] # selects abiotic resources only
temp$consumer<-data_multi_channel$names_species
temp$trophic_group_file<-label_multichannel
temp<-Rmake_table_for_plot(temp,c("consumer","trophic_group_file"),"resource","biomass") # makes diet a long table
consumption<-rbind(consumption,temp)
# size-structured
temp<-read.table(paste(path_results,"consumption_",3,".txt",sep=""),sep=",",header=T)
names(temp)<-data_size_structured$names_tot
temp<-temp[,data_size_structured$names_detritus] # selects abiotic resources only
temp$consumer<-data_size_structured$names_species
temp$trophic_group_file<-label_size_structured
temp<-Rmake_table_for_plot(temp,c("consumer","trophic_group_file"),"resource","biomass") # makes diet a long table
consumption<-rbind(consumption,temp)
# resource type
temp<-data_multi_channel$detritus_traits[,c("type","names")]
names(temp)<-c("type","resource")
consumption<-merge(consumption,temp,by="resource")
# consumer trophic group
temp<-data_multi_channel$species_traits[,c("trophic_group","names")]
names(temp)<-c("trophic_group","consumer")
consumption<-merge(consumption,temp,by="consumer")
rm(temp)
# only consider decomposers
consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores","trophic whales"),]
consumption<-aggregate(biomass~type+trophic_group+trophic_group_file, data=consumption, sum)
consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores","trophic whales"))
levels(consumption$trophic_group)<-c("microbes","macro-detritivores","trophic whales")
consumption$type<-factor(consumption$type,levels=c("faeces","FOM","SOM","DOC"))

p2<-ggplot(data=consumption)+
  geom_bar(aes(type,biomass,fill=trophic_group),stat="identity",position=position_stack(reverse = TRUE))+
  facet_wrap(~trophic_group_file)+
  theme+theme(legend.position="none",
              #legend.position=c(1, 1),
              #legend.justification=c(1, 1),
              axis.title.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank())+
  palette_multi_channel$scale_colour_bar_trophic_groups_detritivores+
  ylab(label_decomposition)

# detritus production
data<-aggregate(production~type+trophic_group_file, data=production, sum)
data$type[data$type%in%c("FOM","DOC")]<-NA
data<-data[data$type%in%c("faeces","SOM"),]
data$type<-droplevels(data$type)
#production$production[production$production==0]=NA

p3<-ggplot(data)+
  geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(production)*0.5,ymax=production,fill=type))+
  facet_wrap(~trophic_group_file)+
  theme+theme(legend.position="none",
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank())+
  palette_multi_channel$scale_colour_bar_detritus_production+
  scale_x_continuous(breaks=1:length(levels(data$type)),
                     labels=levels(data$type))+
  #y_axis_log10+
  ylab(label_production)

# respiration
data<-aggregate(respiration~trophic_group+input_FOM+input_DOC+trophic_group_file, data=data_species, sum)
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")

p4<-ggplot(data)+
  geom_bar(aes(1,respiration,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+
  palette_multi_channel$scale_colour_bar_trophic_groups
legend_1<-get_legend(p4)

p4<-ggplot(data)+
  geom_bar(aes(1,respiration,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  facet_wrap(~trophic_group_file)+
  theme+theme(legend.position="none",
              axis.title.x=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank())+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  scale_y_continuous(labels = scales::percent)+
  ylab(label_respiration_fraction)

data$relative<-0
data$relative[data$trophic_group_file==label_multichannel]<-data$respiration[data$trophic_group_file==label_multichannel]/sum(data$respiration[data$trophic_group_file==label_multichannel])
data$relative[data$trophic_group_file==label_size_structured]<-data$respiration[data$trophic_group_file==label_size_structured]/sum(data$respiration[data$trophic_group_file==label_size_structured])
data$relative<-round(data$relative,3)*100

# Respiration distribution empiric # ----
### Calculation of metabolic parameters
# data from Johnston and Sibly (2018) (10.1038/s41559-018-0648-6)
joule_to_carbon=1/20.1*0.5363
hour_to_day=1/24
celsius_to_kelvin=273.15
kB=8.62e-5 # Bolzmann's constant
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
# metabolic rate : J.h-1
# body mass : mg fresh weight
data_metabo<-read.csv(file=paste(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaMetabolicData.csv",sep=""), header=TRUE)
names(data_metabo)<-c("study","Fauna_group","taxon","metabolic_rate","bodymass","temperature")

# conversion
data_metabo$temperature_arrhenius<--1/(celsius_to_kelvin+data_metabo$temperature)*1/kB
data_metabo$metabolic_rate<-data_metabo$metabolic_rate*joule_to_carbon/hour_to_day
data_metabo$bodymass<-data_metabo$bodymass*fresh_to_dry*dry_to_C
data_metabo$metabolic_rate<-data_metabo$metabolic_rate/data_metabo$bodymass # conversion into mass metabolic rate
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~Fauna_group+log(bodymass):Fauna_group+temperature_arrhenius:Fauna_group)
params<-summary(model)
params
Anova(model,type=2,test.statisticf="F")
params<-(params$coefficients)[,1]
names(params)[1]="Fauna_groupMacrofauna"
params[2]=params[2]+params[1] # total effect for mesofauna
params[3]=params[3]+params[1] # total effect for microbes
params<-as.data.frame(params)
params$names<-row.names(params)
params<-params %>% separate(names,c("Fauna_group","parameter"),sep=":",fill="right")
params<-params %>% separate(Fauna_group,c(NA,"Fauna_group"),sep=11)
params$parameter[is.na(params$parameter)]="R0"
params$parameter<-as.factor(params$parameter)
levels(params$parameter)=c("s_R","R0","E_R")
params<-dcast(params,Fauna_group~parameter,value.var="params")

### Calculation of the respiration of each sub-food web
data_abundance<-read.table(paste0(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaAbundanceData.csv"),header=T,sep=",")
data_bodymass<-read.table(paste0(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv"),header=T,sep=",")
Mkm2_to_m2=1e6*1e6 # million km² to km²
Pg_to_mg=1e18 # petagram to milligram

# mean body mass of the main soil taxa
data_bodymass$bodymass<-data_bodymass$Individual_bodymass_mg_DM*dry_to_C # body
data_bodymass$bodymass_fresh<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry
data_bodymass$trophic_group<-"macro-food web"
data_bodymass$trophic_group[data_bodymass$bodymass_fresh<0.001]="micro-food web"
data_bodymass$trophic_group[data_bodymass$Soil_biota_group=="Bacteria"]<-"microbes"
data_bodymass$Individual_bodymass_mg_DM<-NULL
# abundance of the main soil fauna taxa
data_abundance<-data_abundance[,c("Biome_ecosystem","Soil_biota_group","Biomass_g_m2")]
data_abundance$biomass_density<-data_abundance$Biomass_g_m2*1e3*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance$Biomass_g_m2<-NULL
# abundance of soil microbes
microbes<-read.table("Data/Metabolic_rates_Johnston_2018/Xu_2013_Global_Ecology_and_Biogeography_aggregated_data.csv",header=T,sep=",") # data from the main text of Xu et al 2013 Global Ecology and Biogeography
microbes$biomass_density<-microbes$soil_microbial_C_Pg_0_30_cm*Pg_to_mg/(microbes$Area_Mkm2*Mkm2_to_m2)
microbes$Soil_biota_group="Bacteria"
microbes<-microbes[,names(data_abundance)]

# final data set
data<-rbind(data_abundance,microbes)
rm(data_abundance,microbes)
data<-merge(data,data_bodymass,by=c("Soil_biota_group"))
data$trophic_group<-factor(data$trophic_group,levels=c("microbes","micro-food web","macro-food web"))
data<-data[data$Biome_ecosystem%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data$Biome_ecosystem<-factor(data$Biome_ecosystem,levels=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"))
levels(data$Biome_ecosystem)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")
data<-merge(data,params,by=c("Fauna_group")) # add the metabolic parameters

# respiration
data$metabolism<-exp(data$R0)*(data$bodymass^data$s_R)*exp(-data$E_R/(kB*(celsius_to_kelvin+15)))
data$respiration<-data$metabolism*data$biomass_density
data<-aggregate(respiration~Biome_ecosystem+trophic_group+Soil_biota_group,data,mean)
data<-aggregate(respiration~trophic_group,data,sum)
data$trophic_group_file="Empirical data"

p5<-ggplot(data=data)+
  geom_col(aes(1,respiration,fill=trophic_group),position=position_fill(reverse=TRUE))+
  theme+theme(legend.title=element_blank())+
  scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                    guide = guide_legend(reverse=TRUE))
legend_2<-get_legend(p5)

p5<-ggplot(data=data)+
      geom_col(aes(1,respiration,fill=trophic_group),position=position_fill(reverse=TRUE))+
      facet_wrap(~trophic_group_file)+
      theme+theme(legend.position="none",
                  legend.title=element_blank(),
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                        guide = guide_legend(reverse=TRUE))+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      xlab("Biome")+
      ylab(label_respiration_fraction)

data$percent<-data$respiration/sum(data$respiration)*100

### Final graph #----
h_logo=0.1
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 2+2*h_logo)) +
  draw_plot(p1, 0, 1+2*h_logo, 1.5, 1)+
  draw_plot(p2, 1.5, 1+2*h_logo, 1.5, 1)+
  draw_plot(p3, 0, h_logo, 0.9, 1)+
  draw_plot(p4, 0.9, h_logo, 1, 1)+
  draw_plot(legend_2, 1.9, h_logo, 0.5, 0.4)+
  draw_plot(legend_1, 1.9, 0.4+h_logo, 0.5, 0.4)+
  draw_plot(p5, 2.4, h_logo, 0.6, 1)+
  draw_plot_label(c("A","B","C","D","E"), c(0,1.5,0,0.9,2.3), c(2,2,1,1,1)*(1+h_logo), size = 30)
ggsave(paste(path_figures,"figure_detritus_respiration.pdf",sep=""), graph, width = 14, height = 10, device=cairo_pdf)

# graph<-ggdraw(xlim = c(0, 2.7), ylim = c(0, 2.15)) +
#   draw_plot(p1, 0, 1.15, 1, 1)+
#   draw_plot(p2, 1.05, 1.15, 1, 1)+
#   draw_plot(p3, 2.1, 1.15, 0.6, 1)+
#   draw_plot(p4, 0, 0, 1.5, 1)+
#   draw_plot(p5, 1.5, 0, 1.1, 1)+
#   draw_plot_label(c("A","B","C","D","E"), c(0,0.95,2,0,1.4), c(2.15,2.15,2.15,1.02,1.02), size = 30)
# ggsave(paste(path_figures,"figure_detritus_respiration.pdf",sep=""), graph, width = 12, height = 9, device=cairo_pdf)

# Stock ratios (Table 1) # ----
data<-aggregate(biomass~trophic_group+trophic_group_file, data=data_species, sum)
data<-data[data$trophic_group=="microbes",]
temp<-aggregate(biomass~type+trophic_group_file, data=rbind(data_multi_channel$data_detritus,
                                                            data_size_structured$data_detritus), sum)
temp$trophic_group_file<-as.factor(temp$trophic_group_file)
levels(temp$trophic_group_file)=c(label_multichannel,label_size_structured)
names(temp)[1]="trophic_group"
data<-rbind(data,temp)
rm(temp)
data<-dcast(data, trophic_group_file ~ trophic_group, value.var="biomass")
data$ratio_C<-data$microbes/(data$SOM+data$microbes)*100
data$ratio_detritus<-data$FOM/(data$SOM+data$FOM)*100
ratio<-data[,c(2:ncol(data))]
ratio[1,]<-ratio[1,]/ratio[2,]

# total soil C according to Xu et al. 2013, total value in table 4
area=128.3e6 # km²
area=area*1e6 # conversion to m²
C_mic=16.72 # Pg C
C_mic=C_mic*1e18 # conversion to mg
C_SOM=C_mic/area*(1/0.012-1)

data_FOM<-read.table(paste0(path_data,"Hedenec 2022/litter.csv"),header=T,sep=",") # data of FOM and DOC inputs in the main land ecosystems
data_FOM<-data_FOM[,c(1,4)]
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
data_FOM$litter_stock<-data_FOM$litter_stock_kg_ha*dry_to_C*1e6/1e4 # conversion to mg C per m²
C_FOM=mean(data_FOM$litter_stock)

C_FOM/(C_SOM+C_FOM)*100

# PREDICTION OF BIOMASS DISTRIBUTION # ----
# Biomass distribution model # ----
data<-data_species
data$biomass[data$biomass<1e-7]=NA
data<-aggregate(biomass~trophic_group+input_FOM+input_DOC+trophic_group_file, data=data, sum)
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
# raw biomass per trophic group
p1<-ggplot(data)+
      #geom_bar(aes(trophic_group,biomass,fill=trophic_group),stat="identity")+
      geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
      facet_wrap(~trophic_group_file)+
      theme+theme(#legend.position="bottom",
                  #legend.direction="vertical",
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  panel.grid.minor.y=element_blank(),
                  panel.grid.major.y=element_blank())+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      #guides(fill=guide_legend(nrow=3,reverse=TRUE))+
      scale_x_continuous(labels=levels(data$trophic_group), breaks=c(1:length(levels(data$trophic_group))))+
      coord_flip()+
      scale_y_log10(breaks = 10^seq(-1,3,1),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      ylab(label_biomass)

# relative biomass distribution
data$relative<-0
data$relative[data$trophic_group_file==label_multichannel]<-data$biomass[data$trophic_group_file==label_multichannel]/sum(data$biomass[data$trophic_group_file==label_multichannel])
data$relative[data$trophic_group_file==label_size_structured]<-data$biomass[data$trophic_group_file==label_size_structured]/sum(data$biomass[data$trophic_group_file==label_size_structured])
data$relative<-data$relative*100
p2<-ggplot(data=data)+
      geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_wrap(~trophic_group_file)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  axis.line.y=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      scale_y_continuous(labels = scales::percent)+
      ylab(label_biomass_fraction)

data_model<-data # stores the modelling results for comparison with empirical data

# Biomass distribution empiric # ----
# Johnston and Sibly 2018 Nature Ecology and Evolution
data_abundance<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaAbundanceData.csv",header=T,sep=",")
data_bodymass<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv",header=T,sep=",")
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
Mkm2_to_m2=1e6*1e6 # million km² to km²
Pg_to_mg=1e18 # petagram to milligram

data_bodymass$Individual_bodymass_mg_FM<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry
data_bodymass$trophic_group<-"macro-food web"
data_bodymass$trophic_group[data_bodymass$Individual_bodymass_mg_FM<0.001]="micro-food web"

data_abundance<-merge(data_abundance,data_bodymass,by=c("Soil_biota_group"))
data_abundance$Individual_bodymass_mg_DM<-data_abundance$Individual_bodymass_mg_DM*dry_to_C # conversion into mgC
data_abundance$biomass_density<-data_abundance$Biomass_g_m2*1e3*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance<-data_abundance[,c("Soil_biota_group","Biome_ecosystem","Biomass_g_m2","Individual_bodymass_mg_DM","trophic_group","biomass_density")]

# bacteria<-data.frame(Soil_biota_group="Bacteria",
#                      Biome_ecosystem=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),
#                      Biomass_g_m2=c(22.8,23.5,22.6,22.4,22.3), # ln(Individuals) actually
#                      Individual_bodymass_mg_DM=2.83e-8*dry_to_C, # converted into C body mass
#                      trophic_group="microbes")
# bacteria$biomass_density<-exp(bacteria$Biomass_g_m2)*bacteria$Individual_bodymass_mg_DM

microbes<-read.table("Data/Metabolic_rates_Johnston_2018/Xu_2013_Global_Ecology_and_Biogeography_aggregated_data.csv",header=T,sep=",") # data from the main text of Xu et al 2013 Global Ecology and Biogeography
microbes$biomass_density<-microbes$soil_microbial_C_Pg_0_30_cm*Pg_to_mg/(microbes$Area_Mkm2*Mkm2_to_m2)
microbes$Soil_biota_group="microbes"
microbes$trophic_group="microbes"

data_abundance<-rbind(data_abundance[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")],
                      microbes[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")])
data_abundance$trophic_group<-factor(data_abundance$trophic_group,levels=c("microbes","micro-food web","macro-food web"))
data_abundance<-data_abundance[data_abundance$Biome_ecosystem%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data_abundance$Biome_ecosystem<-factor(data_abundance$Biome_ecosystem,levels=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"))
levels(data_abundance$Biome_ecosystem)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")

data<-aggregate(biomass_density~trophic_group+Soil_biota_group,data_abundance,mean)
data<-aggregate(biomass_density~trophic_group,data,sum)
data$facet="Empirical data"

# p4<-ggplot(data=data)+
#   geom_bar(aes(1,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
#   theme+theme(legend.title=element_blank(),
#               legend.justification = "left")+
#   scale_fill_manual(values=c("lightgrey","darkgrey","black"),
#                     guide = guide_legend(reverse=TRUE))
# legend_1<-get_legend(p4)

p4<-ggplot(data=data)+
      geom_bar(aes(1,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_wrap(~facet)+
      theme+theme(#legend.position="none",
                  legend.title=element_blank(),
                  axis.title=element_blank(),
                  axis.text=element_blank(),
                  axis.ticks=element_blank(),
                  axis.line.y=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                        guide = guide_legend(reverse=TRUE))+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

p4_zoom<-ggplot(data=data)+
            geom_bar(aes(1,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
            theme+theme(legend.position="none",
              legend.title=element_blank(),
              axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x=element_blank(),
              axis.text = element_text(size=15),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
            scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                              guide = guide_legend(reverse=TRUE))+
            scale_y_continuous(breaks=seq(0.945,0.95,0.0025),labels = scales::percent,
                               limits=c(0.945,0.95),oob = rescale_none,position = "right")+
            ylab(label_biomass_fraction)

data$persent<-data$biomass_density/sum(data$biomass_density)*100

# Heděnec et al 2022 Scientific Reports
data_hedenec<-read.table(paste(path_hedenec,"biomass_invertebrates_raw.csv",sep=""),sep=",",header=T)
data_hedenec_trophic_groups<-read.table(paste(path_hedenec,"invertebrates_trophic_groups.csv",sep=""),sep=",",header=T)
data_hedenec<-merge(data_hedenec,data_hedenec_trophic_groups,by=c("Guilds","Fauna"))
rm(data_hedenec_trophic_groups)

data<-aggregate(Biomass~trophic_type+size_class, data=data_hedenec, sum)
data$size_class<-factor(data$size_class,levels=c("micro-foodweb","macro-foodweb"))
levels(data$size_class)=c("micro-food web","macro-food web")
data$trophic_type<-factor(data$trophic_type,levels=c("carnivores","herbivores","detritivores","omnivores","microbivores"))
data<-data[data$size_class=="macro-food web",]
data$Biomass[data$trophic_type=="carnivores"]=data$Biomass[data$trophic_type=="carnivores"]+data$Biomass[data$trophic_type=="omnivores"]
data<-data[-which(data$trophic_type%in%c("herbivores","omnivores")),]
data$facet="Empirical data"

data_model<-data_model[,c("trophic_group","trophic_group_file","biomass")]
names(data_model)=c("trophic_type","facet","Biomass")
data_model<-data_model[data_model$facet=="Multichannel",]
data_model$size_class="macro-food web"
data_model$Biomass[4]<-data_model$Biomass[4]+data_model$Biomass[6]
data_model<-data_model[c(4,5),c(1,4,3,2)]
data_model$trophic_type=c("detritivores","carnivores")
data<-rbind(data,data_model)

p3<-ggplot(data=data[data$size_class=="macro-food web",])+
  geom_bar(aes(size_class,Biomass,fill=trophic_type),stat="identity",position=position_fill())+
  theme+theme(legend.title=element_blank())+
  scale_fill_manual(values=c("red2","slateblue1","seagreen2"))
legend<-get_legend(p3)

p3<-ggplot(data=data[data$size_class=="macro-food web",])+
      geom_bar(aes(size_class,Biomass,fill=trophic_type),stat="identity",position=position_fill())+
      facet_wrap(~facet)+
      theme+theme(legend.position="none",
                  legend.title=element_blank(),
                  axis.ticks.x=element_blank(),
                  #axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("red2","slateblue1","seagreen2"))+
      scale_x_discrete(expand = c(0,0.5))+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      #xlab("Food web")+
      xlab("Macro-food web")+
      ylab(label_biomass_fraction)

### Final graph #----
temp<-plot_grid(p3, NULL, p2, p4, labels = NULL, align = "h", rel_widths = c(0.7, 0.1, 0.55, 0.65),nrow=1)

graph<-ggdraw(xlim = c(0,2), ylim = c(0,2)) +
  draw_plot(legend, 0.05, 1.05, 0.3, 0.3)+
  draw_plot(p1, 0.4, 1, 1.5, 1)+
  draw_plot(temp, 0, 0, 2, 1)+
  #draw_plot(p3, 0, 0, 0.7, 1)+
  #draw_plot(p2, 0.8, 0, 0.55, 1)+
  #draw_plot(p4, 1.33, 0, 0.65, 1)+
  draw_plot(p4_zoom, 1.71, 0.65, 0.25, 0.3)+
  draw_plot_label(c("A","B","C"), c(0.35,0,0.75), c(2,1.03,1.03), size = 25)
ggsave(paste(path_figures,"figure_biomass.pdf",sep=""), graph, width = 12, height = 8, device=cairo_pdf)

# ROBUSTNESS OVER ENVIRONMENTAL GRADIENTS # ----
# Predicted biomass distribution in each biome # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_multi_channel<-Rload_data(path_results,1)
data_size_structured<-Rload_data(path_results,3)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)
data_species<-rbind(data_multi_channel$data_species,
                    data_size_structured$data_species)
data_species$trophic_group_file<-as.factor(data_species$trophic_group_file)
levels(data_species$trophic_group_file)=c(label_multichannel,label_size_structured)
data<-aggregate(biomass~trophic_group+trophic_group_file+input_FOM+input_DOC, data=data_species, sum, na.rm = TRUE)
#data$biomass[data$biomass<1e-7]=NA
data$foodweb<-"micro-food web"
data$foodweb[data$trophic_group%in%c("macro-food web detritivores","macro-food web carnivores","trophic whales")]<-"macro-food web"
input_DOC<-unique(data$input_DOC)
input_FOM<-unique(data$input_FOM)
data_FOM_DOC<-read.table(paste0(path_data,"data_FOM_SOM.csv"),header=T,sep=",") # data of FOM and DOC inputs in the main land ecosystems
for (i in 1:nrow(data_FOM_DOC)){
  data_FOM_DOC$input_DOC[i]<-input_DOC[which.min(abs(input_DOC-data_FOM_DOC$input_DOC[i]))] # place the values by the nearest values used for the simulations
  data_FOM_DOC$input_FOM[i]<-input_FOM[which.min(abs(input_FOM-data_FOM_DOC$input_FOM[i]))]
}
data<-merge(data,data_FOM_DOC,by=c("input_DOC","input_FOM"))
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
data$Biome<-factor(data$Biome,levels=c("Desert","Tundra","Boreal forest","Temperate forest","Mediterranean vegetation","Temperate grassland","Tropical grassland","Tropical forest"))
data<-data[data$Biome%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data$Biome<-droplevels(data$Biome)
levels(data$Biome)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")

p1<-ggplot(data=data)+
  geom_bar(aes(Biome,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+theme(legend.margin=margin(l = 0, unit='pt'))+
  palette_multi_channel$scale_colour_bar_trophic_groups
legend_1<-get_legend(p1)

p1<-ggplot(data=data)+
      geom_bar(aes(Biome,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_wrap(~trophic_group_file)+
      theme+theme(legend.position="none",
                  axis.text.y=element_text(size=15),
                  axis.title.y=element_blank(),
                  axis.title.x=element_text(hjust=1),
                  panel.grid=element_blank(),
                  panel.spacing.x = unit(1.5, "lines"))+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      coord_flip()+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

# Biomass empiric #----
data_abundance<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaAbundanceData.csv",header=T,sep=",")
data_bodymass<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv",header=T,sep=",")
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
Mkm2_to_m2=1e6*1e6 # million km² to km²
Pg_to_mg=1e18 # petagram to milligram

data_bodymass$Individual_bodymass_mg_FM<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry
data_bodymass$trophic_group<-"macro-food web"
data_bodymass$trophic_group[data_bodymass$Individual_bodymass_mg_FM<0.001]="micro-food web"

data_abundance<-merge(data_abundance,data_bodymass,by=c("Soil_biota_group"))
data_abundance$Individual_bodymass_mg_DM<-data_abundance$Individual_bodymass_mg_DM*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance$biomass_density<-data_abundance$Biomass_g_m2*1e3*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance<-data_abundance[,c("Soil_biota_group","Biome_ecosystem","Biomass_g_m2","Individual_bodymass_mg_DM","trophic_group","biomass_density")]

microbes<-read.table("Data/Metabolic_rates_Johnston_2018/Xu_2013_Global_Ecology_and_Biogeography_aggregated_data.csv",header=T,sep=",") # data from the main text of Xu et al 2013 Global Ecology and Biogeography
microbes$biomass_density<-microbes$soil_microbial_C_Pg_0_30_cm*Pg_to_mg/(microbes$Area_Mkm2*Mkm2_to_m2)
microbes$Soil_biota_group="microbes"
microbes$trophic_group="microbes"

data_abundance<-rbind(data_abundance[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")],
                      microbes[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")])
data_abundance$trophic_group<-factor(data_abundance$trophic_group,levels=c("microbes","micro-food web","macro-food web"))
data_abundance<-data_abundance[data_abundance$Biome_ecosystem%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data_abundance$Biome_ecosystem<-factor(data_abundance$Biome_ecosystem,levels=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"))
levels(data_abundance$Biome_ecosystem)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")

data<-aggregate(biomass_density~Biome_ecosystem+trophic_group+Soil_biota_group,data_abundance,mean)
data<-aggregate(biomass_density~Biome_ecosystem+trophic_group,data,sum)
data$facet="Empirical data"

p2<-ggplot(data=data)+
  geom_bar(aes(Biome_ecosystem,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+theme(legend.title=element_blank(),
              legend.margin=margin(l = 0, unit='pt'))+
  scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                    guide = guide_legend(reverse=TRUE))
legend_2<-get_legend(p2)

p2<-ggplot(data=data)+
      geom_bar(aes(Biome_ecosystem,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_wrap(~facet)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.y=element_blank(),
                  axis.line.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  panel.grid=element_blank())+
      scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                        guide = guide_legend(reverse=TRUE))+
      coord_flip()+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

### Final graph #----
temp<-plot_grid(p1,NULL, p2, labels = NULL, align = "h", rel_widths = c(2,0.05, 0.8),nrow=1)

graph<-ggdraw(xlim = c(0,3.6), ylim = c(0,1)) +
  #draw_plot(p1, 0, 0, 2, 1)+
  #draw_plot(p2, 2, 0, 0.9, 1)+
  draw_plot(temp, 0, 0, 2.9, 1)+
  draw_plot(legend_1, 2.9, 0.5, 0.75, 0.3)+
  draw_plot(legend_2, 2.9, 0.1, 0.7, 0.3)
ggsave(paste(path_figures,"figure_biomass_biome.pdf",sep=""), graph, width = 12, height = 4, device=cairo_pdf)

########################## ----
# SUPPORTING INFORMATION # ----
########################## ----
# load data # ----
path_results<-"results_final_main/"#"results_23_02_15/"#"results_23_02_03/"#"results_23_01_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

data_size_structured<-Rload_data(path_results,3)
palette_size_structured<-Rmake_palettes(data_size_structured$species_traits,
                                        data_size_structured$n_tot_species,
                                        data_size_structured$n_faeces)

data_detritivore_spectrum<-Rload_data(path_results,2)
palette_detritivore_spectrum<-Rmake_palettes(data_detritivore_spectrum$species_traits,
                                             data_detritivore_spectrum$n_tot_species,
                                             data_detritivore_spectrum$n_faeces)

params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
# species data
data_species<-rbind(data_multi_channel$data_species,
                    data_size_structured$data_species)
data_species$trophic_group_file<-as.factor(data_species$trophic_group_file)
levels(data_species$trophic_group_file)=c(label_multichannel,label_size_structured)
# detritus data
data_detritus<-data_multi_channel$data_detritus
# detritus production data
production<-data_multi_channel$data_detritus_production

# Food web diagrams # ----
get_data_foodweb<-function(input_FOM,input_DOC,data_results,palette){
  with(data_results,{
    data<-data_species[data_species$input_FOM==input_FOM & data_species$input_DOC==input_DOC,]
    # changes the position of macro-food web carnivores to have them above detritivores
    if (length(data$trophic_group=="macro-food web detritivores")>0){
      data$num[data$trophic_group=="macro-food web carnivores"]=data$num[data$trophic_group=="macro-food web carnivores"]-length(data$trophic_group[data$trophic_group=="macro-food web detritivores"])
    }
    diet<-read.table(paste(path_results,"diet_",data$simu_ID[1],".txt",sep=""),sep=",",header=T)
    # remove irrelevant species (almoast extinct)
    diet[data$num[data$TL==0],]=NA
    diet[,data$num[data$TL==0]]=NA
    diet$pred<-names_species
    diet<-Rmake_table_for_plot(diet,"pred","prey","flow") # makes diet a long table
    levels(diet$prey)<-names_tot
    # add the information on prey
    data_prey<-data[,c("bodymass","TL","names","num")]
    names(data_prey)<-c("bodymass_prey","TL_prey","prey","num_prey")
    diet<-merge(diet,data_prey,by=c("prey"),all=T) # add the information on prey to the data set
    diet$TL_prey[is.na(diet$TL_prey)]=0 # detritus
    diet$bodymass_prey[is.na(diet$bodymass_prey)]=0 # body mass
    names_abio<-c("DOC","SOM",names_faeces,"FOM","N")
    for (i in 1:(n_tot_detritus+n_nutrients)){
      diet$num_prey[diet$prey==names_abio[i]]=i # sets the position of detritus
    }
    plus_SOM=6
    plus_faeces=10
    plus_FOM=18
    plus_N=26
    diet$num_prey[diet$prey=="SOM"]=diet$num_prey[diet$prey=="SOM"]+plus_SOM # to make the dots more readable
    diet$num_prey[which(diet$prey%in%names_faeces)]=diet$num_prey[which(diet$prey%in%names_faeces)]+plus_faeces
    diet$num_prey[diet$prey=="FOM"]=diet$num_prey[diet$prey=="FOM"]+plus_FOM
    # add the information on predators
    data_pred<-data[,which(names(data)%in%c("bodymass","TL","names","num"))]
    names(data_pred)<-c("pred","TL_pred","bodymass_pred","num_pred")
    diet<-merge(diet,data_pred,by=c("pred"),all=T) # add the information on predators to the data set
    diet<-diet[diet$prey!=diet$pred,] # removes self-interactions
    diet<-diet[diet$flow!=0 & is.na(diet$flow)==F,] # removes empty flows
    
    # position of nodes
    nodes<-data[,c("names","TL","num")]
    #nodes<-nodes[order(nodes$names),]
    #nodes$num<-c(1:n_tot_species)
    nodes$TL[nodes$TL==0]=NA
    temp<-data.frame(names=names_abio)
    temp$TL=0
    temp$num<-c(1:(n_tot_detritus+n_nutrients))
    nodes<-rbind(nodes,temp)
    rm(temp)
    nodes$num[nodes$names=="SOM"]=nodes$num[nodes$names=="SOM"]+plus_SOM # to make the dots more readable
    nodes$num[which(nodes$names%in%names_faeces)]=nodes$num[which(nodes$names%in%names_faeces)]+plus_faeces
    nodes$num[nodes$names=="FOM"]=nodes$num[nodes$names=="FOM"]+plus_FOM
    nodes$num[nodes$names=="N"]=nodes$num[nodes$names=="N"]+plus_N
    nodes$names<-factor(nodes$names,levels=names_tot)
    # removes species without interactions
    for (i in 1:n_tot_species){
      if (length(diet$flow[diet$pred==nodes$names[i]])==0){
        nodes$TL[i]=NA
        nodes$num[i]=NA
      }
    }
    # removes empty compartments
    for (i in (n_tot_species+1):n_tot){
      if (data_density$biomass[data_density$input_FOM==input_FOM &
                               data_density$input_DOC==input_DOC &
                               data_density$names==nodes$names[i]]==0){
        nodes$TL[i]=NA
        nodes$num[i]=NA
        nodes$num[(i+1):n_tot]=nodes$num[(i+1):n_tot]-1
        diet$num_prey[diet$prey%in%nodes$names[(i+1):n_tot]]=diet$num_prey[diet$prey%in%nodes$names[(i+1):n_tot]]-1
      }
    }
    return(list(diet=diet,
                nodes=nodes,
                scale_colour_line_total=palette$scale_colour_line_total))
  })
}

# multi-channel
data_foodweb<-get_data_foodweb(300,150,data_multi_channel,palette_multi_channel)
data_foodweb$nodes[data_foodweb$nodes$names=="N",c("TL","num")]=NA
p1<-ggplot()+
  geom_curve(data=data_foodweb$diet,
             aes(x=num_prey, xend=num_pred,
                 y=TL_prey, yend=TL_pred,
                 size=flow, alpha=flow),
             curvature=0)+
  geom_point(data=data_foodweb$nodes,
             aes(num,TL,colour=names),size=7)+
  data_foodweb$scale_colour_line_total+
  scale_size_continuous(name=expression("Flow (mgC m"^-2*" day"^-1*")"))+
  scale_alpha(name=expression("Flow (mgC m"^-2*" day"^-1*")"))+
  guides(colour="none")+
  theme+theme(legend.position=c(0.01, 1),
              legend.justification=c(0, 1),
              legend.background=element_blank(),
              text=element_text(size=30),
              axis.text=element_text(size=30),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x = element_line(arrow = arrow(angle = 15, length = unit(.15,"inches"))),
              panel.grid=element_blank())+
  scale_y_continuous(breaks=c(0:floor(max(data_foodweb$nodes$TL,na.rm=TRUE))),
                     limits=c(-0.5,NA))+
  xlab(label_bodymass_approx)+
  ylab(label_TL)+
  ggtitle(label_multichannel)
#p1<-align_legend(p1,1)

# size structured
data_foodweb<-get_data_foodweb(300,150,data_size_structured,palette_size_structured)
data_foodweb$nodes[data_foodweb$nodes$names=="N",c("TL","num")]=NA
p2<-ggplot()+
  geom_curve(data=data_foodweb$diet,
             aes(x=num_prey, xend=num_pred,
                 y=TL_prey, yend=TL_pred,
                 size=flow, alpha=flow),
             curvature=0)+
  geom_point(data=data_foodweb$nodes,
             aes(num,TL,colour=names),size=7)+
  data_foodweb$scale_colour_line_total+
  scale_size_continuous(name=expression("Flow (mgC m"^-2*" day"^-1*")"))+
  scale_alpha(name=expression("Flow (mgC m"^-2*" day"^-1*")"))+
  guides(colour="none")+
  theme+theme(legend.position=c(0.01, 1),
              legend.justification=c(0, 1),
              legend.background=element_blank(),
              text=element_text(size=30),
              axis.text=element_text(size=30),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              axis.line.x = element_line(arrow = arrow(angle = 15, length = unit(.15,"inches"))),
              panel.grid=element_blank())+
  scale_y_continuous(breaks=c(0:floor(max(data_foodweb$nodes$TL,na.rm=TRUE))),
                     limits=c(-0.5,NA))+
  xlab(label_bodymass_approx)+
  ylab(label_TL)+
  ggtitle(label_size_structured)

# legend
data<-data_species[data_species$trophic_group_file==label_multichannel,c("trophic_group","names")]
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
databis<-data_detritus[data_detritus$trophic_group_file==1,c("type","names")]
names(databis)<-c("trophic_group","names")
data<-rbind(data,databis)
rm(databis)
data$trophic_group<-factor(data$trophic_group,c("DOC","SOM","faeces","FOM","microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales"))

legend<-ggplot(data=data)+
  geom_bar(aes(trophic_group,fill=names),stat="count")+
  theme_void()+theme(legend.position="none",
                     axis.text.y=element_text(family="serif",size=15,hjust=1))+
  palette_multi_channel$scale_colour_bar_total+
  coord_flip()+
  ylab(label_biomass_fraction)

# final graph
space=0.1
graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 2+space)) +
  draw_plot(p1, 0, 1+space, 1, 1)+
  draw_plot(p2, 0, 0, 1, 1)+
  draw_plot(legend, 0.65, 1.5+space, 0.3, 0.5)+
  draw_plot_label(c("A","B"), c(0,0), c(2+space,1), size = 30)
ggsave(paste(path_figures,"supp_foodweb.pdf",sep=""), graph, width = 16, height = 14, device=cairo_pdf)

# Biomass distribution per species # ----
p1<-ggplot(data_multi_channel$data_species)+
  geom_rect(aes(xmin=as.numeric(names)-0.5,xmax=as.numeric(names)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=names))+
  theme+theme(legend.position="none",
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_blank(),
              plot.margin = margin(t=5.5,r=5.5,b=5.5,l=40, "pt"))+ # 5.5 is the default value
  palette_multi_channel$scale_colour_bar_species+
  coord_flip()+
  scale_y_log10(breaks = 10^c(0:3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  ylab(label_biomass)+
  ggtitle(label_multichannel)

p2<-ggplot(data_detritivore_spectrum$data_species)+
  geom_rect(aes(xmin=as.numeric(names)-0.5,xmax=as.numeric(names)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=names))+
  annotate("errorbarh",xmin=24,xmax=32,y=1,height = 0.2)+
  annotate("text",x=29,y=10^1.2,label="extinct detritivores",family="serif",size=6)+
  theme+theme(legend.position="none",
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_blank(),
              plot.margin = margin(t=5.5,r=5.5,b=5.5,l=40, "pt"))+ # 5.5 is the default value
  palette_detritivore_spectrum$scale_colour_bar_species+
  coord_flip()+
  scale_y_log10(breaks = 10^c(0:3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  ylab(label_biomass)+
  ggtitle(label_detritivores_spectrum)

p3<-ggplot(data_size_structured$data_species)+
  geom_rect(aes(xmin=as.numeric(names)-0.5,xmax=as.numeric(names)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=names))+
  theme+theme(legend.position="none",
              axis.title.y=element_blank(),
              axis.text.y=element_blank(),
              axis.ticks.y=element_blank(),
              panel.grid.minor.y=element_blank(),
              panel.grid.major.y=element_blank(),
              plot.margin = margin(t=5.5,r=5.5,b=5.5,l=40, "pt"))+ # 5.5 is the default value
  palette_size_structured$scale_colour_bar_species+
  coord_flip()+
  scale_y_log10(breaks = 10^c(0:3),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  ylab(label_biomass)+
  ggtitle(label_size_structured)

# Biomass distribution compared with Potapov
data_potapov<-read.table(paste(path_potapov,"RawData_MassDensityTrophicLevel.csv",sep=""),sep=",",header=T)
#data_potapov<-data_potapov[,c("SizeClass","Biomass_mg_gC")]
data_potapov<-data_potapov %>% separate(SizeClass, c("bodymass",NA), sep = "-")
data_potapov$bodymass[data_potapov$bodymass==">1000"]=1000
data_potapov$bodymass<-as.numeric(data_potapov$bodymass)
data_potapov<-aggregate(Biomass_mg_gC~bodymass+Layer+TrophicGroup,data=data_potapov,mean) # mean over sites
data_potapov<-aggregate(Biomass_mg_gC~bodymass,data=data_potapov,sum) # sum of trophic groups inside each size class

p4<-ggplot(data_potapov)+
  geom_ribbon(aes(x=bodymass,ymin=min(Biomass_mg_gC), ymax=Biomass_mg_gC),fill="lightgray")+
  theme+theme(legend.position="none")+ 
  palette_multi_channel$scale_colour_line_species+
  coord_flip()+
  x_axis_log10+
  y_axis_log10+
  xlab(label_bodymass)+
  ylab(label_biomass_potapov)+
  ggtitle("Empirical")

graph<-ggdraw(xlim = c(0, 4), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)+
  draw_plot(p4, 3, 0, 1, 1)+
  draw_plot_label(c("A","B","C","D"), c(0,1,2,3), c(1,1,1,1), size = 25)
ggsave(paste(path_figures,"supp_biomass_species.pdf",sep=""), graph, width = 15, height = 5, device=cairo_pdf)

# Biomass multichannel FOM - DOC # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

data<-data_multi_channel$data_species
data<-aggregate(biomass~trophic_group+input_FOM+input_DOC, data=data, sum, na.rm = TRUE)
#data$biomass[data$biomass<1e-7]=NA
data$foodweb<-"micro-food web"
data$foodweb[data$trophic_group%in%c("macro-food web detritivores","macro-food web carnivores","trophic whales")]<-"macro-food web"
input_DOC<-unique(data$input_DOC)
input_FOM<-unique(data$input_FOM)
data_point<-expand_grid(input_DOC=c(input_DOC[1],input_DOC[length(input_DOC)]),
                        input_FOM=c(input_FOM[1],input_FOM[length(input_FOM)]))

levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
p1<-ggplot(data=data[data$input_FOM==data_point$input_FOM[1] & data$input_DOC==data_point$input_DOC[1],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+theme(legend.justification = "left")+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  guides(fill = guide_legend(reverse=TRUE))
legend_bar<-get_legend(p1)

p1<-ggplot(data=data[data$input_FOM==data_point$input_FOM[1] & data$input_DOC==data_point$input_DOC[1],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  facet_grid(input_FOM~input_DOC)+
  theme+theme(legend.position="none",
              axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
  ylab(label_biomass_fraction)

p2<-ggplot(data=data[data$input_FOM==data_point$input_FOM[2] & data$input_DOC==data_point$input_DOC[2],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  facet_grid(input_FOM~input_DOC)+
  theme+theme(legend.position="none",
              axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
  ylab(label_biomass_fraction)

# temp<-data[data$input_FOM==data_point$input_FOM[3] & data$input_DOC==data_point$input_DOC[3],]
# temp$biomass<-temp$biomass/sum(temp$biomass)*100
# temp
p3<-ggplot(data=data[data$input_FOM==data_point$input_FOM[3] & data$input_DOC==data_point$input_DOC[3],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  facet_grid(input_FOM~input_DOC)+
  theme+theme(legend.position="none",
              axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
  ylab(label_biomass_fraction)

p4<-ggplot(data=data[data$input_FOM==data_point$input_FOM[4] & data$input_DOC==data_point$input_DOC[4],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  facet_grid(input_FOM~input_DOC)+
  theme+theme(legend.position="none",
              axis.title=element_blank(),
              axis.text.x=element_blank(),
              axis.ticks.x=element_blank(),
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank(),
              plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
  palette_multi_channel$scale_colour_bar_trophic_groups+
  scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
  ylab(label_biomass_fraction)

data$input_DOC[data$input_DOC==input_DOC[1]]=0 # to have a clean raster
data<-aggregate(biomass~foodweb+input_FOM+input_DOC, data=data, sum)
data<-dcast(data,input_FOM+input_DOC~foodweb,value.var="biomass")
data$ratio<-data$'macro-food web'/(data$'macro-food web'+data$'micro-food web')
data_FOM_DOC<-read.table(paste0(path_data,"data_FOM_SOM.csv"),header=T,sep=",") # data of FOM and DOC inputs in the main land ecosystems
data_FOM_DOC$type="empiric"
data_FOM_DOC<-rbind(data_FOM_DOC,list("Figure 2",300,150,"model"))

p5<-ggplot(data=data)+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  theme_raster+theme(legend.justification = "left")+
  scale_fill_viridis(option="viridis",
                     name="Macro-food web\nbiomass fraction",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)
legend_raster<-get_legend(p5)

p5<-ggplot(data=data)+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  geom_point(data=data_FOM_DOC,aes(input_DOC,input_FOM,colour=type),size=5,shape=18)+
  coord_cartesian(clip = "off")+ # allow the labels to go beyond the edges of the panel
  geom_label_repel(data=data_FOM_DOC,aes(input_DOC,input_FOM,label=Biome,colour=type),
                   point.padding=1,size=4,
                   fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+ # allow the labels to go beyond the edges of the panel
  #annotate(geom="point",x=150,y=300,colour="red",size=5,shape="O")+ # point from the main figure
  geom_point(data=data_point,aes(input_DOC,input_FOM),colour="deepskyblue2",size=5)+
  theme_raster+theme(legend.position="none")+
  scale_fill_viridis(option="viridis",
                     limits=c(0,1))+
  scale_colour_data_type+
  xlab(label_DOC_input)+
  ylab(label_FOM_input)+
  ggtitle(label_multichannel)

label_margin=0.15
top_margin=0.1
vspace=0.05
width_fraction=0.6
width_legend=1
width_tot=width_fraction*2+2+width_legend
height_tot=2+top_margin+vspace
width_fraction=0.6
width_legend=1
width_tot=width_fraction*2+2+width_legend
graph<-ggdraw(xlim = c(0, width_tot), ylim = c(0, height_tot)) +
  draw_plot(p1, label_margin, vspace, width_fraction, 1)+
  draw_plot(p2, label_margin, 1+vspace*2, width_fraction, 1)+
  draw_plot(p5, label_margin+width_fraction, vspace, 2, 2)+
  draw_plot(p3, label_margin+width_fraction+2, vspace, width_fraction, 1)+
  draw_plot(p4, label_margin+width_fraction+2, 1+vspace*2, width_fraction, 1)+
  draw_plot(legend_bar, label_margin+2*width_fraction+2.1, 0, width_legend, 1)+
  draw_plot(legend_raster, label_margin+2*width_fraction+2.1, 1, width_legend, 1)
ggsave(paste(path_figures,"supp_DOC_FOM_multichannel.pdf",sep=""), graph, width = 12, height = 7, device=cairo_pdf)

# Biomass size-structured FOM - DOC # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_size_structured<-Rload_data(path_results,3)
palette_size_structured<-Rmake_palettes(data_size_structured$species_traits,
                                        data_size_structured$n_tot_species,
                                        data_size_structured$n_faeces)

data<-data_size_structured$data_species
data<-aggregate(biomass~trophic_group+input_FOM+input_DOC, data=data, sum, na.rm = TRUE)
#data$biomass[data$biomass<1e-7]=NA
data$foodweb<-"micro-food web"
data$foodweb[data$trophic_group%in%c("macro-food web detritivores","macro-food web carnivores","trophic whales")]<-"macro-food web"
input_DOC<-unique(data$input_DOC)
input_FOM<-unique(data$input_FOM)
data_point<-expand_grid(input_DOC=c(input_DOC[1],input_DOC[length(input_DOC)]),
                        input_FOM=c(input_FOM[1],input_FOM[length(input_FOM)]))

levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-carnivores")
p1<-ggplot(data=data[data$input_FOM==data_point$input_FOM[1] & data$input_DOC==data_point$input_DOC[1],])+
  geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+theme(legend.justification = "left")+
  palette_size_structured$scale_colour_bar_trophic_groups+
  guides(fill = guide_legend(reverse=TRUE))
legend_bar<-get_legend(p1)

p1<-ggplot(data=data[data$input_FOM==data_point$input_FOM[1] & data$input_DOC==data_point$input_DOC[1],])+
      geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_grid(input_FOM~input_DOC)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank(),
                  plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
      palette_size_structured$scale_colour_bar_trophic_groups+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

p2<-ggplot(data=data[data$input_FOM==data_point$input_FOM[2] & data$input_DOC==data_point$input_DOC[2],])+
      geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_grid(input_FOM~input_DOC)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank(),
                  plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
      palette_size_structured$scale_colour_bar_trophic_groups+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

# temp<-data[data$input_FOM==data_point$input_FOM[3] & data$input_DOC==data_point$input_DOC[3],]
# temp$biomass<-temp$biomass/sum(temp$biomass)*100
# temp
p3<-ggplot(data=data[data$input_FOM==data_point$input_FOM[3] & data$input_DOC==data_point$input_DOC[3],])+
      geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_grid(input_FOM~input_DOC)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank(),
                  plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
      palette_size_structured$scale_colour_bar_trophic_groups+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

p4<-ggplot(data=data[data$input_FOM==data_point$input_FOM[4] & data$input_DOC==data_point$input_DOC[4],])+
      geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_grid(input_FOM~input_DOC)+
      theme+theme(legend.position="none",
                  axis.title=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank(),
                  plot.background=element_rect(colour = "black", fill=NA, linetype="solid"))+
      palette_size_structured$scale_colour_bar_trophic_groups+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)

data$input_DOC[data$input_DOC==input_DOC[1]]=0 # to have a clean raster
data<-aggregate(biomass~foodweb+input_FOM+input_DOC, data=data, sum)
data<-dcast(data,input_FOM+input_DOC~foodweb,value.var="biomass")
data$ratio<-data$'macro-food web'/(data$'macro-food web'+data$'micro-food web')
data_FOM_DOC<-read.table(paste0(path_data,"data_FOM_SOM.csv"),header=T,sep=",") # data of FOM and DOC inputs in the main land ecosystems
data_FOM_DOC$type="empiric"
data_FOM_DOC<-rbind(data_FOM_DOC,list("Figure 2",300,150,"model"))

p5<-ggplot(data=data)+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  theme_raster+theme(legend.justification = "left")+
  scale_fill_viridis(option="viridis",
                     name="Macro-food web\nbiomass fraction",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)
legend_raster<-get_legend(p5)

p5<-ggplot(data=data)+
      geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
      geom_point(data=data_FOM_DOC,aes(input_DOC,input_FOM,colour=type),size=5,shape=18)+
      coord_cartesian(clip = "off")+ # allow the labels to go beyond the edges of the panel
      geom_label_repel(data=data_FOM_DOC,aes(input_DOC,input_FOM,label=Biome,colour=type),
                       point.padding=1,size=4,
                       fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+ # allow the labels to go beyond the edges of the panel
      #annotate(geom="point",x=150,y=300,colour="red",size=5,shape="O")+ # point from the main figure
      geom_point(data=data_point,aes(input_DOC,input_FOM),colour="deepskyblue2",size=5)+
      theme_raster+theme(legend.position="none")+
      scale_fill_viridis(option="viridis",
                         limits=c(0,1))+
      scale_colour_data_type+
      xlab(label_DOC_input)+
      ylab(label_FOM_input)+
      ggtitle(label_size_structured)

label_margin=0.15
top_margin=0.1
vspace=0.05
width_fraction=0.6
width_legend=1
width_tot=width_fraction*2+2+width_legend
height_tot=2+top_margin+vspace
width_fraction=0.6
width_legend=1
width_tot=width_fraction*2+2+width_legend
graph<-ggdraw(xlim = c(0, width_tot), ylim = c(0, height_tot)) +
  draw_plot(p1, label_margin, vspace, width_fraction, 1)+
  draw_plot(p2, label_margin, 1+vspace*2, width_fraction, 1)+
  draw_plot(p5, label_margin+width_fraction, vspace, 2, 2)+
  draw_plot(p3, label_margin+width_fraction+2, vspace, width_fraction, 1)+
  draw_plot(p4, label_margin+width_fraction+2, 1+vspace*2, width_fraction, 1)+
  draw_plot(legend_bar, label_margin+2*width_fraction+2.1, 0, width_legend, 1)+
  draw_plot(legend_raster, label_margin+2*width_fraction+2.1, 1, width_legend, 1)
ggsave(paste(path_figures,"supp_DOC_FOM_sizestructured.pdf",sep=""), graph, width = 12, height = 7, device=cairo_pdf)

# Trophic levels detritivore spectrum # ---- 
data<-data_detritivore_spectrum$data_species
data_detritivore_spectrum$data_species$trophic_group<-as.factor(data_detritivore_spectrum$data_species$trophic_group)
levels(data_detritivore_spectrum$data_species$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
data$TL[data$TL==0]=NA

data_potapov<-read.table(paste(path_potapov,"RawData_MassDensityTrophicLevel.csv",sep=""),sep=",",header=T)
data_potapov<-data_potapov[,c("SizeClass","TL")]
data_potapov$TL<-data_potapov$TL-1
data_potapov<-data_potapov %>% separate(SizeClass, c("bodymass",NA), sep = "-")
data_potapov$bodymass[data_potapov$bodymass==">1000"]=1000
data_potapov$bodymass<-as.numeric(data_potapov$bodymass)

p1<-ggplot(data=data)+
      geom_smooth(data=data_potapov, aes(bodymass,TL), se=FALSE, colour="black", linetype="dashed")+
      geom_smooth(aes(bodymass,TL), colour="black", linetype="solid")+
      geom_point(aes(bodymass,TL,colour=names),size=4)+
      annotate("errorbarh",y=1,xmin=1e-7,xmax=1e-4,height = 0.2)+
      annotate("text",x=10^(-5.5),y=0.7,label="extinct detritivores",family="serif",size=6)+
      theme+theme(legend.position="none")+
      palette_detritivore_spectrum$scale_colour_line_species+
      scale_x_log10(minor_breaks = 10^seq(-10,5,1),
                    breaks = 10^(seq(-8,2,2)),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      xlab(label_bodymass)+
      ylab(label_TL)
#ggsave(paste(path_figures,"supp_TL.pdf",sep=""), p2, width = 7, height = 4, device=cairo_pdf)

# legend
p2<-ggplot(data=data_detritivore_spectrum$data_species)+
  geom_bar(aes(trophic_group,fill=names),stat="count")+
  theme_void()+theme(legend.position="none",
                     axis.text.y=element_text(family="serif",size=15,hjust=1))+
  palette_detritivore_spectrum$scale_colour_bar_species+
  coord_flip()+
  ylab(label_biomass_fraction)

graph<-ggdraw(xlim = c(0, 1.4), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0.15, 0.4, 0.7)
ggsave(paste(path_figures,"supp_TL.pdf",sep=""), graph, width = 10, height = 5, device=cairo_pdf)

# Microbes activity # ----
data<-aggregate(active~trophic_group+trophic_group_file, data=data_species, mean)
temp<-aggregate(biomass~trophic_group+trophic_group_file, data=data_species, sum)
data<-merge(data,temp,by=c("trophic_group","trophic_group_file"))
data$biomass[data$biomass<1e-7]=NA
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
# rectangle highlighting the fraction of active microbes
data$rect_active<-data$biomass*data$active
data$rect_active[data$trophic_group!="microbes"]=NA
# text with the  fraction of active microbes
data$text<-paste0("Active (",round(data$active,3)*100,"%)")
data$text[data$trophic_group!="microbes"]=NA
# position of the label
data$text_position<-10^((log10(min(data$biomass,na.rm=TRUE)*0.5)+log10(data$biomass*data$active))/2)
data$text_position[data$trophic_group!="microbes"]=NA

# p1<-ggplot(data)+
#   geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
#   theme+
#   palette_multi_channel$scale_colour_bar_trophic_groups
# legend<-get_legend(p1)

p1<-ggplot(data)+
      geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
      geom_rect(aes(xmin=0.5,xmax=1.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=rect_active),colour="black",fill=alpha(NA,1),linetype="dashed")+
      geom_text(aes(x=1,y=text_position,label=text))+
      facet_wrap(~trophic_group_file)+
      theme+theme(legend.position="none",
                  axis.title.y=element_blank(),
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  panel.grid.minor.y=element_blank(),
                  panel.grid.major.y=element_blank())+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      scale_x_continuous(labels=levels(data$trophic_group), breaks=c(1:length(levels(data$trophic_group))))+
      coord_flip()+
      scale_y_log10(breaks = 10^seq(-1,3,1),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)))+
      ylab(label_biomass)

data$relative[data$trophic_group_file==label_multichannel]<-round(data$biomass[data$trophic_group_file==label_multichannel]*data$active[data$trophic_group_file==label_multichannel]/
                                                                    sum(data$biomass[data$trophic_group_file==label_multichannel]*data$active[data$trophic_group_file==label_multichannel]),3)*100
data$relative[data$trophic_group_file==label_size_structured]<-round(data$biomass[data$trophic_group_file==label_size_structured]*data$active[data$trophic_group_file==label_size_structured]/
                                                                       sum(data$biomass[data$trophic_group_file==label_size_structured]*data$active[data$trophic_group_file==label_size_structured]),3)*100

p2<-ggplot(data=data)+
      geom_bar(aes(1,biomass*active,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      facet_wrap(~trophic_group_file)+
      theme+theme(#legend.position="none",
                  axis.title.x=element_blank(),
                  axis.text.x=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      scale_y_continuous(labels = scales::percent)+
      ylab(label_biomass_fraction_active)

graph<-ggdraw(xlim = c(0,1), ylim = c(0,2)) +
  draw_plot(p1, 0.05, 1, 0.95, 1)+
  draw_plot(p2, 0, 0, 1, 1)+
  #draw_plot(legend, 1.6, 0.25, 0.2, 0.5)+
  draw_plot_label(c("A","B"), c(0,0), c(2,1), size = 30)
ggsave(paste(path_figures,"supp_active.pdf",sep=""), graph, width = 10, height = 8, device=cairo_pdf)

# Decomposition # ----
# detritus stocks
data<-aggregate(biomass~type+input_FOM+input_DOC, data=data_size_structured$data_detritus, sum)
p1<-ggplot(data)+
      geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(biomass)*0.5,ymax=biomass,fill=type))+
      theme+theme(legend.position="none",
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_size_structured$scale_colour_bar_detritus_classes+
      scale_x_continuous(breaks=1:length(levels(data$type)),
                         labels=levels(data$type))+
      y_axis_log10+
      ylab(label_stock)

# decomposition by each decomposer trophic group
data<-data_size_structured$data_species
consumption<-read.table(paste(path_results,"consumption_",data$simu_ID[1],".txt",sep=""),sep=",",header=T)
names(consumption)<-data_size_structured$names_tot
consumption<-consumption[,data_size_structured$names_detritus] # selects abiotic resources only
consumption$consumer<-data_size_structured$names_species
consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
# resource type
temp<-data_size_structured$detritus_traits[,c("type","names")]
names(temp)<-c("type","resource")
consumption<-merge(consumption,temp,by="resource")
# consumer trophic group
temp<-data_size_structured$species_traits[,c("trophic_group","names")]
names(temp)<-c("trophic_group","consumer")
consumption<-merge(consumption,temp,by="consumer")
rm(temp)
# only consider decomposers
consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores","trophic whales"),]
consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores"))
consumption$type<-factor(consumption$type,levels=c("faeces","FOM","SOM","DOC"))

p2<-ggplot(data=consumption)+
      geom_bar(aes(type,biomass,fill=trophic_group),stat="identity",position=position_stack(reverse = TRUE))+
      theme+theme(#legend.position="top",
        legend.position=c(1, 1),
        legend.justification=c(1, 1),
        axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())+
      palette_size_structured$scale_colour_bar_trophic_groups_detritivores+
      ylab(label_decomposition)

# final graph
graph<-ggdraw(xlim = c(0, 2.1), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1.1, 0, 1, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figures,"supp_detritus.pdf",sep=""), graph, width = 11, height = 4, device=cairo_pdf)

# Diet detritivores # ----
data<-data_species[data_species$trophic_group_file==label_multichannel,]
diet<-read.table(paste(path_results,"diet_",data$simu_ID[1],".txt",sep=""),sep=",",header=T)
diet$pred<-data_multi_channel$names_species # add the names of the predators
diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
data<-data[order(data$names),]
diet[data$TL==0,c(1:(ncol(diet)-1))]=NA # removes the diet of extinct carnivores
diet<-diet[which(data$trophic_type%in%c("detritivores")),] # keeps only consumers
diet<-Rmake_table_for_plot(diet,"pred","prey","flow") # makes diet a long table
levels(diet$prey)<-data_multi_channel$names_tot
survivor<-data$names[data$TL>0 & data$trophic_type%in%c("detritivores")]
diet<-diet[diet$pred%in%survivor,]

palette_survivors<-Rmake_colour_palette(data_multi_channel$species_traits,0)[which(data_multi_channel$names_species%in%survivor)]
scale_colour_species_survivors<-scale_colour_manual(values=palette_survivors,
                                                    name=NULL)

p1<-ggplot(data=diet)+
      geom_bar(aes(pred,flow,fill=prey),stat="identity",position="fill")+
      geom_point(aes(pred,-0.05,colour=pred),size=4)+
      theme+theme(legend.position="none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.line.y=element_blank(),
                  panel.grid=element_blank())+
      scale_colour_species_survivors+
      palette_multi_channel$scale_colour_bar_total+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      coord_flip()+
      xlab("Detritivore")+
      ylab(label_diet)

# ### Diet high FOM input
# path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
# data_multi_channel<-Rload_data(path_results,1)
# palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
#                                       data_multi_channel$n_tot_species,
#                                       data_multi_channel$n_faeces)
# params<-read.table(paste(path_results,"parameters_simulation.txt",sep=""),sep=",",header=T)
# data<-data_multi_channel$data_species
# data<-data[which(data$input_DOC==100 & data$input_FOM==600),]
# diet<-read.table(paste(path_results,"diet_",params$simu_ID[which(params$input_DOC==100 & params$input_FOM==600)],".txt",sep=""),sep=",",header=T)
# diet$pred<-data_multi_channel$names_species # add the names of the predators
# diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
# data<-data[order(data$names),]
# diet[data$TL==0,c(1:(ncol(diet)-1))]=NA # removes the diet of extinct carnivores
# diet<-diet[which(data$trophic_type%in%c("microbivores","carnivores")),c(1:data_multi_channel$n_tot_species,ncol(diet))] # keeps only consumers
# diet<-Rmake_table_for_plot(diet,"pred","prey","flow") # makes diet a long table
# levels(diet$prey)<-data_multi_channel$names_species
# survivor<-data$names[data$TL>0 & data$trophic_type%in%c("microbivores","carnivores")]
# diet<-diet[diet$pred%in%survivor,]
# 
# palette_survivors<-Rmake_colour_palette(data_multi_channel$species_traits,0)[which(data_multi_channel$names_species%in%survivor)]
# scale_colour_species_survivors<-scale_colour_manual(values=palette_survivors,
#                                                     name=NULL)
# n_micro<-length(which(data$names[data$trophic_group%in%c("micro-food web microbivores","micro-food web carnivores")]%in%survivor))
# n_macro<-length(which(data$names[data$trophic_group%in%c("macro-food web carnivores")]%in%survivor))
# 
# ggplot(data=diet)+
#   geom_bar(aes(pred,flow,fill=prey),stat="identity",position="fill")+
#   geom_point(aes(pred,-0.05,colour=pred),size=4)+
#   annotate("errorbarh",y=1.03,xmin=0.7,xmax=n_micro+0.3,height = 0.05)+
#   annotate("errorbarh",y=1.03,xmin=n_micro+0.7,xmax=n_micro+n_macro+0.45,height = 0.05)+
#   annotate("text",x=n_micro/2,y=1.08,label="micro-food web",family="serif",size=6,angle=-90)+
#   annotate("text",x=n_micro+n_macro/2,y=1.08,label="macro-food web",family="serif",size=6,angle=-90)+
#   theme+theme(legend.position="none",
#               axis.text.y=element_blank(),
#               axis.ticks.y=element_blank(),
#               axis.line.y=element_blank(),
#               panel.grid.major.y=element_blank())+
#   scale_colour_species_survivors+
#   palette_multi_channel$scale_colour_bar_species+
#   scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
#   coord_flip()+
#   xlab("Predator")+
#   ylab(label_diet)

ggsave(paste(path_figures,"supp_diet.pdf",sep=""), p1, width = 6, height = 6, device=cairo_pdf)

# Energy channels # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)
params<-read.table(paste(path_results,"parameters_simulation.txt",sep=""),sep=",",header=T)
# FOM gradient
input_DOC=150
input_FOM=unique(params$input_FOM)
channel<-matrix(0,length(input_FOM),data_multi_channel$n_tot_species)
for (i in 1:length(input_FOM)){
  data<-data_multi_channel$data_species[which(data_multi_channel$data_species$input_DOC==input_DOC & data_multi_channel$data_species$input_FOM==input_FOM[i]),]
  diet<-read.table(paste(path_results,"diet_",params$simu_ID[which(params$input_DOC==input_DOC & params$input_FOM==input_FOM[i])],".txt",sep=""),sep=",",header=T)
  diet$tot<-rowSums(diet)
  diet$pred<-data_multi_channel$names_species # add the names of the predators
  diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
  data<-data[order(data$names),]
  diet[data$TL==0,c(1:(ncol(diet)-1))]=0 # removes the diet of extinct carnivores
  diet<-diet[,c(1:data_multi_channel$n_tot_species,ncol(diet)+c(-1,0))] # keeps only consumers
  fraction=rep(0,data_multi_channel$n_tot_species) # energy from microbes = 0
  fraction[data$trophic_type=="detritivores"]=1 # energy from detritivores = 1
  
  for (j in which(data$trophic_type%in%c("microbivores","carnivores"))){
    if (diet$tot[j]>0){
      fraction[j]=sum(diet[j,c(1:data_multi_channel$n_tot_species)]*fraction)/diet$tot[j]
    }
  }
  fraction[data$TL==0]=NA
  channel[i,]<-fraction
}
channel<-as.data.frame(channel)
names(channel)<-data_multi_channel$names_species
survivor<-data_multi_channel$species_traits$names[data_multi_channel$species_traits$names!="micro-food web carnivores 1" & data_multi_channel$species_traits$trophic_type%in%c("microbivores","carnivores")]
channel<-channel[,names(channel)%in%survivor]
channel$input_DOC=input_DOC
channel$input_FOM=input_FOM
channel<-Rmake_table_for_plot(channel,c("input_DOC","input_FOM"),"pred","fraction") # makes diet a long table

palette_survivors<-Rmake_colour_palette(data_multi_channel$species_traits,0)[which(data_multi_channel$names_species%in%survivor)]
scale_colour_species_survivors<-scale_colour_manual(values=palette_survivors,
                                                    name=NULL)

n_micro<-length(which(data_multi_channel$species_traits$names[data_multi_channel$species_traits$trophic_group%in%c("micro-food web microbivores","micro-food web carnivores")]%in%survivor))
n_macro<-length(which(data_multi_channel$species_traits$names[data_multi_channel$species_traits$trophic_group%in%c("macro-food web carnivores")]%in%survivor))

p1<-ggplot(data=channel)+
  geom_raster(aes(input_FOM,pred,fill=fraction))+
  theme_raster+
  scale_fill_gradient(low="hotpink2",
                      high="purple2",
                      name="Fraction of energy\nfrom detritivore\nchannel",
                      limits=c(0,1),
                      breaks=seq(0,1,0.25),
                      labels = scales::percent)
legend<-get_legend(p1)

p1<-ggplot(data=channel)+
      geom_raster(aes(input_FOM,pred,fill=fraction))+
      geom_point(aes(-0.05,as.numeric(pred),colour=pred),size=4)+
      annotate("errorbar",x=max(input_FOM)+50,ymin=0.7,ymax=n_micro+0.3,width = 30)+
      annotate("errorbar",x=max(input_FOM)+50,ymin=n_micro+0.7,ymax=n_micro+n_macro+0.45,width = 30)+
      annotate("text",x=max(input_FOM)+80,y=n_micro/2,label="micro-food web",family="serif",size=6,angle=-90)+
      annotate("text",x=max(input_FOM)+80,y=n_micro+n_macro/2,label="macro-food web",family="serif",size=6,angle=-90)+
      theme_raster+theme(legend.position="none",
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         axis.line.y=element_blank(),
                         panel.grid.major.y=element_blank())+
      scale_colour_species_survivors+
      scale_fill_gradient(low="hotpink2",
                          high="purple2",
                         name="Fraction of energy\nfrom detritivore\nchannel",
                         limits=c(0,1),
                         breaks=seq(0,1,0.25),
                         labels = scales::percent)+
      guides(colour = "none")+
      xlab(label_FOM_input)+
      ylab("Predator")

# DOC gradient
input_DOC=unique(params$input_DOC)
input_FOM=300
channel<-matrix(0,length(input_DOC),data_multi_channel$n_tot_species)
for (i in 1:length(input_DOC)){
  data<-data_multi_channel$data_species[which(data_multi_channel$data_species$input_DOC==input_DOC[i] & data_multi_channel$data_species$input_FOM==input_FOM),]
  diet<-read.table(paste(path_results,"diet_",params$simu_ID[which(params$input_DOC==input_DOC[i] & params$input_FOM==input_FOM)],".txt",sep=""),sep=",",header=T)
  diet$tot<-rowSums(diet)
  diet$pred<-data_multi_channel$names_species # add the names of the predators
  diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
  data<-data[order(data$names),]
  diet[data$TL==0,c(1:(ncol(diet)-1))]=0 # removes the diet of extinct carnivores
  diet<-diet[,c(1:data_multi_channel$n_tot_species,ncol(diet)+c(-1,0))] # keeps only consumers
  fraction=rep(0,data_multi_channel$n_tot_species) # energy from microbes = 0
  fraction[data$trophic_type=="detritivores"]=1 # energy from detritivores = 1
  
  for (j in which(data$trophic_type%in%c("microbivores","carnivores"))){
    if (diet$tot[j]>0){
      fraction[j]=sum(diet[j,c(1:data_multi_channel$n_tot_species)]*fraction)/diet$tot[j]
    }
  }
  fraction[data$TL==0]=NA
  channel[i,]<-fraction
}
channel<-as.data.frame(channel)
names(channel)<-data_multi_channel$names_species
survivor<-data_multi_channel$species_traits$names[data_multi_channel$species_traits$names!="micro-food web carnivores 1" & data_multi_channel$species_traits$trophic_type%in%c("microbivores","carnivores")]
channel<-channel[,names(channel)%in%survivor]
channel$input_DOC=input_DOC
channel<-channel[channel$input_DOC>10,]
channel$input_FOM=input_FOM
channel<-Rmake_table_for_plot(channel,c("input_DOC","input_FOM"),"pred","fraction") # makes diet a long table

palette_survivors<-Rmake_colour_palette(data_multi_channel$species_traits,0)[which(data_multi_channel$names_species%in%survivor)]
scale_colour_species_survivors<-scale_colour_manual(values=palette_survivors,
                                                    name=NULL)

n_micro<-length(which(data_multi_channel$species_traits$names[data_multi_channel$species_traits$trophic_group%in%c("micro-food web microbivores","micro-food web carnivores")]%in%survivor))
n_macro<-length(which(data_multi_channel$species_traits$names[data_multi_channel$species_traits$trophic_group%in%c("macro-food web carnivores")]%in%survivor))

p2<-ggplot(data=channel)+
      geom_raster(aes(input_DOC,pred,fill=fraction))+
      geom_point(aes(-0.05,as.numeric(pred),colour=pred),size=4)+
      annotate("errorbar",x=max(input_DOC)+50,ymin=0.7,ymax=n_micro+0.3,width = 30)+
      annotate("errorbar",x=max(input_DOC)+50,ymin=n_micro+0.7,ymax=n_micro+n_macro+0.45,width = 30)+
      annotate("text",x=max(input_DOC)+70,y=n_micro/2,label="micro-food web",family="serif",size=6,angle=-90)+
      annotate("text",x=max(input_DOC)+70,y=n_micro+n_macro/2,label="macro-food web",family="serif",size=6,angle=-90)+
      theme_raster+theme(legend.position="none",
                         axis.text.y=element_blank(),
                         axis.ticks.y=element_blank(),
                         axis.line.y=element_blank(),
                         panel.grid.major.y=element_blank())+
      scale_colour_species_survivors+
      scale_fill_gradient(low="hotpink2",
                          high="purple2",
                          name="Fraction of energy\nfrom detritivore\nchannel",
                          limits=c(0,1),
                          breaks=seq(0,1,0.25),
                          labels = scales::percent)+
      guides(colour = "none")+
      xlab(label_DOC_input)+
      ylab("Predator")

graph<-ggdraw(xlim = c(0, 2.5), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2.15, 0.25, 0.15, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figures,"supp_channel.pdf",sep=""), graph, width = 12, height = 5, device=cairo_pdf)

# Time series # ----
# multichannel
data_TS<-read.table(paste(path_results,"TS_density_",1,".txt",sep=""),sep=",",header=T)
data_TS<-data_TS[,1:(data_multi_channel$n_tot_species+1)]
names(data_TS)[2:ncol(data_TS)]<-data_multi_channel$names_species
#data_TS<-data_TS[data_TS$c0_prey>5e3 & data_TS$c0_pred>5e4,]
data_TS<-Rmake_table_for_plot(data_TS,c("time"),"names","biomass")
data_TS$trophic_group_file=1
# size-structured
temp<-read.table(paste(path_results,"TS_density_",3,".txt",sep=""),sep=",",header=T)
temp<-temp[,1:(data_multi_channel$n_tot_species+1)]
names(temp)[2:ncol(temp)]<-data_multi_channel$names_species
#temp<-temp[temp$c0_prey>5e3 & temp$c0_pred>5e4,]
temp<-Rmake_table_for_plot(temp,c("time"),"names","biomass")
temp$trophic_group_file=3
# merge
data_TS<-rbind(data_TS,temp)
data_TS$trophic_group_file<-as.factor(data_TS$trophic_group_file)
levels(data_TS$trophic_group_file)=c(label_multichannel,label_size_structured)
data_TS$biomass[data_TS$biomass<1e-6]=NA
# time range for plot x_axis
time_range<-max(data_TS$time)-min(data_TS$time)
time_range<-floor(min(data_TS$time)+time_range*c(0.25,0.75))
close_order<-floor(log10(time_range))-1
time_range<-round(time_range*10^(-close_order))*10^close_order

p1<-ggplot(data_TS)+
  geom_line(aes(time,biomass,colour=names),linewidth=1.5)+
  facet_wrap(~trophic_group_file)+
  theme+theme(legend.position="none",
              panel.spacing.x = unit(2, "lines"),
              plot.margin = margin(t=5.5,r=15,b=5.5,l=5.5, "pt"))+ # 5.5 is the default value
  palette_multi_channel$scale_colour_line_species+
  scale_x_continuous(breaks = time_range,
                     name = label_time)+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                name = label_biomass)
ggsave(paste(path_figures,"supp_time_series.pdf",sep=""), p1, width = 8, height = 4, device=cairo_pdf)

# Biomass empiric #----
data_abundance<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaAbundanceData.csv",header=T,sep=",")
data_bodymass<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv",header=T,sep=",")
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
Mkm2_to_m2=1e6*1e6 # million km² to km²
Pg_to_mg=1e18 # petagram to milligram

data_bodymass$Individual_bodymass_mg_FM<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry
data_bodymass$trophic_group<-"macro-food web"
data_bodymass$trophic_group[data_bodymass$Individual_bodymass_mg_FM<0.001]="micro-food web"

data_abundance<-merge(data_abundance,data_bodymass,by=c("Soil_biota_group"))
data_abundance$Individual_bodymass_mg_DM<-data_abundance$Individual_bodymass_mg_DM*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance$biomass_density<-data_abundance$Biomass_g_m2*1e3*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance<-data_abundance[,c("Soil_biota_group","Biome_ecosystem","Biomass_g_m2","Individual_bodymass_mg_DM","trophic_group","biomass_density")]

# bacteria<-data.frame(Soil_biota_group="Bacteria",
#                      Biome_ecosystem=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),
#                      Biomass_g_m2=c(22.8,23.5,22.6,22.4,22.3), # ln(Individuals) actually
#                      Individual_bodymass_mg_DM=2.83e-8*dry_to_C, # converted into C body mass
#                      trophic_group="microbes")
# bacteria$biomass_density<-exp(bacteria$Biomass_g_m2)*bacteria$Individual_bodymass_mg_DM

microbes<-read.table("Data/Metabolic_rates_Johnston_2018/Xu_2013_Global_Ecology_and_Biogeography_aggregated_data.csv",header=T,sep=",") # data from the main text of Xu et al 2013 Global Ecology and Biogeography
microbes$biomass_density<-microbes$soil_microbial_C_Pg_0_30_cm*Pg_to_mg/(microbes$Area_Mkm2*Mkm2_to_m2)
microbes$Soil_biota_group="microbes"
microbes$trophic_group="microbes"

data_abundance<-rbind(data_abundance[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")],
                      microbes[,c("Soil_biota_group","Biome_ecosystem","trophic_group","biomass_density")])
data_abundance$trophic_group<-factor(data_abundance$trophic_group,levels=c("microbes","micro-food web","macro-food web"))
data_abundance<-data_abundance[data_abundance$Biome_ecosystem%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data_abundance$Biome_ecosystem<-factor(data_abundance$Biome_ecosystem,levels=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"))
levels(data_abundance$Biome_ecosystem)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")

data<-aggregate(biomass_density~Biome_ecosystem+trophic_group+Soil_biota_group,data_abundance,mean)
data<-aggregate(biomass_density~Biome_ecosystem+trophic_group,data,sum)
p1<-ggplot(data=data)+
      geom_col(aes(Biome_ecosystem,biomass_density,fill=trophic_group),position="dodge")+
      theme+theme(legend.title=element_blank(),
                  axis.title.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                        guide = guide_legend(reverse=TRUE))+
      y_axis_log10+
      xlab("Biome")+
      ylab(label_biomass)

p2<-ggplot(data=data)+
      geom_bar(aes(Biome_ecosystem,biomass_density,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      theme+theme(legend.title=element_blank(),
                  axis.title.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                        guide = guide_legend(reverse=TRUE))+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction_2lines)

# Heděnec et al 2022 Scientific Reports
data_hedenec<-read.table(paste(path_hedenec,"biomass_invertebrates_raw.csv",sep=""),sep=",",header=T)
data_hedenec_trophic_groups<-read.table(paste(path_hedenec,"invertebrates_trophic_groups.csv",sep=""),sep=",",header=T)
data_hedenec<-merge(data_hedenec,data_hedenec_trophic_groups,by=c("Guilds","Fauna"))
rm(data_hedenec_trophic_groups)

data<-aggregate(Biomass~trophic_type+size_class+Biom, data=data_hedenec, sum)
data$size_class<-factor(data$size_class,levels=c("micro-foodweb","macro-foodweb"))
levels(data$size_class)=c("micro","macro")
data$Biom<-factor(data$Biom,levels=c("Polar region","Taiga","Desert","Mediterranean vegetation","Temperate forest",
                                     "Temperate grassland","Tropical forest","Tropical grassland"))
levels(data$Biom)[4]="Mediterranean"
data$trophic_type<-factor(data$trophic_type,levels=c("carnivores","herbivores","detritivores","omnivores","microbivores"))

# p3<-ggplot(data=data)+
#       geom_bar(aes(1,Biomass,fill=trophic_type),stat="identity")+
#       facet_wrap(~Biom, nrow=2)+
#       theme+theme(legend.title=element_blank(),
#                   axis.ticks.x=element_blank(),
#                   panel.grid.minor.x=element_blank(),
#                   panel.grid.major.x=element_blank())+
#       scale_fill_manual(values=c("red2","olivedrab3","purple2","lightsalmon2","seagreen2"))+
#       xlab("Food web")+
#       ylab(label_biomass_fraction)

p3<-ggplot(data=data)+
      geom_bar(aes(size_class,Biomass,fill=trophic_type),stat="identity",position=position_fill())+
      facet_wrap(~Biom, nrow=2)+
      theme+theme(legend.title=element_blank(),
                  axis.ticks.x=element_blank(),
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      scale_fill_manual(values=c("red2","olivedrab3","purple2","lightsalmon2","seagreen2"))+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      xlab("Food web")+
      ylab(label_biomass_fraction)

graph<-ggdraw(xlim = c(0, 1), ylim = c(0, 2.8)) +
  draw_plot(p1, 0, 1.8, 1, 1)+
  draw_plot(p2, 0, 1, 1, 0.8)+
  draw_plot(p3, 0, 0, 1, 1)+
  draw_plot_label(c("A","B","C"), c(0,0,0), c(2.8,1.8,1.1), size = 30)
ggsave(paste(path_figures,"supp_biomass_empiric.pdf",sep=""), graph, width = 11, height = 12, device=cairo_pdf)

# Biomass FOM - DOC # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

data<-data_multi_channel$data_species
data<-aggregate(biomass~trophic_group+input_FOM+input_DOC, data=data, sum, na.rm = TRUE)
#data$biomass[data$biomass<1e-7]=NA
data$foodweb<-"micro-food web"
data$foodweb[data$trophic_group%in%c("macro-food web detritivores","macro-food web carnivores","trophic whales")]<-"macro-food web"
input_DOC<-unique(data$input_DOC)
input_FOM<-unique(data$input_FOM)
data_FOM_DOC<-read.table(paste0(path_data,"data_FOM_SOM.csv"),header=T,sep=",") # data of FOM and DOC inputs in the main land ecosystems
for (i in 1:nrow(data_FOM_DOC)){
  data_FOM_DOC$input_DOC[i]<-input_DOC[which.min(abs(input_DOC-data_FOM_DOC$input_DOC[i]))] # place the values by the nearest values used for the simulations
  data_FOM_DOC$input_FOM[i]<-input_FOM[which.min(abs(input_FOM-data_FOM_DOC$input_FOM[i]))]
}
data<-merge(data,data_FOM_DOC,by=c("input_DOC","input_FOM"))
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
data$Biome<-factor(data$Biome,levels=c("Desert","Tundra","Boreal forest","Temperate forest","Mediterranean vegetation","Temperate grassland","Tropical grassland","Tropical forest"))

p1<-ggplot(data=data)+
      geom_bar(aes(Biome,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
      theme+theme(axis.title.y=element_blank(),
                  panel.grid=element_blank())+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      coord_flip()+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      ylab(label_biomass_fraction)
ggsave(paste(path_figures,"supp_biomass_fraction_DOC_FOM.pdf",sep=""), p1, width = 9, height = 4, device=cairo_pdf)

# Decomposition FOM - DOC # ----
path_results<-"results_final_FOM_SOM/"#"results_23_02_06/"
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

params_data<-read.table(paste(path_results,"parameters_simulation.txt",sep=""),sep=",",header=T)
data<-NULL
for (i in 1:nrow(params_data)){
  consumption<-read.table(paste(path_results,"consumption_",params_data$simu_ID[i],".txt",sep=""),sep=",",header=T)
  names(consumption)<-data_multi_channel$names_tot
  consumption<-consumption[,data_multi_channel$names_detritus] # selects abiotic resources only
  consumption$consumer<-data_multi_channel$names_species
  consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
  # resource type
  temp<-data_multi_channel$detritus_traits[,c("type","names")]
  names(temp)<-c("type","resource")
  consumption<-merge(consumption,temp,by="resource")
  # consumer trophic group
  temp<-data_multi_channel$species_traits[,c("trophic_group","names")]
  names(temp)<-c("trophic_group","consumer")
  consumption<-merge(consumption,temp,by="consumer")
  rm(temp)
  # only consider decomposers
  consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores"),]
  consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
  consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores"))
  consumption$type<-factor(consumption$type,levels=c("faeces","FOM","SOM","DOC"))
  consumption$input_FOM<-params_data$input_FOM[i]
  consumption$input_DOC<-params_data$input_DOC[i]
  data<-rbind(data,consumption)
}

data<-dcast(data,type+input_FOM+input_DOC~trophic_group,value.var="biomass")
data$ratio<-data$'macro-food web detritivores'/(data$microbes+data$'macro-food web detritivores')
input_DOC<-unique(data$input_DOC)
input_FOM<-unique(data$input_FOM)
data_point<-expand_grid(input_DOC=c(input_DOC[1],input_DOC[length(input_DOC)]),
                        input_FOM=c(input_FOM[1],input_FOM[length(input_FOM)]))
data$input_DOC[data$input_DOC==input_DOC[1]]=0 # to have a clean raster

data_FOM_DOC<-read.table(paste0(path_data,"data_FOM_SOM.csv"),header=T,sep=",")
data_FOM_DOC$type="empiric"
data_FOM_DOC<-rbind(data_FOM_DOC,list("Figure 2",300,150,"model"))

p1<-ggplot(data=data[data$type=="FOM",])+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  theme_raster+
  scale_fill_viridis(option="viridis",
                     name="Macro-food web\nfraction of\ndecomposition",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)
legend<-get_legend(p1)

p1<-ggplot(data=data[data$type=="FOM",])+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  geom_point(data=data_FOM_DOC,aes(input_DOC,input_FOM,colour=type),size=5,shape=18)+
  coord_cartesian(clip = "off")+ # allow the labels to go beyond the edges of the panel
  geom_label_repel(data=data_FOM_DOC,aes(input_DOC,input_FOM,label=Biome,colour=type),
                   point.padding=1,size=4,
                   fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+ # allow the labels to go beyond the edges of the panel
  #geom_point(data=data_point,aes(input_DOC,input_FOM),colour="red",size=5)+
  theme_raster+theme(legend.position="none")+
  scale_fill_viridis(option="viridis",
                     name="Macro-food web\nfraction of\ndecomposition",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)+
  scale_colour_data_type+
  xlab(label_DOC_input)+
  ylab(label_FOM_input)+
  ggtitle("FOM")

p2<-ggplot(data=data[data$type=="faeces",])+
  geom_raster(aes(input_DOC,input_FOM,fill=ratio))+
  geom_point(data=data_FOM_DOC,aes(input_DOC,input_FOM,colour=type),size=5,shape=18)+
  coord_cartesian(clip = "off")+ # allow the labels to go beyond the edges of the panel
  geom_label_repel(data=data_FOM_DOC,aes(input_DOC,input_FOM,label=Biome,colour=type),
                   point.padding=1,size=4,
                   fill = "white", xlim = c(-Inf, Inf), ylim = c(-Inf, Inf))+ # allow the labels to go beyond the edges of the panel
  #geom_point(data=data_point,aes(input_DOC,input_FOM),colour="red",size=5)+
  theme_raster+theme(legend.position="none")+
  scale_fill_viridis(option="viridis",
                     name="Macro-food web\nfraction of\ndecomposition",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)+
  scale_colour_data_type+
  xlab(label_DOC_input)+
  ylab(label_FOM_input)+
  ggtitle("faeces")

graph<-ggdraw(xlim = c(0,2.4), ylim = c(0,1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(legend, 2, 0.25, 0.4, 0.5)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figures,"supp_decomposition_ratio.pdf",sep=""), graph, width = 12, height = 5, device=cairo_pdf)

# Respiration empiric # ----
### Calculation of metabolic parameters
# data from Johnston and Sibly (2018) (10.1038/s41559-018-0648-6)
joule_to_carbon=1/20.1*0.5363
hour_to_day=1/24
celsius_to_kelvin=273.15
kB=8.62e-5 # Bolzmann's constant
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
# metabolic rate : J.h-1
# body mass : mg fresh weight
data_metabo<-read.csv(file=paste(path_data,"Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaMetabolicData.csv",sep=""), header=TRUE)
names(data_metabo)<-c("study","Fauna_group","taxon","metabolic_rate","bodymass","temperature")

# original values before conversion (same results as Johnston)
data_metabo$temperature_arrhenius<--1/(celsius_to_kelvin+data_metabo$temperature)*1/kB
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~Fauna_group+log(bodymass):Fauna_group+temperature_arrhenius:Fauna_group)
par(mfrow=c(2,2))
plot(model)
summary(model)
Anova(model,type=2,test.statisticf="F")

ggplot(data=data_metabo)+
  geom_point(aes(bodymass,metabolic_rate/bodymass,colour=Fauna_group))+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE))+
  theme+theme(legend.title=element_blank())+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Fresh body mass (mg)"))+
  ylab(expression("Respiration rate (J hr"^-1*")"))

# conversion
data_metabo$metabolic_rate<-data_metabo$metabolic_rate*joule_to_carbon/hour_to_day
data_metabo$bodymass<-data_metabo$bodymass*fresh_to_dry*dry_to_C
data_metabo$metabolic_rate<-data_metabo$metabolic_rate/data_metabo$bodymass # conversion into mass metabolic rate
model=lm(data=data_metabo,
         formula=log(metabolic_rate)~Fauna_group+log(bodymass):Fauna_group+temperature_arrhenius:Fauna_group)
par(mfrow=c(2,2))
plot(model)
params<-summary(model)
params
Anova(model,type=2,test.statisticf="F")
params<-(params$coefficients)[,1]
names(params)[1]="Fauna_groupMacrofauna"
params[2]=params[2]+params[1] # total effect for mesofauna
params[3]=params[3]+params[1] # total effect for microbes
params<-as.data.frame(params)
params$names<-row.names(params)
params<-params %>% separate(names,c("Fauna_group","parameter"),sep=":",fill="right")
params<-params %>% separate(Fauna_group,c(NA,"Fauna_group"),sep=11)
params$parameter[is.na(params$parameter)]="R0"
params$parameter<-as.factor(params$parameter)
levels(params$parameter)=c("s_R","R0","E_R")
params<-dcast(params,Fauna_group~parameter,value.var="params")

### Calculation of the respiration of each sub-food web
data_abundance<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilBiotaAbundanceData.csv",header=T,sep=",")
data_bodymass<-read.table("Data/Metabolic_rates_Johnston_2018/JohnstonSiblySoilGroupIndividualBodymass.csv",header=T,sep=",")
Mkm2_to_m2=1e6*1e6 # million km² to km²
Pg_to_mg=1e18 # petagram to milligram

# mean body mass of the main soil taxa
data_bodymass$bodymass<-data_bodymass$Individual_bodymass_mg_DM*dry_to_C # body
data_bodymass$bodymass_fresh<-data_bodymass$Individual_bodymass_mg_DM/fresh_to_dry
data_bodymass$trophic_group<-"macro-food web"
data_bodymass$trophic_group[data_bodymass$bodymass_fresh<0.001]="micro-food web"
data_bodymass$trophic_group[data_bodymass$Soil_biota_group=="Bacteria"]<-"microbes"
data_bodymass$Individual_bodymass_mg_DM<-NULL
# abundance of the main soil fauna taxa
data_abundance<-data_abundance[,c("Biome_ecosystem","Soil_biota_group","Biomass_g_m2")]
data_abundance$biomass_density<-data_abundance$Biomass_g_m2*1e3*fresh_to_dry*dry_to_C # conversion into mgC
data_abundance$Biomass_g_m2<-NULL
# abundance of soil microbes
microbes<-read.table("Data/Metabolic_rates_Johnston_2018/Xu_2013_Global_Ecology_and_Biogeography_aggregated_data.csv",header=T,sep=",") # data from the main text of Xu et al 2013 Global Ecology and Biogeography
microbes$biomass_density<-microbes$soil_microbial_C_Pg_0_30_cm*Pg_to_mg/(microbes$Area_Mkm2*Mkm2_to_m2)
microbes$Soil_biota_group="Bacteria"
microbes<-microbes[,names(data_abundance)]

# final data set
data<-rbind(data_abundance,microbes)
rm(data_abundance,microbes)
data<-merge(data,data_bodymass,by=c("Soil_biota_group"))
data$trophic_group<-factor(data$trophic_group,levels=c("microbes","micro-food web","macro-food web"))
data<-data[data$Biome_ecosystem%in%c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"),]
data$Biome_ecosystem<-factor(data$Biome_ecosystem,levels=c("Tundra","Boreal forest","Temperate forest","Temperate grassland","Tropical forest"))
levels(data$Biome_ecosystem)<-c("Tundra","Boreal\nforest","Temperate\nforest","Temperate\ngrassland","Tropical\nforest")
data<-merge(data,params,by=c("Fauna_group")) # add the metabolic parameters

# respiration
data$metabolism<-exp(data$R0)*(data$bodymass^data$s_R)*exp(-data$E_R/(kB*(celsius_to_kelvin+15)))
data$respiration<-data$metabolism*data$biomass_density

ggplot(data=data)+
  geom_point(aes(bodymass_fresh,metabolism,colour=Fauna_group))+
  scale_colour_viridis(discrete=T,
                       guide=guide_legend(reverse=TRUE))+
  theme+theme(legend.title=element_blank())+
  x_axis_log10+
  y_axis_log10+
  xlab(expression("Fresh body mass (mg)"))+
  ylab(label_respiration)

databis<-aggregate(respiration~Biome_ecosystem+trophic_group+Soil_biota_group,data,mean)
databis<-aggregate(respiration~trophic_group,databis,sum)

ggplot(data=databis)+
  geom_col(aes(1,respiration,fill=trophic_group),stat="identity",position=position_fill(reverse=TRUE))+
  theme+theme(legend.title=element_blank(),
              axis.title.x=element_blank(),)+
  scale_fill_manual(values=c("lightgrey","darkgrey","black"),
                    guide = guide_legend(reverse=TRUE))+
  xlab("Biome")+
  ylab(label_respiration)

# Detritus size # ----
data_detritus<-data_multi_channel$data_detritus

ggplot(data=data_detritus)+
  geom_point(aes(names,bodymass,colour=names),size=2)+
  palette_multi_channel$scale_colour_line_detritus+
  theme+
  y_axis_log10

ggplot(data=data_detritus)+
  geom_point(aes(names,radius,colour=names),size=2)+
  palette_multi_channel$scale_colour_line_detritus+
  theme+
  y_axis_log10

ggplot(data=data_detritus)+
  geom_rect(aes(xmin=as.numeric(names)-0.5,xmax=as.numeric(names)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=names))+
  palette_multi_channel$scale_colour_bar_detritus+
  theme+
  y_axis_log10

# SENSITIVITY ANALYSIS DECOMPOSITION # ----
### load data # ----
path_results<-"results_final_decomposition/"#"results_23_02_22_2/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_sensitivity_decomposition<-Rload_data(path_results,1)
palette_sensitivity_decomposition<-Rmake_palettes(data_sensitivity_decomposition$species_traits,
                                                  data_sensitivity_decomposition$n_tot_species,
                                                  data_sensitivity_decomposition$n_faeces)
data_species<-NULL
data_detritus<-NULL
data_TS<-NULL
for (i in 1:nrow(params)){
  # species data
  data_sensitivity_decomposition<-Rload_data(path_results,i)
  temp<-data_sensitivity_decomposition$data_species
  data_species<-rbind(data_species,temp)
  # detritus data
  temp<-data_sensitivity_decomposition$data_detritus
  data_detritus<-rbind(data_detritus,temp)
  # density time series data
  temp<-read.table(paste(path_results,"TS_density_",i,".txt",sep=""),sep=",",header=T)
  temp<-temp[,1:(data_sensitivity_decomposition$n_tot_species+1)]
  temp$a_phi_DOC=params$a_phi_DOC[i]
  temp$a_phi_FOM=params$a_phi_FOM[i]
  temp$a_phi_SOM=params$a_phi_SOM[i]
  data_TS<-rbind(data_TS,temp)
}
rm(temp)

### biomass distribution # ----
data_species$biomass[data_species$biomass<1e-7]=NA
data<-aggregate(biomass~trophic_group+a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data_species, sum)
a_phi_SOM=unique(data$a_phi_SOM)
biomass_range<-floor(log10(c(min(data$biomass,na.rm=T),max(data$biomass,na.rm=T)))) # rescale the biomass axis
biomass_range<-10^rev(seq(biomass_range[2],biomass_range[1],-2))

# biomass distribution
p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data)+
    geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
    facet_grid(a_phi_DOC~a_phi_FOM)+
    theme+theme(legend.position="none",
                axis.title.y.left=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_blank())+
    palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
    scale_x_continuous(labels=levels(data$trophic_group), breaks=c(1:length(levels(data$trophic_group))),
                       sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = biomass_range,#trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
    coord_flip()+
    ylab(label_biomass)+
    ggtitle(parse(text=paste0(label_a_phi_SOM,'*" = "*',a_phi_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_biomass.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

# update of the parameter space to remove combination without coexistence
data<-data[-which(data$a_phi_FOM==10 | (data$a_phi_FOM==100 & data$a_phi_DOC==1) | data$a_phi_SOM==100),]
data_detritus<-data_detritus[-which(data_detritus$a_phi_FOM==10 | (data_detritus$a_phi_FOM==100 & data_detritus$a_phi_DOC==1)),]
a_phi_SOM=unique(data$a_phi_SOM)

# relative biomass distribution
p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data=data[data$a_phi_SOM==a_phi_SOM[i],])+
    geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
    facet_grid(a_phi_DOC~a_phi_FOM)+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
    scale_x_continuous(sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
    scale_y_continuous(breaks=c(0,0.5,1), labels = scales::percent,
                       sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
    ylab(label_biomass_fraction)+
    ggtitle(parse(text=paste0(label_a_phi_SOM,'*" = "*',a_phi_SOM[i])))
}

graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_biomass_relative.pdf",sep=""), graph, width = 12, height = 5, device=cairo_pdf)

### distribution stock and decomposition # ----
# detritus stocks
data<-aggregate(biomass~type+a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data_detritus, sum)
biomass_range<-floor(log10(c(min(data$biomass[data$biomass>0],na.rm=T),max(data$biomass,na.rm=T)))) # rescale the biomass axis
biomass_range<-10^rev(seq(biomass_range[2],biomass_range[1],-2))
a_phi_SOM=unique(data$a_phi_SOM)

p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data=data[data$a_phi_SOM==a_phi_SOM[i],])+
    geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(biomass)*0.5,ymax=biomass,fill=type))+
    facet_grid(a_phi_DOC~a_phi_FOM)+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_decomposition$scale_colour_bar_detritus_classes+
    scale_x_continuous(breaks=1:length(levels(data$type)),
                       labels=levels(data$type),
                       sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = biomass_range,#trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
    ylab(label_stock)+
    ggtitle(parse(text=paste0(label_a_phi_SOM,'*" = "*',a_phi_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_detritus.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

# update of the parameter space to remove combination without coexistence
data_species<-data_species[-which(data_species$a_phi_FOM==10 | (data_species$a_phi_FOM==100 & data_species$a_phi_DOC==1) | data_species$a_phi_SOM==100),]
data_detritus<-data_detritus[-which(data_detritus$a_phi_SOM==100),]
data_TS<-data_TS[-which(data_TS$a_phi_FOM==10 | (data_TS$a_phi_FOM==100 & data_TS$a_phi_DOC==1) | data_TS$a_phi_SOM==100),]
a_phi_SOM=unique(data_species$a_phi_SOM)

# decomposition by each decomposer trophic group
detritus_order<-c("FOM","faeces","SOM","DOC")
consumption_temp=NULL
for (i in 1:nrow(params)){
  consumption<-read.table(paste(path_results,"consumption_",params$simu_ID[i],".txt",sep=""),sep=",",header=T)
  names(consumption)<-data_sensitivity_decomposition$names_tot
  consumption<-consumption[,data_sensitivity_decomposition$names_detritus] # selects abiotic resources only
  consumption$consumer<-data_sensitivity_decomposition$names_species
  consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
  # resource type
  temp<-data_sensitivity_decomposition$detritus_traits[,c("type","names")]
  names(temp)<-c("type","resource")
  consumption<-merge(consumption,temp,by="resource")
  # consumer trophic group
  temp<-data_sensitivity_decomposition$species_traits[,c("trophic_group","names")]
  names(temp)<-c("trophic_group","consumer")
  consumption<-merge(consumption,temp,by="consumer")
  rm(temp)
  # only consider decomposers
  consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores"),]
  consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
  consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores"))
  consumption$type<-factor(consumption$type,levels=detritus_order)
  consumption$a_phi_DOC=params$a_phi_DOC[i]
  consumption$a_phi_FOM=params$a_phi_FOM[i]
  consumption$a_phi_SOM=params$a_phi_SOM[i]
  consumption_temp<-rbind(consumption_temp,consumption)
}
consumption<-consumption_temp
rm(consumption_temp)
consumption<-consumption[-which(consumption$a_phi_FOM==10 | (consumption$a_phi_FOM==100 & consumption$a_phi_DOC==1) | consumption$a_phi_SOM==100),]

p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data=consumption[consumption$a_phi_SOM==a_phi_SOM[i],])+
            geom_bar(aes(as.numeric(type),biomass,fill=trophic_group),stat="identity")+
            facet_grid(a_phi_DOC~a_phi_FOM)+
            theme+theme(legend.position="bottom",
                        axis.title.x.bottom=element_blank(),
                        axis.text.x=element_text(angle=90,vjust=0.5),
                        panel.grid.minor.x=element_blank(),
                        panel.grid.major.x=element_blank())+
            palette_sensitivity_decomposition$scale_colour_bar_trophic_groups_detritivores+
            scale_x_continuous(breaks=1:length(levels(consumption$type)),
                               labels=levels(consumption$type),
                               sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
            scale_y_continuous(sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
            ylab(label_decomposition)+
            ggtitle(parse(text=paste0(label_a_phi_SOM,'*" = "*',a_phi_SOM[i])))
}

graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_decomposition.pdf",sep=""), graph, width = 12, height = 6, device=cairo_pdf)

### decomposition rate # ----
data<-consumption[consumption$trophic_group=="microbes",]
names(data)[names(data)=="biomass"]="decomposition"
temp<-aggregate(biomass~a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data_species[data_species$trophic_group=="microbes",], sum)
data<-merge(data,temp,by=c("a_phi_DOC","a_phi_FOM","a_phi_SOM"))
temp<-aggregate(biomass~type+a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data_detritus, sum)
names(temp)[names(temp)=="biomass"]="stock"
data<-merge(data,temp,by=c("a_phi_DOC","a_phi_FOM","a_phi_SOM","type"))
data$decomposition_rate<-data$decomposition/data$biomass
data$a_phi_SOM<-as.factor(data$a_phi_SOM)
levels(data$a_phi_SOM)<-paste0(" ",levels(data$a_phi_SOM)," ")
rm(temp)
decomposition_range<-floor(log10(c(min(data$decomposition_rate,na.rm=T),max(data$decomposition_rate,na.rm=T)))) # rescale the biomass axis
decomposition_range<-10^rev(seq(decomposition_range[2],decomposition_range[1],-1))

p1<-ggplot(data=data[data$type=="FOM",])+
      geom_rect(aes(xmin=as.numeric(a_phi_SOM)-0.45,xmax=as.numeric(a_phi_SOM)+0.45,ymin=min(decomposition_rate,na.rm=TRUE)*0.5,ymax=decomposition_rate))+
      facet_grid(a_phi_DOC~a_phi_FOM)+
      theme+theme(legend.position="none",
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
      scale_x_continuous(breaks=1:length(levels(data$a_phi_SOM)),
                         labels=levels(data$a_phi_SOM),
                         sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
      scale_y_log10(breaks = decomposition_range,#trans_breaks("log10", function(x) 10^x),
                    labels = scales::trans_format("log10", scales::math_format(10^.x)),
                    sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
      xlab(label_a_phi_SOM)+
      ylab(label_decomposition_rate)+
      ggtitle("FOM")
  
p2<-ggplot(data=data[data$type=="faeces",])+
      geom_rect(aes(xmin=as.numeric(a_phi_SOM)-0.45,xmax=as.numeric(a_phi_SOM)+0.45,ymin=min(decomposition_rate,na.rm=TRUE)*0.5,ymax=decomposition_rate))+
      facet_grid(a_phi_DOC~a_phi_FOM)+
      theme+theme(legend.position="none",
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
      scale_x_continuous(breaks=1:length(levels(data$a_phi_SOM)),
                         labels=levels(data$a_phi_SOM),
                         sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
      scale_y_continuous(sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
      xlab(label_a_phi_SOM)+
      ylab(label_decomposition_rate)+
      ggtitle("faeces")

p3<-ggplot(data=data[data$type=="SOM",])+
      geom_rect(aes(xmin=as.numeric(a_phi_SOM)-0.45,xmax=as.numeric(a_phi_SOM)+0.45,ymin=min(decomposition_rate,na.rm=TRUE)*0.5,ymax=decomposition_rate))+
      facet_grid(a_phi_DOC~a_phi_FOM)+
      theme+theme(legend.position="none",
                  panel.grid.minor.x=element_blank(),
                  panel.grid.major.x=element_blank())+
      palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
      scale_x_continuous(breaks=1:length(levels(data$a_phi_SOM)),
                         labels=levels(data$a_phi_SOM),
                         sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
      scale_y_continuous(sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
      xlab(label_a_phi_SOM)+
      ylab(label_decomposition_rate)+
      ggtitle("SOM")

data$weighted_decomposition<-data$decomposition_rate*data$stock
data<-data[data$type%in%c("FOM","faeces"),]
data<-aggregate(cbind(weighted_decomposition,stock)~a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data, sum)
data$weighted_decomposition<-data$weighted_decomposition/data$stock

ggplot(data=data)+
  geom_rect(aes(xmin=as.numeric(a_phi_SOM)-0.45,xmax=as.numeric(a_phi_SOM)+0.45,ymin=min(weighted_decomposition,na.rm=TRUE)*0.5,ymax=weighted_decomposition))+
  facet_grid(a_phi_DOC~a_phi_FOM)+
  theme+theme(legend.position="none",
              panel.grid.major.x=element_blank())+
  palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
  scale_x_continuous(breaks=1:length(levels(data$a_phi_SOM)),
                     labels=levels(data$a_phi_SOM),
                     sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
  scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                labels = scales::trans_format("log10", scales::math_format(10^.x)),
                sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
  xlab(label_a_phi_SOM)+
  ylab(label_decomposition_rate)+
  ggtitle("FOM + faeces")
  
graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1, 1)+
  draw_plot(p3, 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_decomposition_rate.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

### respiration #----
data_species$respiration[data_species$respiration<1e-7]=NA
data<-aggregate(respiration~trophic_group+a_phi_DOC+a_phi_FOM+a_phi_SOM, data=data_species, sum)
a_phi_SOM=unique(data$a_phi_SOM)

p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data)+
    geom_bar(aes(1,respiration,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
    facet_grid(a_phi_DOC~a_phi_FOM)+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_decomposition$scale_colour_bar_trophic_groups+
    scale_x_continuous(sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
    scale_y_continuous(breaks=c(0,0.5,1), labels = scales::percent,
                       sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
    ylab(label_respiration_fraction)+
    ggtitle(parse(text=paste0(label_a_phi_SOM,'*" = "*',a_phi_SOM[i])))
}

graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_respiration.pdf",sep=""), graph, width = 12, height = 5, device=cairo_pdf)

### equilibrium #----
names(data_TS)[2:(ncol(data_TS)-3)]<-data_sensitivity_decomposition$names_species
data_TS<-Rmake_table_for_plot(data_TS,c("time","a_phi_DOC","a_phi_FOM","a_phi_SOM"),"names","biomass")
data_TS$biomass[data_TS$biomass<1e-7]=NA

time_range<-max(data_TS$time)-min(data_TS$time)
time_range<-floor(min(data_TS$time)+time_range*c(0.25,0.75)) # first and last ticks on the time axis
close_order<-floor(log10(time_range))-1 # to round the ticks to hundreds
time_range<-round(time_range*10^(-close_order))*10^close_order
time_range<-time_range/100 # shorter breaks

biomass_range<-floor(log10(c(min(data_TS$biomass,na.rm=T),max(data_TS$biomass,na.rm=T))))
biomass_range<-10^seq(biomass_range[1],biomass_range[2],3)

p<-as.list(array(0,length(a_phi_SOM)))
for (i in 1:length(a_phi_SOM)){
  p[[i]]<-ggplot(data=data_TS[data_TS$a_phi_SOM==a_phi_SOM[i],])+
            geom_line(aes(time/100,biomass,colour=names),size=1.5)+
            facet_grid(a_phi_DOC~a_phi_FOM)+
            theme+theme(legend.position="none",
                        panel.spacing.x = unit(1, "lines"))+
            palette_sensitivity_decomposition$scale_colour_line_species+
            scale_x_continuous(breaks = time_range, name = label_time_short,
                               sec.axis = sec_axis(~ . , name = label_a_phi_FOM, breaks = NULL, labels = NULL))+
            scale_y_log10(breaks = biomass_range,#trans_breaks("log10", function(x) 10^x),
                          labels = scales::trans_format("log10", scales::math_format(10^.x)),
                          name = label_biomass,
                          sec.axis = sec_axis(~ . , name = label_a_phi_DOC, breaks = NULL, labels = NULL))+
            ggtitle(parse(text=paste0(label_a_phi_SOM_short,'*" = "*',a_phi_SOM[i])))
}
graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_decomposition_equilibrium.pdf",sep=""), graph, width = 12, height = 5, device=cairo_pdf)

# SENSITIVITY ANALYSIS HALF-SATURATION # ----
### load data # ----
path_results<-"results_final_half_saturation/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_sensitivity_half_saturation<-Rload_data(path_results,1)
palette_sensitivity_half_saturation<-Rmake_palettes(data_sensitivity_half_saturation$species_traits,
                                                    data_sensitivity_half_saturation$n_tot_species,
                                                    data_sensitivity_half_saturation$n_faeces)
data_species<-NULL
data_detritus<-NULL
data_TS<-NULL
for (i in 1:nrow(params)){
  # species data
  data_sensitivity_half_saturation<-Rload_data(path_results,i)
  temp<-data_sensitivity_half_saturation$data_species
  data_species<-rbind(data_species,temp)
  # detritus data
  temp<-data_sensitivity_half_saturation$data_detritus
  data_detritus<-rbind(data_detritus,temp)
  # density time series data
  temp<-read.table(paste(path_results,"TS_density_",i,".txt",sep=""),sep=",",header=T)
  temp<-temp[,1:(data_sensitivity_half_saturation$n_tot_species+1)]
  temp$a_K_DOC=params$a_K_DOC[i]
  temp$a_K_FOM=params$a_K_FOM[i]
  temp$a_K_SOM=params$a_K_SOM[i]
  data_TS<-rbind(data_TS,temp)
}
rm(temp)

### biomass distribution # ----
data_species$biomass[data_species$biomass<1e-7]=NA
data<-aggregate(biomass~trophic_group+a_K_DOC+a_K_FOM+a_K_SOM, data=data_species, sum)
a_K_SOM=unique(data$a_K_SOM)
biomass_range<-floor(log10(c(min(data$biomass,na.rm=T),max(data$biomass,na.rm=T)))) # rescale the biomass axis
biomass_range<-10^c(2,3)#rev(seq(biomass_range[2],biomass_range[1],-2))

# biomass distribution
p<-as.list(array(0,length(a_K_SOM)))
for (i in 1:length(a_K_SOM)){
  p[[i]]<-ggplot(data)+
    geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
    facet_grid(a_K_DOC~a_K_FOM)+
    theme+theme(legend.position="none",
                axis.title.y.left=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_blank())+
    palette_sensitivity_half_saturation$scale_colour_bar_trophic_groups+
    scale_x_continuous(labels=levels(data$trophic_group), breaks=c(1:length(levels(data$trophic_group))),
                       sec.axis = sec_axis(~ . , name = label_a_K_DOC, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = biomass_range,#trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  sec.axis = sec_axis(~ . , name = label_a_K_FOM, breaks = NULL, labels = NULL))+
    coord_flip()+
    ylab(label_biomass)+
    ggtitle(parse(text=paste0(label_a_K_SOM,'*" = "*',a_K_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_half_saturation_biomass.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

# update of the parameter space to remove combination without coexistence
#data<-data[-which(data$a_K_FOM==10 | (data$a_K_FOM==100 & data$a_K_DOC==1) | data$a_K_SOM==100),]
#data_detritus<-data_detritus[-which(data_detritus$a_K_FOM==10 | (data_detritus$a_K_FOM==100 & data_detritus$a_K_DOC==1)),]
a_K_SOM=unique(data$a_K_SOM)

# relative biomass distribution
p<-as.list(array(0,length(a_K_SOM)))
for (i in 1:length(a_K_SOM)){
  p[[i]]<-ggplot(data=data[data$a_K_SOM==a_K_SOM[i],])+
    geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
    facet_grid(a_K_DOC~a_K_FOM)+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_half_saturation$scale_colour_bar_trophic_groups+
    scale_x_continuous(sec.axis = sec_axis(~ . , name = label_a_K_FOM, breaks = NULL, labels = NULL))+
    scale_y_continuous(breaks=c(0,0.5,1), labels = scales::percent,
                       sec.axis = sec_axis(~ . , name = label_a_K_DOC, breaks = NULL, labels = NULL))+
    ylab(label_biomass_fraction)+
    ggtitle(parse(text=paste0(label_a_K_SOM,'*" = "*',a_K_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_half_saturation_biomass_relative.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

### distribution stock and decomposition # ----
# detritus stocks
data<-aggregate(biomass~type+a_K_DOC+a_K_FOM+a_K_SOM, data=data_detritus, sum)
biomass_range<-floor(log10(c(min(data$biomass[data$biomass>0],na.rm=T),max(data$biomass,na.rm=T)))) # rescale the biomass axis
biomass_range<-10^rev(seq(biomass_range[2],biomass_range[1],-2))
a_K_SOM=unique(data$a_K_SOM)

p<-as.list(array(0,length(a_K_SOM)))
for (i in 1:length(a_K_SOM)){
  p[[i]]<-ggplot(data=data[data$a_K_SOM==a_K_SOM[i],])+
    geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(biomass)*0.5,ymax=biomass,fill=type))+
    facet_grid(a_K_DOC~a_K_FOM)+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_half_saturation$scale_colour_bar_detritus_classes+
    scale_x_continuous(breaks=1:length(levels(data$type)),
                       labels=levels(data$type),
                       sec.axis = sec_axis(~ . , name = label_a_K_FOM, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = biomass_range,#trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  sec.axis = sec_axis(~ . , name = label_a_K_DOC, breaks = NULL, labels = NULL))+
    ylab(label_stock)+
    ggtitle(parse(text=paste0(label_a_K_SOM,'*" = "*',a_K_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_half_saturation_detritus.pdf",sep=""), graph, width = 18, height = 5, device=cairo_pdf)

# update of the parameter space to remove combination without coexistence
#data_species<-data_species[-which(data_species$a_K_FOM==10 | (data_species$a_K_FOM==100 & data_species$a_K_DOC==1) | data_species$a_K_SOM==100),]
#data_detritus<-data_detritus[-which(data_detritus$a_K_SOM==100),]
#data_TS<-data_TS[-which(data_TS$a_K_FOM==10 | (data_TS$a_K_FOM==100 & data_TS$a_K_DOC==1) | data_TS$a_K_SOM==100),]
#a_K_SOM=unique(data_species$a_K_SOM)

# decomposition by each decomposer trophic group
detritus_order<-c("FOM","faeces","SOM","DOC")
consumption_temp=NULL
for (i in 1:nrow(params)){
  consumption<-read.table(paste(path_results,"consumption_",params$simu_ID[i],".txt",sep=""),sep=",",header=T)
  names(consumption)<-data_sensitivity_half_saturation$names_tot
  consumption<-consumption[,data_sensitivity_half_saturation$names_detritus] # selects abiotic resources only
  consumption$consumer<-data_sensitivity_half_saturation$names_species
  consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
  # resource type
  temp<-data_sensitivity_half_saturation$detritus_traits[,c("type","names")]
  names(temp)<-c("type","resource")
  consumption<-merge(consumption,temp,by="resource")
  # consumer trophic group
  temp<-data_sensitivity_half_saturation$species_traits[,c("trophic_group","names")]
  names(temp)<-c("trophic_group","consumer")
  consumption<-merge(consumption,temp,by="consumer")
  rm(temp)
  # only consider decomposers
  consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores"),]
  consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
  consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores"))
  consumption$type<-factor(consumption$type,levels=detritus_order)
  consumption$a_K_DOC=params$a_K_DOC[i]
  consumption$a_K_FOM=params$a_K_FOM[i]
  consumption$a_K_SOM=params$a_K_SOM[i]
  consumption_temp<-rbind(consumption_temp,consumption)
}
consumption<-consumption_temp
rm(consumption_temp)
#consumption<-consumption[-which(consumption$a_K_FOM==10 | (consumption$a_K_FOM==100 & consumption$a_K_DOC==1) | consumption$a_K_SOM==100),]

p<-as.list(array(0,length(a_K_SOM)))
for (i in 1:length(a_K_SOM)){
  p[[i]]<-ggplot(data=consumption[consumption$a_K_SOM==a_K_SOM[i],])+
    geom_bar(aes(as.numeric(type),biomass,fill=trophic_group),stat="identity")+
    facet_grid(a_K_DOC~a_K_FOM)+
    theme+theme(legend.position="bottom",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_text(angle=90,vjust=0.5),
                panel.grid.major.x=element_blank())+
    palette_sensitivity_decomposition$scale_colour_bar_trophic_groups_detritivores+
    scale_x_continuous(breaks=1:length(levels(consumption$type)),
                       labels=levels(consumption$type),
                       sec.axis = sec_axis(~ . , name = label_a_K_FOM, breaks = NULL, labels = NULL))+
    scale_y_continuous(sec.axis = sec_axis(~ . , name = label_a_K_DOC, breaks = NULL, labels = NULL))+
    ylab(label_decomposition)+
    ggtitle(parse(text=paste0(label_a_K_SOM,'*" = "*',a_K_SOM[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_sensitivity_half_saturation_decomposition.pdf",sep=""), graph, width = 12, height = 6, device=cairo_pdf)

# SENSITIVITY ANALYSIS SELF-REGULATION PLUS ALLOMETRIC EXPONENT # ----
### load data # ----
path_results<-"results_final_self_regulation/"#"results_23_02_20/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)
data_biomass<-NULL
data_TS<-NULL
for (i in 1:nrow(params)){
  # biomass
  temp<-Rload_data(path_results,i)
  temp<-temp$data_species[,c("c0_prey","c0_pred","s_c","trophic_group","names","biomass")]
  data_biomass<-rbind(data_biomass,temp)
  # times series
  temp<-read.table(paste(path_results,"TS_density_",i,".txt",sep=""),sep=",",header=T)
  temp<-temp[,1:(data_multi_channel$n_tot_species+1)]
  temp$c0_prey=params$c0_prey[i]
  temp$c0_pred=params$c0_pred[i]
  temp$s_c=params$s_c[i]
  data_TS<-rbind(data_TS,temp)
}

### Biomass # ----
#data_biomass<-data_biomass[data_biomass$c0_prey>5e3 & data_biomass$c0_pred>5e4,]
data_biomass$c0_pred_label<-as.factor(log10(data_biomass$c0_pred))
for (i in 1:length(levels(data_biomass$c0_pred_label))){
  levels(data_biomass$c0_pred_label)[i]=paste("10^",levels(data_biomass$c0_pred_label)[i])
}
data_biomass$c0_prey_label<-as.factor(log10(data_biomass$c0_prey))
for (i in 1:length(levels(data_biomass$c0_prey_label))){
  levels(data_biomass$c0_prey_label)[i]=paste("10^",levels(data_biomass$c0_prey_label)[i])
}
data_biomass$biomass[data_biomass$biomass<1e-7]=NA
data_biomass$s_c_label<-as.factor(data_biomass$s_c)
s_c=levels(data_biomass$s_c_label)

# species biomass
p<-as.list(array(0,length(s_c)))
for (i in 1:length(s_c)){
  p[[i]]<-ggplot(data=data_biomass[data_biomass$s_c==s_c[i],])+
    geom_rect(aes(xmin=as.numeric(names)-0.5,xmax=as.numeric(names)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=names))+
    facet_grid(c0_prey_label~c0_pred_label, labeller=label_parsed, scale="free")+
    theme+theme(legend.position="none",
                axis.title.y.left=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_blank())+
    palette_multi_channel$scale_colour_bar_species+
    coord_flip()+
    scale_x_continuous(sec.axis = sec_axis(~ . , name = label_c0_prey, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = 10^seq(-4,4,2),#trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  name = label_biomass,
                  sec.axis = sec_axis(~ . , name = label_c0_pred, breaks = NULL, labels = NULL))+
    ggtitle(parse(text=paste0('italic(s[c])*" = "*',s_c[i])))
}

graph<-ggdraw(xlim = c(0,3), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)+
  draw_plot(p[[3]], 2, 0, 1, 1)
ggsave(paste(path_figures,"supp_interference_biomass_species.pdf",sep=""), graph, width = 20, height = 8, device=cairo_pdf)

# update the parameter space
data_biomass<-data_biomass[data_biomass$s_c>-0.4,]
data_TS<-data_TS[data_TS$s_c>-0.4,]
data_biomass$s_c_label<-droplevels(data_biomass$s_c_label)
s_c=levels(data_biomass$s_c_label)

# raw biomass per trophic group
data<-aggregate(biomass~trophic_group+c0_prey_label+c0_pred_label+s_c_label, data=data_biomass, sum)
p<-as.list(array(0,length(s_c)))
for (i in 1:length(s_c)){
  p[[i]]<-ggplot(data=data[data$s_c==s_c[i],])+
    geom_rect(aes(xmin=as.numeric(trophic_group)-0.5,xmax=as.numeric(trophic_group)+0.5,ymin=min(biomass,na.rm=TRUE)*0.5,ymax=biomass,fill=trophic_group))+
    facet_grid(c0_prey_label~c0_pred_label, labeller=label_parsed, scale="free")+
    theme+theme(legend.position="none",
                axis.title.y.left=element_blank(),
                axis.text.y=element_blank(),
                axis.ticks.y=element_blank(),
                panel.grid.minor.y=element_blank(),
                panel.grid.major.y=element_blank())+
    palette_multi_channel$scale_colour_bar_trophic_groups+
    coord_flip()+
    scale_x_continuous(labels=levels(data$trophic_group), breaks=c(1:length(levels(data$trophic_group))),
                       sec.axis = sec_axis(~ . , name = label_c0_prey, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  name = label_biomass,
                  sec.axis = sec_axis(~ . , name = label_c0_pred, breaks = NULL, labels = NULL))+
    ggtitle(parse(text=paste0('italic(s[c])*" = "*',s_c[i])))
}

graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_interference_biomass.pdf",sep=""), graph, width = 20, height = 11, device=cairo_pdf)

# relative biomass distribution
temp<-aggregate(biomass~c0_prey_label+c0_pred_label+s_c_label, data=data, sum)
names(temp)<-c("c0_prey_label","c0_pred_label","s_c_label","relative")
data<-merge(data,temp,by=c("c0_prey_label","c0_pred_label","s_c_label"))
data$relative<-data$biomass/data$relative

p<-as.list(array(0,length(s_c)))
for (i in 1:length(s_c)){
  p[[i]]<-ggplot(data=data[data$s_c_label==s_c[i],])+
    geom_bar(aes(1,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
    facet_grid(c0_prey_label~c0_pred_label, labeller=label_parsed, scale="free")+
    theme+theme(legend.position="none",
                axis.title.x.bottom=element_blank(),
                axis.text.x=element_blank(),
                axis.ticks.x=element_blank(),
                panel.grid.minor.x=element_blank(),
                panel.grid.major.x=element_blank())+
    palette_multi_channel$scale_colour_bar_trophic_groups+
    scale_x_continuous(sec.axis = sec_axis(~ . , name = label_c0_pred, breaks = NULL, labels = NULL))+
    scale_y_continuous(name=label_biomass_fraction, labels = scales::percent,
                       sec.axis = sec_axis(~ . , name = label_c0_prey, breaks = NULL, labels = NULL))+
    ggtitle(parse(text=paste0('italic(s[c])*" = "*',s_c[i])))
}

graph<-ggdraw(xlim = c(0,2), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)+
  draw_plot(p[[2]], 1, 0, 1, 1)
ggsave(paste(path_figures,"supp_interference_fraction.pdf",sep=""), graph, width = 20, height = 8, device=cairo_pdf)

### Time series # ----
# update the parameter space
data_TS<-data_TS[data_TS$s_c==-0.3,]
# formate data
names(data_TS)[2:(ncol(data_TS)-3)]<-data_multi_channel$names_species
#data_TS<-data_TS[data_TS$c0_prey>5e3 & data_TS$c0_pred>5e4,]
data_TS<-Rmake_table_for_plot(data_TS,c("time","c0_prey","c0_pred","s_c"),"names","biomass")
data_TS$c0_pred<-as.factor(log10(data_TS$c0_pred))
for (i in 1:length(levels(data_TS$c0_pred))){
  levels(data_TS$c0_pred)[i]=paste("10^",levels(data_TS$c0_pred)[i])
}
data_TS$c0_prey<-as.factor(log10(data_TS$c0_prey))
for (i in 1:length(levels(data_TS$c0_prey))){
  levels(data_TS$c0_prey)[i]=paste("10^",levels(data_TS$c0_prey)[i])
}
# set the breaks for the time axis
time_range<-max(data_TS$time)-min(data_TS$time)
time_range<-floor(min(data_TS$time)+time_range*c(0.25,0.75))
close_order<-floor(log10(time_range))-1
time_range<-round(time_range*10^(-close_order))*10^close_order
time_range<-time_range/100 # shorter breaks
data_TS$biomass[data_TS$biomass<1e-6]=NA
data_TS$s_c<-as.factor(data_TS$s_c)

p<-as.list(array(0,length(s_c)))
for (i in 1:length(s_c)){
  p[[i]]<-ggplot(data=data_TS[data_TS$s_c==s_c[i],])+
    geom_line(aes(time/100,biomass,colour=names),linewidth=1.5)+
    facet_grid(c0_prey~c0_pred, labeller=label_parsed, scale="free")+
    theme+theme(legend.position="none",
                panel.spacing.x = unit(1, "lines"))+
    palette_multi_channel$scale_colour_line_species+
    scale_x_continuous(breaks = time_range,
                       name = label_time_short,
                       sec.axis = sec_axis(~ . , name = label_c0_pred, breaks = NULL, labels = NULL))+
    scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                  labels = scales::trans_format("log10", scales::math_format(10^.x)),
                  name = label_biomass,
                  sec.axis = sec_axis(~ . , name = label_c0_prey, breaks = NULL, labels = NULL))+
    ggtitle(parse(text=paste0('italic(s[c])*" = "*',s_c[i])))
}

graph<-ggdraw(xlim = c(0,1), ylim = c(0,1)) +
  draw_plot(p[[1]], 0, 0, 1, 1)
ggsave(paste(path_figures,"supp_interference_TS.pdf",sep=""), graph, width = 20, height = 8, device=cairo_pdf)

# SENSITIVITY ANALYSIS DORMANCY ASYMMETRY # ----
path_results<-"results_final_dormancy/"#"results_23_02_21_2/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_multi_channel<-Rload_data(path_results,1)
palette_multi_channel<-Rmake_palettes(data_multi_channel$species_traits,
                                      data_multi_channel$n_tot_species,
                                      data_multi_channel$n_faeces)

# Fraction of active microbes
data_active<-NULL
for (i in 1:nrow(params)){
  temp<-Rload_data(path_results,i)
  temp<-temp$data_species[temp$data_species$trophic_group=="microbes",c("Q_asy","names","active")]
  data_active<-rbind(data_active,temp)
}
data<-aggregate(active~+Q_asy, data=data_active, mean)

p1<-ggplot(data)+
      geom_line(aes(Q_asy,active),linewidth=2)+
      theme+
      x_axis_log10+
      scale_y_continuous(labels = scales::percent)+
      xlab(label_Q_asy)+
      ylab(label_active)

# Relative biomass
data_biomass<-NULL
for (i in 1:nrow(params)){
  temp<-Rload_data(path_results,i)
  temp<-temp$data_species[,c("Q_asy","trophic_group","names","biomass")]
  data_biomass<-rbind(data_biomass,temp)
}
data_biomass$biomass[data_biomass$biomass<1e-7]=NA
data<-aggregate(biomass~trophic_group+Q_asy, data=data_biomass, sum)
temp<-aggregate(biomass~+Q_asy, data=data, sum)
names(temp)<-c("Q_asy","relative")
data<-merge(data,temp,by=c("Q_asy"))
data$relative<-data$biomass/data$relative

p2<-ggplot(data)+
  geom_area(aes(Q_asy,relative,fill=trophic_group),position = position_stack(reverse = TRUE))+
  theme+theme(legend.position="bottom",
              legend.direction="vertical")+
  guides(fill=guide_legend(nrow=3,reverse=TRUE))+
  palette_multi_channel$scale_colour_bar_trophic_groups
legend<-get_legend(p2)

p2<-ggplot(data)+
      geom_area(aes(Q_asy,relative,fill=trophic_group),position = position_stack(reverse = TRUE))+
      theme+theme(legend.position="none")+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      x_axis_log10+
      scale_y_continuous(labels = scales::percent)+
      xlab(label_Q_asy)+
      ylab(label_biomass_fraction)

# Relative respiration
data_respiration<-NULL
for (i in 1:nrow(params)){
  temp<-Rload_data(path_results,i)
  temp<-temp$data_species[,c("Q_asy","trophic_group","names","biomass","respiration")]
  data_respiration<-rbind(data_respiration,temp)
}
data_respiration$respiration[data_respiration$biomass<1e-7]=NA
data<-aggregate(respiration~trophic_group+Q_asy, data=data_respiration, sum)
temp<-aggregate(respiration~+Q_asy, data=data, sum)
names(temp)<-c("Q_asy","relative")
data<-merge(data,temp,by=c("Q_asy"))
data$relative<-data$respiration/data$relative

p3<-ggplot(data)+
      geom_area(aes(Q_asy,relative,fill=trophic_group),position = position_stack(reverse = TRUE))+
      theme+theme(legend.position="none")+
      palette_multi_channel$scale_colour_bar_trophic_groups+
      x_axis_log10+
      scale_y_continuous(labels = scales::percent)+
      xlab(label_Q_asy)+
      ylab(label_respiration_fraction)

h=0.2
graph<-ggdraw(xlim = c(0, 3), ylim = c(0, 1.1+h)) +
  draw_plot(p1, 0, h, 1, 1)+
  draw_plot(p2, 1, h, 1, 1)+
  draw_plot(p3, 2, h, 1, 1)+
  draw_plot(legend, 1.5, 0, 1.5, h)+
  draw_plot_label(c("A","B","C"), c(0,1,2), c(1,1,1)+h+0.1, size = 30)
ggsave(paste(path_figures,"supp_dormancy.pdf",sep=""), graph, width = 15, height = 5.8, device=cairo_pdf)

# EMPIRIC DATA # ----
# Brose et al. 2014
data_abundance<-read.table("Data/Brose 2014 special issue/17471_2_Dataset/17471_2_data.csv",header=T,sep=",")
data<-data_abundance[,c(2,12:ncol(data_abundance))]
temp1<-data[,c(1,2:(ncol(data)/2+1))]
temp2<-data[,c(1,(ncol(temp1)+1):ncol(data))]
names(temp2)<-names(temp1)
temp1<-melt(temp1, id.vars = c("Number"), # density per square meter
            variable.name = "plot", 
            value.name = "abundance")
temp2<-melt(temp2, id.vars = c("Number"),
            variable.name = "plot", 
            value.name = "bodymass") # body mass (mg fresh weight)
data<-data_abundance[,1:11]
data<-merge(data,temp1,by=c("Number"))
data<-merge(data,temp2,by=c("Number","plot"))
rm(temp1,temp2)
fresh_to_dry=0.2 # Makarieva et al. (2005), Ehnes et al. (2011) Johnston and Sibly (2018)
dry_to_C=0.42 # Andrieux et al. (2021) (for animals and litter)
data$biomass<-data$abundance*data$bodymass*fresh_to_dry*dry_to_C
data_abundance<-data

data_microbes<-read.table("Data/Brose 2014 special issue/Xu_2013_Global_Ecology_and_Biogeography_germany.csv",header=T,sep=",")
data_microbes<-mean(data_microbes[,4]) #mmolC kg-1 soil (dw)
molC_to_mgC=12.011e3 # conversion to mg !!!!!
volume_to_surface=0.1 # 1x1x0.1 m parallelepiped
soil_mass_to_volume= 1/1.4e3 #1/1e3 # kg m-3 Bowden et al. (2014) or 1.4 g.cm-3 Heathman et al. (2003)
data_microbes<-data_microbes*1e-3*molC_to_mgC/(soil_mass_to_volume*volume_to_surface)

data<-aggregate(biomass~plot+feeding_type_general, data=data_abundance, sum)
#data<-separate(data, funct_group, c("size_class", NA), sep = "_")
data<-aggregate(biomass~feeding_type_general, data=data, mean)
data<-rbind(data,c("microbes",data_microbes))
data$feeding_type_general<-as.factor(data$feeding_type_general)
data$biomass<-as.numeric(data$biomass)

ggplot(data=data)+
  geom_bar(aes(feeding_type_general,biomass,fill=feeding_type_general),stat="identity")+
  theme+
  coord_flip()

# SENSITIVITY ANALYSIS RADIUS FOM # ----
path_results<-"results_final_radius_FOM/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_sensitivity_radius_FOM<-Rload_data(path_results,1)
palette_sensitivity_radius_FOM<-Rmake_palettes(data_sensitivity_radius_FOM$species_traits,
                                               data_sensitivity_radius_FOM$n_tot_species,
                                               data_sensitivity_radius_FOM$n_faeces)

data<-aggregate(biomass~type+radius_FOM, data=data_sensitivity_radius_FOM$data_detritus, sum)
ggplot(data)+
  geom_rect(aes(xmin=as.numeric(type)-0.45,xmax=as.numeric(type)+0.45,ymin=min(biomass)*0.5,ymax=biomass,fill=type))+
  facet_wrap(~radius_FOM)+
  theme_vbar+theme(legend.position="none")+
  palette_sensitivity_radius_FOM$scale_colour_bar_detritus_classes+
  scale_x_continuous(breaks=1:length(levels(data$type)),
                     labels=levels(data$type))+
  y_axis_log10+
  ylab(label_stock)

# decomposition by each decomposer trophic group
data<-NULL
for (i in 1:nrow(params)){
  consumption<-read.table(paste(path_results,"consumption_",params$simu_ID[i],".txt",sep=""),sep=",",header=T)
  names(consumption)<-data_sensitivity_radius_FOM$names_tot
  consumption<-consumption[,data_sensitivity_radius_FOM$names_detritus] # selects abiotic resources only
  consumption$consumer<-data_sensitivity_radius_FOM$names_species
  consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
  # resource type
  temp<-data_sensitivity_radius_FOM$detritus_traits[,c("type","names")]
  names(temp)<-c("type","resource")
  consumption<-merge(consumption,temp,by="resource")
  # consumer trophic group
  temp<-data_sensitivity_radius_FOM$species_traits[,c("trophic_group","names")]
  names(temp)<-c("trophic_group","consumer")
  consumption<-merge(consumption,temp,by="consumer")
  rm(temp)
  # only consider decomposers
  consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores","trophic whales"),]
  consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
  consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores","trophic whales"))
  levels(consumption$trophic_group)<-c("microbes","macro-detritivores","trophic whales")
  consumption$type<-factor(consumption$type,levels=c("faeces","FOM","SOM","DOC"))
  consumption$radius_FOM<-params$radius_FOM[i]
  data<-rbind(data,consumption)
}

ggplot(data=data)+
  geom_bar(aes(type,biomass,fill=trophic_group),stat="identity",position=position_stack(reverse = TRUE))+
  facet_wrap(~radius_FOM)+
  theme_vbar+theme(#legend.position="top",
    legend.position=c(1, 1),
    legend.justification=c(1, 1))+
  palette_sensitivity_radius_FOM$scale_colour_bar_trophic_groups_detritivores+
  ylab(label_decomposition)

# SENSITIVITY ANALYSIS PREF SIGMA # ----
path_results<-"results_final_pref_sigma/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_sensitivity_pref_sigma<-Rload_data(path_results,1)
palette_sensitivity_pref_sigma<-Rmake_palettes(data_sensitivity_pref_sigma$species_traits,
                                               data_sensitivity_pref_sigma$n_tot_species,
                                               data_sensitivity_pref_sigma$n_faeces)
data_species<-NULL
for (i in 1:nrow(params)){
  # species data
  data_sensitivity_pref_sigma<-Rload_data(path_results,i)
  temp<-data_sensitivity_pref_sigma$data_species
  data_species<-rbind(data_species,temp)
}
rm(temp)

# biomass
data<-data_species
data$biomass[data$biomass<1e-7]=NA
data<-aggregate(biomass~trophic_group+pref_sigma, data=data, sum)
levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
data$relative<-0
for (i in unique(data$pref_sigma)){
  data$relative[data$pref_sigma==i]<-data$biomass[data$pref_sigma==i]/sum(data$biomass[data$pref_sigma==i])
}
data$relative<-data$relative*100
data$pref_sigma<-as.factor(data$pref_sigma)

p1<-ggplot(data=data)+
  geom_bar(aes(pref_sigma,biomass,fill=trophic_group),stat="identity",position=position_fill(reverse = TRUE))+
  theme+theme(legend.position="top",
              panel.grid.minor.x=element_blank(),
              panel.grid.major.x=element_blank())+
  palette_sensitivity_pref_sigma$scale_colour_bar_trophic_groups+
  guides(fill=guide_legend(nrow=3,reverse=TRUE))+
  xlab(expression(italic("\u03C3"^M)))+
  scale_y_continuous(labels = scales::percent)+
  ylab(label_biomass_fraction)

# TL
data<-data_species
data$TL[data$TL==0]=NA

data_potapov<-read.table(paste(path_potapov,"RawData_MassDensityTrophicLevel.csv",sep=""),sep=",",header=T)
data_potapov<-data_potapov[,c("SizeClass","TL")]
data_potapov$TL<-data_potapov$TL-1
data_potapov<-data_potapov %>% separate(SizeClass, c("bodymass",NA), sep = "-")
data_potapov$bodymass[data_potapov$bodymass==">1000"]=1000
data_potapov$bodymass<-as.numeric(data_potapov$bodymass)

data$pref_sigma<-as.factor(data$pref_sigma)
levels(data$pref_sigma)<-paste0(expression(italic("\u03C3"^M)),'*" = "*',levels(data$pref_sigma))

p2<-ggplot(data=data)+
  geom_point(aes(bodymass,TL,colour=names),size=4)+
  geom_smooth(data=data_potapov, aes(bodymass,TL), se=FALSE, colour="black", linetype="dashed")+
  geom_smooth(aes(bodymass,TL), colour="black", linetype="solid")+
  facet_wrap(~pref_sigma,labeller=label_parsed)+
  theme+theme(legend.position="none")+
  palette_sensitivity_pref_sigma$scale_colour_line_species+
  scale_x_log10(minor_breaks = 10^seq(-10,5,1),
                breaks = 10^(seq(-8,2,2)),
                labels = scales::trans_format("log10", scales::math_format(10^.x)))+
  xlab(label_bodymass)+
  ylab(label_TL)

# diet
data_diet<-NULL
for (i in 1:nrow(params)){
  data_sensitivity_pref_sigma<-Rload_data(path_results,i)
  data<-data_sensitivity_pref_sigma$data_species
  diet<-read.table(paste(path_results,"diet_",data$simu_ID[i],".txt",sep=""),sep=",",header=T)
  diet$pred<-data_multi_channel$names_species # add the names of the predators
  diet$pred<-factor(diet$pred,levels=data_multi_channel$names_species)
  data<-data[order(data$names),]
  diet[data$TL==0,c(1:(ncol(diet)-1))]=NA # removes the diet of extinct carnivores
  diet<-diet[which(data$trophic_type%in%c("microbivores","carnivores")),c(1:data_multi_channel$n_tot_species,ncol(diet))] # keeps only consumers
  diet<-Rmake_table_for_plot(diet,"pred","prey","flow") # makes diet a long table
  levels(diet$prey)<-data_multi_channel$names_species
  survivor<-data$names[data$TL>0 & data$trophic_type%in%c("microbivores","carnivores")]
  diet<-diet[diet$pred%in%survivor,]
  diet$pref_sigma<-params$pref_sigma[i]
  data_diet<-rbind(data_diet,diet)
}
data_diet$pref_sigma<-as.factor(data_diet$pref_sigma)
levels(data_diet$pref_sigma)<-paste0(expression(italic("\u03C3"^M)),'*" = "*',levels(data_diet$pref_sigma))

palette_survivors<-Rmake_colour_palette(data_sensitivity_pref_sigma$species_traits,0)[which(data_sensitivity_pref_sigma$names_species%in%survivor)]
scale_colour_species_survivors<-scale_colour_manual(values=palette_survivors,
                                                    name=NULL)

p3<-ggplot(data=data_diet)+
      geom_bar(aes(pred,flow,fill=prey),stat="identity",position="fill")+
      geom_point(aes(pred,-0.05,colour=pred),size=4)+
      facet_wrap(~pref_sigma,labeller=label_parsed)+
      theme+theme(legend.position="none",
                  axis.text.y=element_blank(),
                  axis.ticks.y=element_blank(),
                  axis.line.y = element_line(arrow = arrow(angle = 15, length = unit(.15,"inches"),type = "closed")),
                  #axis.line.y=element_blank(),
                  panel.grid=element_blank())+
      scale_colour_species_survivors+
      palette_multi_channel$scale_colour_bar_species+
      scale_y_continuous(breaks=seq(0,1,0.25),labels = scales::percent)+
      coord_flip()+
      xlab("Increasing body size of predators")+
      ylab(label_diet)

levels(data$trophic_group)<-c("microbes","microbivores","micro-carnivores","macro-detritivores","macro-carnivores","trophic whales")
legend<-ggplot(data=data)+
  geom_bar(aes(trophic_group,fill=names),stat="count")+
  theme_void()+theme(legend.position="none",
                     axis.text.y=element_text(family="serif",size=15,hjust=1))+
  palette_multi_channel$scale_colour_bar_species+
  coord_flip()+
  ylab(label_biomass_fraction)

graph<-ggdraw(xlim = c(0, 2), ylim = c(0, 2)) +
  draw_plot(p1, 0, 1, 0.6, 1)+
  draw_plot(p2, 0.6, 1, 1.4, 1)+
  draw_plot(p3, 0, 0, 1.5, 1)+
  draw_plot(legend, 1.5, 0.25, 0.5, 0.5)+
  draw_plot_label(c("A","B","C"), c(0,0.59,0), c(2,2,1.05), size = 30)
ggsave(paste(path_figures,"supp_sensitivity_pref_sigma_biomass_TL.pdf",sep=""), graph, width = 17, height = 10, device=cairo_pdf)

# SENSITIVITY ANALYSIS RADIUS FOM AND DECOMPOSITION # ----
path_results<-"results_final_radius_FOM_decomposition/"
params<-read.table(paste0(path_results,"parameters_simulation.txt"),sep=",",header=TRUE)
data_sensitivity_radius_FOM_decomposition<-Rload_data(path_results,1)
palette_sensitivity_radius_FOM_decomposition<-Rmake_palettes(data_sensitivity_radius_FOM_decomposition$species_traits,
                                                             data_sensitivity_radius_FOM_decomposition$n_tot_species,
                                                             data_sensitivity_radius_FOM_decomposition$n_faeces)

# decomposition by each decomposer trophic group
data<-NULL
biomass<-NULL
for (i in 1:nrow(params)){
  # consumption
  consumption<-read.table(paste(path_results,"consumption_",params$simu_ID[i],".txt",sep=""),sep=",",header=T)
  names(consumption)<-data_sensitivity_radius_FOM_decomposition$names_tot
  consumption<-consumption[,data_sensitivity_radius_FOM_decomposition$names_detritus] # selects abiotic resources only
  consumption$consumer<-data_sensitivity_radius_FOM_decomposition$names_species
  consumption<-Rmake_table_for_plot(consumption,"consumer","resource","biomass") # makes diet a long table
  # resource type
  temp<-data_sensitivity_radius_FOM_decomposition$detritus_traits[,c("type","names")]
  names(temp)<-c("type","resource")
  consumption<-merge(consumption,temp,by="resource")
  # consumer trophic group
  temp<-data_sensitivity_radius_FOM_decomposition$species_traits[,c("trophic_group","names")]
  names(temp)<-c("trophic_group","consumer")
  consumption<-merge(consumption,temp,by="consumer")
  rm(temp)
  # only consider decomposers
  consumption<-consumption[consumption$trophic_group%in%c("microbes","macro-food web detritivores","trophic whales"),]
  consumption<-aggregate(biomass~type+trophic_group, data=consumption, sum)
  consumption$trophic_group<-factor(consumption$trophic_group,levels=c("microbes","macro-food web detritivores","trophic whales"))
  levels(consumption$trophic_group)<-c("microbes","macro-detritivores","trophic whales")
  consumption$type<-factor(consumption$type,levels=c("faeces","FOM","SOM","DOC"))
  consumption$total<-0
  for (j in 1:length(levels(consumption$type))){
    consumption$total[consumption$type==levels(consumption$type)[j]]=sum(consumption$biomass[consumption$type==levels(consumption$type)[j]])
  }
  consumption$radius_FOM<-params$radius_FOM[i]
  consumption$a_phi_FOM<-params$a_phi_FOM[i]
  data<-rbind(data,consumption)
  # biomass
  data_sensitivity_radius_FOM_decomposition<-Rload_data(path_results,i)
  temp<-data_sensitivity_radius_FOM_decomposition$data_species
  fraction<-sum(temp$biomass[temp$trophic_group=="microbes"])/sum(temp$biomass)
  biomass<-rbind(biomass,data.frame(radius_FOM=params$radius_FOM[i],
                              a_phi_FOM=params$a_phi_FOM[i],
                              fraction=fraction))
}

p1<-ggplot(data=biomass)+
  geom_raster(aes(radius_FOM,a_phi_FOM,fill=fraction))+
  theme_raster+
  scale_fill_viridis(option="viridis",
                     name="Fraction of\ntotal biomass",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)+
  x_axis_log10+
  xlab(label_radius_FOM)+
  ylab(label_a_phi_FOM)+
  ggtitle("Microbes")

p2<-ggplot(data=data[data$type%in%c("faeces","FOM")&data$trophic_group=="microbes",])+
  geom_raster(aes(radius_FOM,a_phi_FOM,fill=biomass/total))+
  facet_wrap(~type)+
  theme_raster+
  scale_fill_viridis(option="viridis",
                     name="Fraction of detritus\nconsumed by microbes",
                     limits=c(0,1),
                     breaks=seq(0,1,0.25),
                     labels = scales::percent)+
  x_axis_log10+
  xlab(label_radius_FOM)+
  ylab(label_a_phi_FOM)

graph<-ggdraw(xlim = c(0, 2.8), ylim = c(0, 1)) +
  draw_plot(p1, 0, 0, 1, 1)+
  draw_plot(p2, 1, 0, 1.8, 1)+
  draw_plot_label(c("A","B"), c(0,1), c(1,1), size = 30)
ggsave(paste(path_figures,"supp_sensitivity_FOM_radius_decomposition.pdf",sep=""), graph, width = 20, height = 6, device=cairo_pdf)

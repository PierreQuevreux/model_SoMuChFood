library(rstudioapi) # to set the working directory
library(reshape2) # to format data frames

setwd(dirname(getActiveDocumentContext()$path)) # Set working directory to source file location
path_input_data<-"input_data/"

# varied parameters
input_FOM = 300 # input of fresh organic matter
input_DOC = 150 # input of dissolved organic carbon
input_N = 500 # input of mineral nitrogen
CN_ini_FOM = 1 # C:N ratio of FOM
n_trophic_group_file = 1 # number of communities tested
trophic_group_file = 1
radius_FOM = 10^c(3,4,5,6)

var_params<-expand.grid(input_FOM=input_FOM,
                        input_DOC=input_DOC,
                        input_N=input_N,
                        CN_ini_FOM=CN_ini_FOM,
                        radius_FOM=radius_FOM,
                        trophic_group_file=trophic_group_file)
var_params$TS_record="yes" # record or not time series

sim_param<-read.table(paste(path_input_data,"parameters_simulation.csv",sep=""),sep=",",header=F, colClasses = "character", row.names=1)
sim_param<-as.data.frame(t(sim_param))
sim_param<-sim_param[,-which(names(sim_param)%in%names(var_params))]
sim_param<-merge(sim_param,var_params)
sim_param$simu_ID<-c(1:nrow(sim_param))
# write the simulation file
write.table(sim_param,paste(path_input_data,"parameters_simulation.txt",sep=""),sep=",",row.names=F,col.names=T)

# convert the tables of parameters of each trophic group into the proper format
for (i in 1:n_trophic_group_file){
  trophic_group_parameters<-read.table(paste(path_input_data,"parameters_trophic_groups_",i,".csv",sep=""),sep=",",header=F,colClasses = "character", row.names=1)
  trophic_group_parameters<-as.data.frame(t(trophic_group_parameters))
  trophic_group_parameters<-type.convert(trophic_group_parameters, as.is = TRUE) 
  write.table(trophic_group_parameters,paste(path_input_data,"parameters_trophic_groups_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
}

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
a_K_DOC=10^seq(-3,-1) # VARIATIONS OF TROPHIC_GROUP PARAMETERS
a_K_FOM=10^seq(-4,-2) # VARIATIONS OF TROPHIC_GROUP PARAMETERS
a_K_SOM=10^seq(-4,-2) # VARIATIONS OF TROPHIC_GROUP PARAMETERS

var_params<-expand.grid(input_FOM=input_FOM,
                        input_DOC=input_DOC,
                        input_N=input_N,
                        CN_ini_FOM=CN_ini_FOM,
                        a_K_DOC=a_K_DOC,
                        a_K_FOM=a_K_FOM,
                        a_K_SOM=a_K_SOM,
                        trophic_group_file=trophic_group_file)
var_params$TS_record="yes" # record or not time series
n_trophic_group_file = nrow(var_params) # number of communities tested # VARIATIONS OF TROPHIC_GROUP PARAMETERS
var_params$trophic_group_file = c(1:n_trophic_group_file) # VARIATIONS OF TROPHIC_GROUP PARAMETERS

sim_param<-read.table(paste(path_input_data,"parameters_simulation.csv",sep=""),sep=",",header=F, colClasses = "character", row.names=1)
sim_param<-as.data.frame(t(sim_param))
sim_param<-sim_param[,-which(names(sim_param)%in%names(var_params))]
sim_param<-merge(sim_param,var_params)
sim_param$simu_ID<-c(1:nrow(sim_param))
# write the simulation file
write.table(sim_param,paste(path_input_data,"parameters_simulation.txt",sep=""),sep=",",row.names=F,col.names=T)

# generate trophic_group_parameters files
trophic_group_parameters<-read.table(paste(path_input_data,"parameters_trophic_groups_1.csv",sep=""),sep=",",header=F,colClasses = "character", row.names=1)
trophic_group_parameters<-as.data.frame(t(trophic_group_parameters))
trophic_group_parameters<-type.convert(trophic_group_parameters, as.is = TRUE)
for (i in 1:n_trophic_group_file){
  trophic_group_parameters$a_K_DOC[trophic_group_parameters$trophic_group=="microbes"]=var_params$a_K_DOC[i]
  trophic_group_parameters$a_K_FOM[trophic_group_parameters$trophic_group=="microbes"]=var_params$a_K_FOM[i]
  trophic_group_parameters$a_K_SOM[trophic_group_parameters$trophic_group=="microbes"]=var_params$a_K_SOM[i]
  write.table(trophic_group_parameters,paste(path_input_data,"parameters_trophic_groups_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
}

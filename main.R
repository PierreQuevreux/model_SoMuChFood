##################
### LIBRARIRES ### ----
##################

library(doParallel) # for parallel computation
library(reshape2) # to format data frames

source("internal_functions_R.R")

path_input_data<-"input_data/"
path_results<-"results/"

###################
### SIMULATIONS ### ----
###################

sim_param<-read.table(paste(path_input_data,"parameters_simulation.txt",sep=""),sep=",",header=T)

# simulations
start_time<-Sys.time()
print(start_time)
n_simu=dim(sim_param)[1]
no_cores <- detectCores()# - 1
cl <- makeCluster(no_cores)  
registerDoParallel(cl)
results<-foreach(i=1:n_simu) %dopar% Rfoodweb_simulation(as.list(sim_param[i,]))
stopCluster(cl)
end_time<-Sys.time()
print(end_time)
end_time-start_time

Rsave_data_tables(results,sim_param,n_simu,path_results) # write the results in .txt files
library(stats) # for the integrate function

#####################
### MAIN FUNCTION ### ----
#####################

Rfoodweb_simulation<-function(foodweb_inputs){
  #attach(foodweb_inputs) # FOR TEST
  #attach(trophic_group_parameters) # FOR TEST
  with(foodweb_inputs,{
    #################################
    ### IMPORTING food web_inputs ### ----
    #################################
    # food web_inputs contains several lists with all the parameters # ----
    ### sim_param: parameters specific of the simulation
    ### trophic_group_parameters: parameters of species
    # simulation parameters (sim_param)
    ### dimension: number of dimensions describing the food web (trophic group, body mass, C:N)
    ### n_detritus_pool: number of detritus pool (faeces excluded), which are FOM, SOM and DOC
    ### n_nutrients: number of different nutrients
    ### growth_model: model of primary production
    ### predation_model: model of predation
    ### deltaT: time step expressed in days
    ### n_step_dynamics: total number of time steps
    ### n_step_output: number of time steps (-1) between the recorded time steps
    n_parameters=length(foodweb_inputs)
    trophic_group_parameters<-read.table(paste(path_input_data,"parameters_trophic_groups_",trophic_group_file,".txt",sep=""),sep=",",header=T,stringsAsFactors=T)
    n_trophic_group=dim(trophic_group_parameters)[1] # number of trophic groups (organisms from a trophic group have the same physiological parameters)
    foodweb_inputs$n_trophic_group=n_trophic_group
    #print("OK")
    # input files checking # ----
    Rcheck_species_presence(trophic_group_parameters,n_trophic_group) # check the values of n_mass and n_CN
    Rcheck_age_structured(trophic_group_parameters,n_trophic_group) # check is age structured species have more than 1 mass class
    Rcheck_time_steps(T_simu,T_trans,T_step,n_step_output) # check the validity of the simulated and recorded time steps
    # body mass and C:N distributions # ----
    # trophic species are defined by three dimensions: trophic group, body mas and C:N (and a trophic type is associated to each trophic group)
    # body mass (fresh body mass which is converted into carbon for allometric parameter calculation)
    n_mass_list=array(0,n_trophic_group) # array containing the number of mass intervals for each trophic group
    species_logbodymass_list=vector("list",n_trophic_group) # list containing the array of minimal values of each body mass interval
    species_bodymass_average_list=vector("list",n_trophic_group) # list containing the array of average body masses of each interval
    # stoichiometry
    n_CN_list=array(0,n_trophic_group) # array containing the number of C:N for each trophic group
    species_CN_list=vector("list",n_trophic_group) # list containing the array of minimal values of each CN interval
    species_CN_average_list=vector("list",n_trophic_group) # list containing the array of C:N for each CN interval
    # initialisation of the arrays contained by lists
    for (i in 1:n_trophic_group){
      # body mass (fresh) # ----
      n_mass_list[i]=trophic_group_parameters$n_mass[i] # number of mass intervals for each trophic group 
      species_logbodymass_list[[i]]=array(0,(n_mass_list[i]+1)) # array of minimal values of each body mass interval (+1 because the width of the interval is stored in the last element of the vector)
      if (n_mass_list[i]>1){
        species_logbodymass_list[[i]][(n_mass_list[i]+1)]=(trophic_group_parameters$logmass_max[i]-trophic_group_parameters$logmass_min[i])/(n_mass_list[i]-1.0) # body mass interval
      } else{
        species_logbodymass_list[[i]][(n_mass_list[i]+1)]=1.0 # if no subdivision into intervals
      }
      species_bodymass_average_list[[i]]=array(0,n_mass_list[i]) # array of average body masses of each interval
      for (j in 1:n_mass_list[i]){
        species_logbodymass_list[[i]][j]=trophic_group_parameters$logmass_min[i]+(j-1)*species_logbodymass_list[[i]][(n_mass_list[i]+1)] # lower value of each body mass interval
      }
      # C:N # ----
      n_CN_list[i]=trophic_group_parameters$n_CN[i] # number of C:N intervals for each trophic group
      species_CN_list[[i]]=array(0,(n_CN_list[i]+1)) # array of minimal values of each C:N interval (+1 because the width of the interval is stored in the last element of the vector)
      if (n_CN_list[i]>1){
        species_CN_list[[i]][(n_CN_list[i]+1)]=(trophic_group_parameters$CN_max[i]-trophic_group_parameters$CN_min[i])/(n_CN_list[i]-1.0)
      } else{
        species_CN_list[[i]][(n_CN_list[i]+1)]=1.0 # if no subdivision into intervals
      }
      species_CN_average_list[[i]]=array(0,n_CN_list[i]) # array of average C:N
      for (j in 1:n_CN_list[i]){
        species_CN_list[[i]][j]=trophic_group_parameters$CN_min[i]+(j-1)*species_CN_list[[i]][(n_CN_list[i]+1)] # lower value of each C:N interval
      }
    }
    # computation of average (fresh) body mass and C:N
    species_bodymass_average_list=Rmake_average_bodymass_list(species_logbodymass_list,n_mass_list,n_trophic_group) # computes the average body mass if n_mass_list[i]>1
    species_CN_average_list=Rmake_average_CN_list(species_CN_list,n_CN_list,trophic_group_parameters$CN_min,trophic_group_parameters$CN_max,trophic_group_parameters$CN_sigma,n_trophic_group)
    n_tot_species=sum(n_mass_list*n_CN_list) # total number of trophic species
    foodweb_inputs$n_tot_species=n_tot_species
    
    # defines detritus compartments # ----
    n_faeces=sum(n_mass_list) # temporary number of faeces pools
    faeces_pos=as.data.frame(matrix(0,n_faeces,4)) #
    pos=0 # index
    for(i in 1:n_trophic_group){
      for(j in 1:n_mass_list[i]){
        if(n_mass_list[i]>0){
        pos=pos+1
        faeces_pos[pos,1]=species_bodymass_average_list[[i]][j] # size of organisms producing the faeces
        faeces_pos[pos,2]=i # trophic group index
        faeces_pos[pos,3]=j # body mass index
        }
      }
    }
    faeces_pos[,4]=faeces_pos[,1] # value of faeces size
    faeces_pos[,1]=as.factor(faeces_pos[,1]) # first column serving as ID to identify faeces of the same size class
    n_faeces=length(levels(faeces_pos[,1])) # updates the number of faeces pools (faeces produced by different species but with the same size belong to the same pool)
    levels(faeces_pos[,1])=seq(1:n_faeces) # attribute an index to the detritus with the same size
    n_tot_detritus=n_faeces+n_detritus_pool # total number of detritus (faeces+FOM+SOM+DOC)
    # index of the various detritus and nutrient compartments
    n_tot=n_tot_species+n_tot_detritus+n_nutrients # total number of compartments in the ecosystem
    if (n_faeces>0){
      faeces_index_detritus=1:n_faeces
    }else(faeces_index_detritus=0)
    FOM_index_detritus=n_faeces+1 # position of FOM in the detritus vector
    SOM_index_detritus=n_faeces+2 # position of SOM in the detritus vector
    DOC_index_detritus=n_faeces+3 # position of DOC in the detritus vector
    faeces_index=n_tot_species+faeces_index_detritus # position of faeces
    FOM_index=n_tot_species+FOM_index_detritus # position of FOM
    SOM_index=n_tot_species+SOM_index_detritus # position of SOM
    DOC_index=n_tot_species+DOC_index_detritus # position of DOC
    detritus_index=c(n_tot_species+(1:n_tot_detritus)) # position of detritus
    nutrients_index=c(n_tot_species+n_tot_detritus+(1:n_nutrients)) # position of nutrients
    N_index=n_tot_species+n_tot_detritus+1 # position of nitrogen
    # update the parameter list
    foodweb_inputs$n_tot=n_tot
    foodweb_inputs$n_faeces=n_faeces
    foodweb_inputs$n_tot_detritus=n_tot_detritus
    foodweb_inputs$faeces_index_detritus=faeces_index_detritus
    foodweb_inputs$FOM_index_detritus=FOM_index_detritus
    foodweb_inputs$SOM_index_detritus=SOM_index_detritus
    foodweb_inputs$DOC_index_detritus=DOC_index_detritus
    foodweb_inputs$faeces_index=faeces_index # position of faeces
    foodweb_inputs$FOM_index=FOM_index
    foodweb_inputs$SOM_index=SOM_index
    foodweb_inputs$DOC_index=DOC_index
    foodweb_inputs$detritus_index=detritus_index # position of detritus
    foodweb_inputs$nutrients_index=nutrients_index # position of nutrients
    foodweb_inputs$N_index=N_index
    ## detritus vectors (+ mineral nutrients)
    detritus_traits=as.data.frame(matrix(0,n_tot_detritus,3))
    names(detritus_traits)=c("bodymass","radius","type")
    detritus_traits$bodymass=c(sort(unique(faeces_pos[,4])),size_FOM,size_SOM,0) # faeces ordered by increasing size (ordering also done by as.factor)
    detritus_traits$radius=c(Rmake_faeces_radius(detritus_traits$bodymass[faeces_index_detritus],r_int,r_slope),radius_FOM,radius_SOM,0) # compute the radius of detritus pellets according to their "body mass"
    detritus_traits$bodymass[1:n_faeces]=Rmake_faeces_mass(detritus_traits$bodymass[1:n_faeces],faeces_bodymass_ratio)
    detritus_traits$type=c(rep("faeces",n_faeces),"FOM","SOM","DOC")
    ## fraction available
    detritus_fraction_available<-vector("list",n_trophic_group)
    for (i in 1:n_trophic_group){
      if(trophic_group_parameters$trophic_type[i]=="microbes"){
        detritus_fraction_available[[i]]=Rmake_fraction_available_detritus(detritus_traits$radius,trophic_group_parameters$radius_available[i]) # fraction of detritus available for decomposition
        detritus_fraction_available[[i]][DOC_index_detritus]=0 # DOC is dissolved (no surface)
      }
    }
    
    # defines species # ----
    species_traits<-as.data.frame(matrix(nrow=n_tot_species,ncol=9))
    names(species_traits)=c("trophic_group","trophic_group_index","trophic_type","recycling","population_structure","bodymass","CN","faeces_pos","fresh_to_C")
    species_traits$trophic_group=factor(species_traits$trophic_group,levels=levels(trophic_group_parameters$trophic_group))
    species_traits$trophic_type=factor(species_traits$trophic_type,levels=levels(trophic_group_parameters$trophic_type))
    species_traits$population_structure=factor(species_traits$population_structure,levels=levels(trophic_group_parameters$population_structure))
    #recycling=vector(0,n_tot_species)
    pos=0 # index
    pos_f=0 # index in the faeces position matrix
    for (i in 1:n_trophic_group){
      for (j in 1:n_mass_list[i]){
        if(n_mass_list[i]>0){
          pos_f=pos_f+1
          for (k in 1:n_CN_list[i]){
            pos=pos+1
            species_traits$trophic_group[pos]=trophic_group_parameters$trophic_group[i] # microbes, micro-food web, macro-food web...
            species_traits$trophic_group_index[pos]=i # index of the trophic group to navigate in trophic_group_parameters
            species_traits$trophic_type[pos]=trophic_group_parameters$trophic_type[i] # microbe, detritivore, carnivore
            species_traits$recycling[pos]=trophic_group_parameters$recycling[i] # fraction of direct recycling
            species_traits$population_structure[pos]=trophic_group_parameters$population_structure[i] # age structured or unstructured
            species_traits$bodymass[pos]=species_bodymass_average_list[[i]][j] # average body mass
            species_traits$CN[pos]=species_CN_average_list[[i]][k] # average C:N
            species_traits$faeces_pos[pos]=as.numeric(faeces_pos[pos_f,1]) # associated faeces pool index
            species_traits$fresh_to_C[pos]=trophic_group_parameters$fresh_to_dry[i]*trophic_group_parameters$dry_to_C[i] # conversion from fresh body mass to carbon
          }
        }
      }
    }
    # index of the various species classes
    if (n_tot_species>0){
      species_index=c(1:n_tot_species) # position of species
    }else{species_index=integer(0)}
    microbes_index=which(species_traits$trophic_type=="microbes") # position of microbes
    microbivores_index=which(species_traits$trophic_type=="microbivores") # position of microbivores
    detritivores_index=which(species_traits$trophic_type=="detritivores") # position of detritivores
    carnivores_index=which(species_traits$trophic_type=="carnivores") # position of carnivores
    animals_index=species_index[which((species_index%in%microbes_index)==FALSE)] # position of non-microbe species
    age_structured_index=which(trophic_group_parameters$population_structure=="age_structured") # position of age structured trophic groups
    if (length(age_structured_index)>0){
      for (i in 1:length(age_structured_index)){
        if (n_mass_list[age_structured_index[i]]<2){
          age_structured_index<-age_structured_index[-i] # does not consider the species as age structured if there is only one size class
        }
      }
    }
    trophic_group_index_list=vector("list",n_trophic_group)
    for (i in 1:n_trophic_group){
      trophic_group_index_list[[i]]=which(species_traits$trophic_group_index==i) # vector with the ID of each species of trophic group i
    }
    # update the parameter list
    foodweb_inputs$species_index=species_index # position of species
    foodweb_inputs$microbes_index=microbes_index # position of microbes
    foodweb_inputs$microbivores_index=microbivores_index # position of microbes
    foodweb_inputs$detritivores_index=detritivores_index # position of detritivores
    foodweb_inputs$carnivores_index=carnivores_index # position of carnivores
    foodweb_inputs$animals_index=animals_index # position of non-microbe species
    foodweb_inputs$age_structured_index=age_structured_index # position of age structured trophic groups
    foodweb_inputs$trophic_group_index_list=trophic_group_index_list
    
    # vectors of variables # ----
    # density of biomass, detritus and nutrients
    density_vector=array(0,n_tot)
    pos=0 # index
    for (i in 1:n_trophic_group){
      for (j in 1:n_mass_list[i]){
        if(n_mass_list[i]>0){
          for (k in 1:n_CN_list[i]){
            pos=pos+1
            density_vector[pos]=trophic_group_parameters$initial_density[i]/(n_mass_list[i]*n_CN_list[i]) # initial biomass density split into the different species
          }
        }
      }
    }
    density_vector[faeces_index]=density_ini_faeces/n_faeces # initial density of faeces
    density_vector[FOM_index]=density_ini_FOM # initial density of FOM
    density_vector[SOM_index]=density_ini_SOM # initial density of SOM
    density_vector[DOC_index]=density_ini_DOC # initial density of DOC
    density_vector[N_index]=density_ini_N # initial density of nitrogen
    
    # C:N vector
    CN_vector=array(0,n_tot)
    CN_vector[species_index]=species_traits$CN # C:N of species
    CN_vector[faeces_index]=CN_ini_faeces # initial C:N of faeces
    CN_vector[FOM_index]=CN_ini_FOM # initial C:N of FOM
    CN_vector[SOM_index]=CN_ini_SOM # initial C:N of SOM
    CN_vector[DOC_index]=Inf # DOC is pure carbon
    
    # fraction of the population of each species that is active (feeding and respiration)
    active_vector=array(1,n_tot_species) # the entire population is set active
    
    # matrix of species with age structure
    age_structure_list=Rmake_age_structure_list(n_mass_list,n_CN_list,trophic_group_parameters,foodweb_inputs) # matrix with age structure of populations (each column for on species and each row for a life stage)
    age_factors_vector=Rmake_age_factors_vector(n_mass_list,species_logbodymass_list,trophic_group_parameters,foodweb_inputs)
    
    # vectors of the mean values of the recorded variables
    TL_vector=array(0,n_tot) # vector with the trophic level of each species (even nonliving compartments)
    TL_vector[species_index]=1 # initialisation of species (0 for nonliving compartments)
    TL_vector[microbivores_index]=2
    density_mean=array(0,n_tot) # mean density of each compartment
    active_mean=array(0,n_tot_species) # mean activity coefficient
    CN_mean=array(0,n_tot) # mean C:N ratio of each compartment
    limitation_index_mean=array(0,n_tot_species) # mean quantification of limitation by C or N
    consumption_mean=matrix(0,n_tot_species,n_tot) # mean biomass of each resource that is consumed
    diet_mean=matrix(0,n_tot_species,n_tot) # mean biomass of each resource that is assimilated
    decomposition_mean=array(0,n_tot_detritus) # mean decomposition of each detritus compartment
    respiration_mean=array(0,n_tot_species) # mean respiration of each species
    N_mineralisation_mean=array(0,n_tot_species) # mean quantity of nutrients mineralised by each species
    C_burrowing_mean=0 # mean quantity of carbon definitively stored in deep soil
    detritus_production_mean=array(0,4) # mean production rate of faeces, FOM, SOM and DOC by organisms
    
    ################################################
    ### COMPUTING INSTANTANEOUS BIOLOGICAL RATES ### ----
    ## that are independent of population densities
    ## and thus only needs to be computed once
    ################################################
    
    metabolic_rate_vector=Rmake_metabolic_rate_vector(species_traits,trophic_group_parameters,foodweb_inputs)
    mortality_rate_vector=Rmake_mortality_rate_vector(species_traits,trophic_group_parameters,foodweb_inputs)
    self_regulation_rate_vector=Rmake_self_regulation_rate_vector(species_traits,trophic_group_parameters,foodweb_inputs)
    shared_self_regulation_rate_vector=Rmake_shared_self_regulation_rate_vector(species_traits,trophic_group_parameters,foodweb_inputs)
    interference_rate_vector=Rmake_interference_rate_vector(species_traits,trophic_group_parameters,foodweb_inputs)
    attack_rate=Rmake_attack_rate(species_traits,detritus_traits,trophic_group_parameters,foodweb_inputs)
    attack_rate_base_matrix=attack_rate$attack_rate_base_matrix
    preference_bodymass_matrix=attack_rate$preference_bodymass_matrix
    preference_CN_matrix=attack_rate$preference_CN_matrix
    attack_rate_matrix=attack_rate$attack_rate_matrix
    handling_time_matrix=Rmake_handling_time_matrix(species_traits,trophic_group_parameters,foodweb_inputs)
    max_uptake_matrix=Rmake_max_uptake_matrix(trophic_group_parameters,species_traits,foodweb_inputs)
    half_saturation_matrix=Rmake_half_saturation_matrix(trophic_group_parameters,species_traits,foodweb_inputs)
    assimilation_efficicency_matrix=Rmake_assimilation_efficicency_matrix(species_traits,trophic_group_parameters,foodweb_inputs)
    dormancy_rate_vector=Rmake_dormancy_rate(species_traits,trophic_group_parameters,foodweb_inputs)
    dormancy_asymmetry_vector=Rmake_dormancy_asymmetry(species_traits,trophic_group_parameters,foodweb_inputs)
    light_interaction_list=Rmake_light_interaction_list(attack_rate_matrix,foodweb_inputs)
    reverse_interaction_list=Rmake_reverse_interaction_list(attack_rate_matrix,foodweb_inputs)
    
    # adds the biological rates and parameter structures to foodweb_inputs
    biological_rates=list(metabolic_rate_vector=metabolic_rate_vector,
                          mortality_rate_vector=mortality_rate_vector,
                          self_regulation_rate_vector=self_regulation_rate_vector,
                          shared_self_regulation_rate_vector=shared_self_regulation_rate_vector,
                          interference_rate_vector=interference_rate_vector,
                          attack_rate_base_matrix=attack_rate_base_matrix,
                          preference_bodymass_matrix=preference_bodymass_matrix,
                          preference_CN_matrix=preference_CN_matrix,
                          attack_rate_matrix=attack_rate_matrix,
                          handling_time_matrix=handling_time_matrix,
                          max_uptake_matrix=max_uptake_matrix,
                          half_saturation_matrix=half_saturation_matrix,
                          assimilation_efficicency_matrix=assimilation_efficicency_matrix,
                          dormancy_rate_vector=dormancy_rate_vector,
                          dormancy_asymmetry_vector=dormancy_asymmetry_vector,
                          light_interaction_list=light_interaction_list,
                          reverse_interaction_list=reverse_interaction_list,
                          age_structure_list=age_structure_list,
                          age_factors_vector=age_factors_vector,
                          species_logbodymass_list=species_logbodymass_list,
                          detritus_fraction_available=detritus_fraction_available,
                          n_mass_list=n_mass_list,
                          n_CN_list=n_CN_list)
    
    #######################
    ## FOOD WEB DYNAMICS ## ----
    #######################
    # creates output lists # ----
    output_step=list(density_vector=density_vector,
                     CN_vector=CN_vector,
                     active_vector=active_vector)
    if (T_ini>0){
      simu_step_ini=seq(1,floor(T_ini/T_step_ini)) # number of time steps in the very initial dynamics (to avoid violent variations)
      simu_t_steps=c(seq(T_step_ini,T_ini,T_step_ini),seq(T_ini+T_step,T_simu,T_step)) # time steps
    }else{
      simu_step_ini=integer(0)
      simu_t_steps=seq(T_step,T_simu,T_step) # time steps
    }
    simu_step=c(simu_step_ini,seq(1,floor((T_simu-T_ini)/T_step))+length(simu_step_ini)) # total number of time steps after the initial dynamics
    deltaT=rep(T_step,length(simu_step)) # regular deltaT
    deltaT[simu_step_ini]=T_step_ini # deltaT during the initial phase
    simu_step_record=length(simu_t_steps[simu_t_steps>T_trans]) # number of recorded time steps
    if (simu_step_record<n_step_output){
      n_step_output=simu_step_record # records all the times steps (number bellow the value specified in the inputs)
    }
    if (T_trans>T_ini){ # no initial dynamics are recorded
      recorded_time_step=seq(T_trans+T_step,T_simu,(T_simu-T_trans)/(n_step_output)) # recorded time steps
    }else{ # recording of time steps of initial dynamics
      if (n_step_output>(T_simu-T_trans)/T_step){
        n_step_output=(T_simu-T_trans)/T_step # records less time steps if not enough time steps are calculated
      }
      recorded_time_step=seq(T_trans+T_step_ini,T_simu,(T_simu-T_trans)/(n_step_output)) # recorded time steps
    }
    recorded_time_step=c(recorded_time_step,T_simu+1) # this last time step is not meant to be recorded
    times_series=vector('list',n_step_output) # final output for analysis
    time=0 # time
    rec=1 # recording counter
    # simulation
    #deltaT=deltaT[1] # FOR TEST
    #for (i in 1:180){ # FOR TEST
    # temporal dynamics # ----
    for (i in simu_step){
      time=simu_t_steps[i]
      output_step=Rcomputation_foodweb_step_dynamics(output_step$density_vector,
                                                     output_step$CN_vector,
                                                     output_step$active_vector,
                                                     output_step$TL_vector,
                                                     species_traits,
                                                     detritus_traits,
                                                     biological_rates,
                                                     deltaT[i],
                                                     trophic_group_parameters,
                                                     foodweb_inputs) # compute the dynamics at one time step
      #print(i)
      if (time>=recorded_time_step[rec]-deltaT[i]/2){ # R is not very precise and equality is not systematically met
        #print(time)
        output_step$time=time
        TL_vector=Rcomputation_TL(TL_vector,output_step$diet,foodweb_inputs) # computes trophic levels recursively
        density_mean=density_mean+output_step$density_vector
        CN_mean=CN_mean+output_step$CN_vector
        active_mean=active_mean+output_step$active_vector
        limitation_index_mean=limitation_index_mean+output_step$limitation_index
        consumption_mean=consumption_mean+output_step$consumption_matrix
        diet_mean=diet_mean+output_step$diet_matrix
        decomposition_mean=decomposition_mean+output_step$decomposition
        respiration_mean=respiration_mean+output_step$respiration
        N_mineralisation_mean=N_mineralisation_mean+output_step$N_mineralisation
        C_burrowing_mean=C_burrowing_mean+output_step$C_burrowing
        detritus_production_mean=detritus_production_mean+output_step$detritus_production
        if (TS_record=="yes"){
          times_series[[rec]]=output_step
        }
        rec=rec+1
      }
    }
    # FOR TEST #
    #density_vector=output_step$density_vector # FOR TEST
    #CN_vector=output_step$CN_vector # FOR TEST
    #active_vector=output_step$active_vector # FOR TEST
    #TL_vector=output_step$TL_vector # FOR TEST
    
    # formatting of the outputs # ----
    if (n_tot_species>0){
      species_traits$ID=1
    }else{species_traits$ID=integer(0)}
    if (n_tot_species>=2){
      for (i in 2:n_tot_species){
        if (species_traits$trophic_group[i]==species_traits$trophic_group[i-1]){
          species_traits$ID[i]=species_traits$ID[i-1]+1
        }
      }
    }
    species_traits$names<-paste(species_traits$trophic_group,species_traits$ID,sep=" ")
    if (n_faeces>0){
      names_faeces=paste("faeces_",seq(1:n_faeces),sep="")
    }else{names_faeces=character(0)}
    names_tot=c(species_traits$names,names_faeces,"FOM","SOM","DOC","N")
    names_detritus<-names_tot[detritus_index]
    names_nutrients<-names_tot[nutrients_index]
    detritus_traits$names<-names_detritus
    
    density_mean=density_mean/n_step_output
    CN_mean=CN_mean/n_step_output
    active_mean=active_mean/n_step_output
    limitation_index_mean=limitation_index_mean/n_step_output
    consumption_mean=consumption_mean/n_step_output
    diet_mean=diet_mean/n_step_output
    decomposition_mean=decomposition_mean/n_step_output
    respiration_mean=respiration_mean/n_step_output
    N_mineralisation_mean=N_mineralisation_mean/n_step_output
    C_burrowing_mean=C_burrowing_mean/n_step_output
    detritus_production_mean=detritus_production_mean/n_step_output
    
    TL_vector<-TL_vector[species_index]
    
    if (TS_record=="yes"){
      TS_density=Rmake_time_series(times_series,n_step_output,"density_vector",names_tot)
      TS_CN=Rmake_time_series(times_series,n_step_output,"CN_vector",names_tot)
      TS_active=Rmake_time_series(times_series,n_step_output,"active_vector",species_traits$names)
      TS_limitation_index=Rmake_time_series(times_series,n_step_output,"limitation_index",species_traits$names)
      TS_decomposition=Rmake_time_series(times_series,n_step_output,"decomposition",c(names_faeces,"FOM","SOM","DOC"))
      TS_respiration=Rmake_time_series(times_series,n_step_output,"respiration",species_traits$names)
      TS_N_mineralisation=Rmake_time_series(times_series,n_step_output,"N_mineralisation",species_traits$names)
    }else{
      TS_density=NULL
      TS_CN=NULL
      TS_active=NULL
      TS_limitation_index=NULL
      TS_decomposition=NULL
      TS_respiration=NULL
      TS_N_mineralisation=NULL
    }
    
    info=list(n_tot=n_tot,
              n_tot_species=n_tot_species,
              n_faeces=n_faeces,
              n_detritus_pool=n_detritus_pool,
              n_tot_detritus=n_tot_detritus,
              n_nutrients=n_nutrients,
              n_parameters=n_parameters,
              names_species=species_traits$names,
              names_faeces=names_faeces,
              names_detritus=names_detritus,
              names_nutrients=names_nutrients,
              names_tot=names_tot)
    
    output_final=list(TS_density=TS_density,
                      TS_CN=TS_CN,
                      TS_active=TS_active,
                      TS_limitation_index=TS_limitation_index,
                      TS_decomposition=TS_decomposition,
                      TS_respiration=TS_respiration,
                      TS_N_mineralisation=TS_N_mineralisation,
                      density_mean=density_mean,
                      CN_mean=CN_mean,
                      TL_vector=TL_vector,
                      active_mean=active_mean,
                      limitation_index_mean=limitation_index_mean,
                      consumption_mean=consumption_mean,
                      diet_mean=diet_mean,
                      decomposition_mean=decomposition_mean,
                      respiration_mean=respiration_mean,
                      N_mineralisation_mean=N_mineralisation_mean,
                      C_burrowing_mean=C_burrowing_mean,
                      detritus_production_mean=detritus_production_mean,
                      species_traits=species_traits,
                      detritus_traits=detritus_traits,
                      info=info)
    return(output_final)
  })
}

###############################################################################
### FUNCTIONS OF BIOLOGICAL RATES (computed before time series simulations) ### ----
###############################################################################
# BODY MASS # ----
# average body mass of each size class of a body mass distribution
Rmake_average_bodymass<-function(species_logbodymass,n_mass){
  if (n_mass>1){
    if (species_logbodymass[(n_mass+1)]==0){
      return(rep(10^species_logbodymass[1],n_mass)) # if all species have the same size
    } else{
    return(10^species_logbodymass[1:n_mass]*(10^species_logbodymass[(n_mass+1)]-1)/(species_logbodymass[(n_mass+1)]*log(10))) # array with the average body mass of each interval
    }
  } else{
    return(10^species_logbodymass[1])
  }
}
# average body mass of each class of an entire interval for each trophic group
Rmake_average_bodymass_list<-function(species_logbodymass_list,n_mass_list,n_trophic_group){
  average_bodymass_list=vector("list", n_trophic_group)
  for (i in 1:n_trophic_group){
    average_bodymass_list[[i]]=Rmake_average_bodymass(species_logbodymass_list[[i]],n_mass_list[i])
  }
  return(average_bodymass_list)
}

# C:N # ----
# distribution of C:N
Rdistribution_CN<-function(x,mean,sd){
  return(x*dnorm(x,mean,sd)) # normal distribution of C:N (mean determined by C:N_min and C:N_max)
}
# average C:N of a stoichiometric class
Rmake_average_CN<-function(species_CN,n_CN,mean_CN,CN_sigma){
  if (n_CN>1){
    average_CN=array(0,n_CN) # array containing the average C:N of each interval
    for(j in 1:n_CN){
      average_CN[j]=integrate(Rdistribution_CN,
                       lower=species_CN[j],
                       upper=(species_CN[j]+species_CN[n_CN+1]),
                       mean=mean_CN,sd=CN_sigma)$value/
        integrate(dnorm, # rescaling by the probabilities of the considered interval
                  lower=species_CN[j],
                  upper=(species_CN[j]+species_CN[n_CN+1]),
                  mean=mean_CN,sd=CN_sigma)$value
    }
    return(average_CN)
  } else{
    return(species_CN[1])
  }
}
# average CN of each class of an entire interval for each trophic group
Rmake_average_CN_list<-function(species_CN_list,n_CN_list,CN_min,CN_max,CN_sigma,n_trophic_group){
  average_CN_list=vector("list", n_trophic_group)
  mean_CN=0
  for (i in 1:n_trophic_group){
    mean_CN=(CN_max[i]+CN_min[i])/2
    average_CN_list[[i]]=Rmake_average_CN(species_CN_list[[i]],n_CN_list[i],mean_CN,CN_sigma[i])
  }
  return(average_CN_list)
}

# FAECES MASS # ----
# average "body mass" of the faeces produced by each size class
Rmake_faeces_mass<-function(species_bodymass_average,faeces_bodymass_ratio){
    faeces_mass=species_bodymass_average*faeces_bodymass_ratio
    return(faeces_mass)
}
# radius of faeces available for microbes
Rmake_faeces_radius<-function(species_bodymass_average,r_int,r_slope){
    faeces_radius=r_int+r_slope*species_bodymass_average
    return(faeces_radius)
}
# fraction of detritus available for microbes
Rmake_fraction_available_detritus<-function(radius_faeces,radius_available){
  fraction_available=array(0,length(radius_faeces))
  for (i in 1:length(radius_faeces)){
    fraction_available[i]=min(1-((radius_faeces[i] - radius_available) / radius_faeces[i])^3,1)
  }
    return(fraction_available)
}

# AGE CLASS MATRIX # ----
# matrix with the index of each age class of each trophic species
Rmake_age_structure_matrix<-function(n_mass,n_CN,pos){
  age_structure_matrix=matrix(0,n_mass,n_CN) # body mass in rows and C:N in columns
  # pos : position of each "species" in the density vector
  for (i in 1:n_mass){
    if(n_mass>0){
      for (j in 1:n_CN){
        pos=pos+1 # species are ordered by trophic group, body mass class and C:N
        age_structure_matrix[i,j]=pos # each row corresponds to an age class and columns to the different species (defined by different C:N)
      }
    }
  }
  return(age_structure_matrix)
}
# list with the position matrix of each trophic group with age structured populations
Rmake_age_structure_list<-function(n_mass_list,n_CN_list,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    age_structure_list=vector("list",n_trophic_group)
    pos=0 # position in the density vector
    for (i in 1:n_trophic_group){
      if(trophic_group_parameters$population_structure[i]=="age_structured"){
        age_structure_list[[i]]=Rmake_age_structure_matrix(n_mass_list[i],n_CN_list[i],pos)
      }
      pos=pos+n_mass_list[i]*n_CN_list[i]
    }
    return(age_structure_list)
  })
}
# transition factors between age classes
Rmake_age_factors_vector<-function(n_mass_list,species_logbodymass_list,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    age_factors_vector=array(0,n_trophic_group)
    for (i in 1:n_trophic_group){
      if(trophic_group_parameters$population_structure[i]=="age_structured"){
        age_factors_vector[i]=10^(species_logbodymass_list[[i]][n_mass_list[i]+1])/(10^(species_logbodymass_list[[i]][n_mass_list[i]+1])-1)
      }
    }
  return(age_factors_vector)
  })
}

# METABOLISM # ----
# instantaneous metabolic rate of a trophic species
Rmake_metabolic_rate<-function(bodymass,kB,temperature,parameters){
  with(parameters,{
    conversion=fresh_to_dry*dry_to_C
    return(R0*exp(-E_R/(kB*temperature))*(bodymass*conversion)^s_R)
  })
}
# instantaneous metabolic rate of each class for each trophic group
Rmake_metabolic_rate_vector<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    metabolic_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      metabolic_rate_vector[i]=Rmake_metabolic_rate(species_traits$bodymass[i],kB,temperature,parameters)
    }
    return(metabolic_rate_vector)
  })
}

# MORTALITY # ----
# instantaneous density-independent mortality rate of a trophic species
Rmake_mortality_rate<-function(bodymass,kB,temperature,parameters){
  with(parameters,{
    conversion=fresh_to_dry*dry_to_C
    return(mu0*exp(-E_mu/(kB*temperature))*(bodymass*conversion)^s_mu)
  })
}
# instantaneous mortality rate of each class for each trophic group
Rmake_mortality_rate_vector<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    mortality_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      mortality_rate_vector[i]=Rmake_mortality_rate(species_traits$bodymass[i],kB,temperature,parameters)
    }
    return(mortality_rate_vector)
  })
}

# SELF-REGULATION # ----
# instantaneous density-dependent mortality (self-regulation) rate of a trophic species
Rmake_self_regulation_rate<-function(bodymass,kB,temperature,parameters){
    with(parameters,{
      conversion=fresh_to_dry*dry_to_C
      return(D0*exp(-E_D/(kB*temperature))*(bodymass*conversion)^s_D)
  })
}
# instantaneous self-regulation of each class for each trophic group
Rmake_self_regulation_rate_vector<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    self_regulation_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      self_regulation_rate_vector[i]=Rmake_self_regulation_rate(species_traits$bodymass[i],kB,temperature,parameters)
    }
    return(self_regulation_rate_vector)
  })
}
# shared self-regulation for each trophic group   ???????????? A adapter
Rmake_shared_self_regulation_rate_vector<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    shared_self_regulation_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      #shared_self_regulation_rate_vector[i]=Rmake_self_regulation_rate(species_traits$bodymass[i],parameters)
    }
    return(shared_self_regulation_rate_vector)
  })
}

# INTERFERENCE # ----
# interference between consumers
Rmake_interference_rate<-function(bodymass,kB,temperature,parameters){
  with(parameters,{
    conversion=fresh_to_dry*dry_to_C
    return(c0*exp(-E_c/(kB*temperature))*(bodymass*conversion)^s_c)
  })
}
# interference of each class for each trophic group
Rmake_interference_rate_vector<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    interference_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      interference_rate_vector[i]=Rmake_interference_rate(species_traits$bodymass[i],kB,temperature,parameters)
    }
    return(interference_rate_vector)
  })
}

# ATTACK RATE (CARNIVORES AND DETRITIVORES) # ----
# prey-independent basal attack rate of predators
Rmake_attack_rate_base<-function(bodymass,kB,temperature,parameters){
  with(parameters,{
    conversion=fresh_to_dry*dry_to_C
    return(a_a*a0*exp(-E_a/(kB*temperature))*(bodymass*conversion)^s_a) # a_a is a tuning parameter
  })
}
# prey-independent basal attack rate of predators of each class for each trophic group
Rmake_attack_rate_base_matrix<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    attack_rate_base_matrix=matrix(0,n_tot_species,n_tot)
    # predatory interactions
    for (i in c(microbivores_index,carnivores_index)){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      attack_rate_base_matrix[i,species_index]=Rmake_attack_rate_base(species_traits$bodymass[i],kB,temperature,parameters)
    }
    # detritivore interactions (faeces and FOM consumption)
    for (i in detritivores_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      attack_rate_base_matrix[i,c(faeces_index,FOM_index)]=Rmake_attack_rate_base(species_traits$bodymass[i],kB,temperature,parameters)
    }
    return(attack_rate_base_matrix)
  })
}
# prey preference depending on body mass ratio (for species and detritus)
Rmake_preference_bodymass<-function(prey_bodymass,pred_bodymass,parameters){
  with(parameters,{
    return(exp(-((log(prey_bodymass / pred_bodymass) - log(pref_theta)) / pref_sigma)^2)) # use of fresh body mass
  })
}
# matrix with all prey preference depending on body mass ratio (for species and detritus)
Rmake_preference_bodymass_matrix<-function(species_traits,detritus_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    preference_bodymass_matrix=matrix(0,n_tot_species,n_tot) # matrix with all trophic species and detritus
    # microbivory interactions
    for (i in microbivores_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      preference_bodymass_matrix[i,microbes_index]=1#Rmake_preference_bodymass(species_traits$bodymass[microbes_index],species_traits$bodymass[i],parameters)
    }
    # predatory interactions
    for (i in carnivores_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      preference_bodymass_matrix[i,species_index]=Rmake_preference_bodymass(species_traits$bodymass[species_index],species_traits$bodymass[i],parameters)
    }
    # detritivore interactions
    for (i in detritivores_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      preference_bodymass_matrix[i,c(faeces_index,FOM_index)]=Rmake_preference_bodymass(detritus_traits$bodymass[c(1:n_faeces,FOM_index_detritus)],species_traits$bodymass[i],parameters)
      #preference_bodymass_matrix[i,FOM_index]=1 # no size for FOM
    }
    pos<-which(preference_bodymass_matrix<pref_threshold,arr.ind=TRUE) # elements of the matrix bellow a threshold
    preference_bodymass_matrix[pos]=0 # neglects the interactions bellow the threshold
    # for (i in species_index){
    #   preference_bodymass_matrix[i,i]=0 # neglects cannibalism
    # }
    return(preference_bodymass_matrix)
  })
}
# prey preference depending prey stoichiometry content (for animals and detritus)
Rmake_preference_CN<-function(prey_CN,pred_CN,pref_CN_sigma){
  return(exp(-((prey_CN - pred_CN) / pref_CN_sigma)^2))
}
# prey preference depending on prey stoichiometry content (for animals and detritus) ????????????? a actualiser car le C:N des détritus est une variable
Rmake_preference_CN_matrix<-function(species_traits,trophic_group_parameters,preference_bodymass_matrix,foodweb_inputs){
  with(foodweb_inputs,{
    preference_CN_matrix=matrix(0,n_tot_species,n_tot) # matrix with all trophic species and detritus
    pos<-which(preference_bodymass_matrix>0,arr.ind=TRUE) # nonnull element of the body mass preference matrix
    preference_bodymass_matrix[pos]=1 # only considers nonnegligible interactions
    denominator=0 # sum of the preference for the rescaling
    # predatory interactions
    for (i in c(microbivores_index,carnivores_index)){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      preference_CN_matrix[i,1:n_tot_species]=Rmake_preference_CN(species_traits$CN[1:n_tot_species],species_traits$CN[i],trophic_group_parameters$pref_CN_sigma[species_traits$trophic_group_index[i]]) # preference for prey depending on C:N
      preference_CN_matrix[i,1:n_tot_species]=preference_CN_matrix[i,1:n_tot_species]*preference_bodymass_matrix[i,1:n_tot_species] # only keeps non negligible interactions
      denominator=sum(preference_CN_matrix[i,1:n_tot_species])
      for (j in 1:n_tot_species){
        if(preference_CN_matrix[i,j]>0){
          preference_CN_matrix[i,j]=preference_CN_matrix[i,j]/denominator # rescaling
        }
      }
    }
    # detritivore interactions
    for (i in detritivores_index){
      preference_CN_matrix[i,c(faeces_index,FOM_index)]=1 # this is recalculated dynamically at each time step
    }
    return(preference_CN_matrix)
  })
}
# effective attack rate resulting from basal attack rate, body mass and C:N ratios
Rmake_attack_rate<-function(species_traits,detritus_traits,trophic_group_parameters,foodweb_inputs){
  attack_rate_base_matrix=Rmake_attack_rate_base_matrix(species_traits,trophic_group_parameters,foodweb_inputs) # basal attack rate (prey-independent)
  preference_bodymass_matrix=Rmake_preference_bodymass_matrix(species_traits,detritus_traits,trophic_group_parameters,foodweb_inputs) # preference for prey depending on body mass ratio
  preference_CN_matrix=Rmake_preference_CN_matrix(species_traits,trophic_group_parameters,preference_bodymass_matrix,foodweb_inputs) # preference for prey depending on C:N ratio
  attack_rate_matrix=attack_rate_base_matrix*preference_bodymass_matrix*preference_CN_matrix
  return(list(attack_rate_base_matrix=attack_rate_base_matrix,
              preference_bodymass_matrix=preference_bodymass_matrix,
              preference_CN_matrix=preference_CN_matrix,
              attack_rate_matrix=attack_rate_matrix))
}

# HANDLING TIME # ----
# handling time for a predator-prey couple
Rmake_handling_time<-function(bodymass,kB,temperature,parameters){
  with(parameters,{
    conversion=fresh_to_dry*dry_to_C
    return(a_h*h0*exp(-E_h/(kB*temperature))*(bodymass*conversion)^s_h) # a_h is a tuning parameter
  })
}
# matrix with all handling times of carnivores and decomposers
Rmake_handling_time_matrix<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    handling_time_matrix=matrix(0,n_tot_species,n_tot)
    # predatory interactions
    for (i in c(microbivores_index,carnivores_index)){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      handling_time_matrix[i,1:n_tot_species]=Rmake_handling_time(species_traits$bodymass[i],kB,temperature,parameters)
    }
    # detritivore interactions (faeces and FOM consumption)
    for (i in detritivores_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      handling_time_matrix[i,c(faeces_index,FOM_index)]=Rmake_handling_time(species_traits$bodymass[i],kB,temperature,parameters) # no dependency of detritus particle size in Ott et al. (2012)
    }
    return(handling_time_matrix)
  })
}

# MICROBE MAXIMAL RESOURCE UPTAKE # ----
# maximum nitrogen uptake by microbes
Rmake_N_max_uptake<-function(temperature,parameters){
  with(parameters,{
    phi_N_max=a_phi_N*phi0_N*exp(s_phi_N*temperature)
    return(phi_N_max)
  })
}
# maximum FOM uptake by microbes
Rmake_FOM_max_uptake<-function(temperature,parameters){
  with(parameters,{
    phi_FOM_max=a_phi_FOM*phi0_FOM*exp(s_phi_FOM*temperature)
    return(phi_FOM_max)
  })
}
# maximum SOM uptake by microbes
Rmake_SOM_max_uptake<-function(temperature,parameters){
  with(parameters,{
    phi_SOM_max=a_phi_SOM*phi0_SOM*exp(s_phi_SOM*temperature)
    return(phi_SOM_max)
  })
}
# maximum DOC uptake by microbes
Rmake_DOC_max_uptake<-function(temperature,parameters){
  with(parameters,{
    phi_DOC_max=a_phi_DOC*phi0_DOC*exp(s_phi_DOC*temperature)
    return(phi_DOC_max)
  })
}
# maximum uptake matrix
Rmake_max_uptake_matrix<-function(trophic_group_parameters,species_traits,foodweb_inputs){
  with(foodweb_inputs,{
    max_uptake_matrix=matrix(0,n_tot_species,n_tot_detritus+n_nutrients)
    for (i in microbes_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      max_uptake_matrix[i,1:n_faeces]=Rmake_FOM_max_uptake(temperature,parameters) # faeces
      max_uptake_matrix[i,FOM_index_detritus]=Rmake_FOM_max_uptake(temperature,parameters) # FOM
      max_uptake_matrix[i,SOM_index_detritus]=Rmake_SOM_max_uptake(temperature,parameters) # SOM
      max_uptake_matrix[i,DOC_index_detritus]=Rmake_DOC_max_uptake(temperature,parameters) # DOC
      max_uptake_matrix[i,n_tot_detritus+1]=Rmake_N_max_uptake(temperature,parameters) # mineral nitrogen
    }
    return(max_uptake_matrix)
  })
}

# MICROBE HALF-SATURATION # ----
# half-saturation of nitrogen uptake by microbes
Rmake_N_half_saturation<-function(temperature,parameters){
  with(parameters,{
    K_N=a_K_N*K0_N*exp(s_K_N*temperature)
    return(K_N)
  })
}
# half-saturation of FOM uptake by microbes
Rmake_FOM_half_saturation<-function(temperature,parameters){
  with(parameters,{
    K_FOM=a_K_FOM*K0_FOM*exp(s_K_FOM*temperature)
    return(K_FOM)
  })
}
# half-saturation of SOM uptake by microbes
Rmake_SOM_half_saturation<-function(temperature,parameters){
  with(parameters,{
    K_SOM=a_K_SOM*K0_SOM*exp(s_K_SOM*temperature)
    return(K_SOM)
  })
}
# half-saturation of DOC uptake by microbes
Rmake_DOC_half_saturation<-function(temperature,parameters){
  with(parameters,{
    K_DOC=a_K_DOC*K0_DOC*exp(s_K_DOC*temperature)
    return(K_DOC)
  })
}
# half-saturation matrix
Rmake_half_saturation_matrix<-function(trophic_group_parameters,species_traits,foodweb_inputs){
  with(foodweb_inputs,{
    half_saturation_matrix=matrix(0,n_tot_species,n_tot_detritus+n_nutrients)
    for (i in microbes_index){
      parameters=as.list(trophic_group_parameters[species_traits$trophic_group_index[i],])
      half_saturation_matrix[i,1:n_faeces]=Rmake_FOM_half_saturation(temperature,parameters) # faeces
      half_saturation_matrix[i,FOM_index_detritus]=Rmake_FOM_half_saturation(temperature,parameters) # FOM
      half_saturation_matrix[i,SOM_index_detritus]=Rmake_SOM_half_saturation(temperature,parameters) # SOM
      half_saturation_matrix[i,DOC_index_detritus]=Rmake_DOC_half_saturation(temperature,parameters) # DOC
      half_saturation_matrix[i,n_tot_detritus+1]=Rmake_N_half_saturation(temperature,parameters) # mineral nitrogen
    }
    return(half_saturation_matrix)
  })
}

# ASSIMILATION EFFICIENCY # ----
Rmake_assimilation_efficicency_matrix<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    assimilation_efficicency_matrix=matrix(0,n_tot_species,n_tot)
    for (i in species_index){
      for (j in species_index){
        assimilation_efficicency_matrix[i,j]=trophic_group_parameters$epsilon[species_traits$trophic_group_index[j]] # assimilation of species resources
      }
      assimilation_efficicency_matrix[i,c(detritus_index)]=epsilon_detritus # assimilation of detritus
    }
    return(assimilation_efficicency_matrix)
  })
}

# DORMANCY SWITCH RATE # ----
# dormancy switch rate
Rmake_dormancy_rate<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    dormancy_rate_vector=array(0,n_tot_species)
    for (i in species_index){
      dormancy_rate_vector[i]=trophic_group_parameters$Q[species_traits$trophic_group_index[i]]
    }
    return(dormancy_rate_vector)
  })
}
# asymmetry between dormancy and awakening
Rmake_dormancy_asymmetry<-function(species_traits,trophic_group_parameters,foodweb_inputs){
  with(foodweb_inputs,{
    dormancy_asymmetry_vector=array(0,n_tot_species)
    for (i in species_index){
      dormancy_asymmetry_vector[i]=trophic_group_parameters$Q_asy[species_traits$trophic_group_index[i]]
    }
    return(dormancy_asymmetry_vector)
  })
}

# LIGHT INTERACTION MATRIX # ----
# list containing arrays with the index of each resource for each species
# microbes consumed by detritivores are not considered here
Rmake_light_interaction_list<-function(attack_rate_matrix,foodweb_inputs){
  with(foodweb_inputs,{
    light_interaction_list=vector("list",n_tot_species)
    for (i in microbes_index){
      light_interaction_list[[i]]=c(detritus_index,nutrients_index) # microbes interact with all detritus and mineral nutrient pools
    }
    for (i in c(detritivores_index,microbivores_index,carnivores_index)){
      light_interaction_list[[i]]=which(attack_rate_matrix[i,]>0,arr.ind=TRUE) # index of resource significantly consumed by trophic species i
    }
    return(light_interaction_list)
  })
}
# lists with the consumers of each resource ?????????????????????? à bien vérifier avant d'utiliser !!!!!
Rmake_reverse_interaction_list<-function(attack_rate_matrix,foodweb_inputs){
  with(foodweb_inputs,{
    reverse_interaction_list=vector("list", n_tot)
    for (j in 1:n_tot){
      reverse_interaction_list[[j]]=which(attack_rate_matrix[,j]>0,arr.ind=TRUE) # index of consumers significantly consuming resource j
    }
    return(reverse_interaction_list)
  })
}

# SPECIES INFORMATION # ----
Rmake_species_information<-function(species_traits,foodweb_inputs){
  with(foodweb_inputs,{
    species_information<-as.data.frame(matrix(nrow=n_tot_species,ncol=6))
    names(species_information)=c("species","trophic_group","trophic_type","population_structure","bodymass","CN")
    species_information$species=c(1:n_tot_species)
    species_information$trophic_group=species_traits$trophic_group
    species_information$trophic_type=species_traits$trophic_type
    species_information$population_structure=species_traits$population_structure
    species_information$bodymass=species_traits$bodymass
    species_information$CN=species_traits$CN
    return(species_information)
  })
}

########################################################################
### FUNCTIONS OF DEMOGRAPGIC PROCESSES (computed at each time step)  ### ----
########################################################################
# REAL RATE # ----
# real biomass loss or gain over a time step deltaT
Rcomputation_real_rate<-function(deltaT,instantaneous_rate){
 return(1-exp(-deltaT*instantaneous_rate))
}

# DETRITUS SURFACE # ----
# fraction of the total surface of detritus hosting microbes for each class of detritus
Rcomputation_detritus_surface<-function(detritus_density,detritus_fraction_available,n_trophic_group){ # detritus_density and detritus_fraction_available must be arrays
  detritus_surface<-vector("list",n_trophic_group)
  for (i in which(!sapply(detritus_fraction_available, is.null))){
    detritus_surface[[i]]=detritus_density*detritus_fraction_available[[i]]
    detritus_surface[[i]]=detritus_surface[[i]]/sum(detritus_surface[[i]])
    detritus_surface[[i]][is.na(detritus_surface[[i]])]=0 # no surface of no detritus
  }
  return(detritus_surface)
}

# METABOLISM - RESPIRATION # ----
Rcomputation_respiration<-function(density_vector,active_vector,metabolic_rate_vector,deltaT,foodweb_inputs){
  with(foodweb_inputs,{
    respiration=array(0,n_tot_species)
    for (i in species_index){
      respiration[i]=active_vector[i]*density_vector[i]*Rcomputation_real_rate(deltaT,metabolic_rate_vector[i])
    }
    return(respiration)
  })
}

# CONSUMPTION # ----
# prey preference depending on detritus stoichiometry for animal decomposers
Rmake_preference_CN_detritus_vector<-function(density_vector,CN_vector,CN_detritivore,light_interaction,pref_CN_sigma,pref_threshold){
  #CN_detritivore=CN_vector[i] # for test
  #light_interaction=light_interaction_list[[i]] # for test
  #CN_sigma=trophic_group_parameters$CN_sigma[species_traits$trophic_group_index[i]] # for test
  
  n=length(light_interaction)
  preference_CN_vector=array(0,n) # vector with the preference for each class of consumed detritus according to stoichiometry
  denominator=0 # for rescaling
  for (j in 1:n){
    if (density_vector[light_interaction[j]]>0){
      preference_CN_vector[j]=Rmake_preference_CN(CN_vector[light_interaction[j]],CN_detritivore,pref_CN_sigma) # preference for prey depending on C:N
      denominator=denominator+preference_CN_vector[j]
    }
  }
  if (denominator>0){
    for (j in 1:n){
      preference_CN_vector[j]=preference_CN_vector[j]/denominator # rescaling
      if (preference_CN_vector[j]<pref_threshold){
        preference_CN_vector[j]=0
      }
    }
  }else{
    preference_CN_vector=array(0,n)
  }
  return(preference_CN_vector)
}
# consumption of prey by predators and total consumption of prey
Rcomputation_instant_consumption_matrix<-function(density_vector,
                                                  CN_vector,
                                                  active_vector,
                                                  interference_rate_vector,
                                                  attack_rate_matrix,
                                                  handling_time_matrix,
                                                  light_interaction_list,
                                                  max_uptake_matrix,
                                                  half_saturation_matrix,
                                                  detritus_fraction_available,
                                                  species_traits,
                                                  trophic_group_parameters,
                                                  foodweb_inputs){
  with(foodweb_inputs,{
    instant_consumption_matrix=matrix(0,n_tot_species,n_tot) # instantaneous consumption of each resource by each consumer
    sum_consumption_array=array(0,n_tot) # total consumption of each resource
    denominator=0 # temporary variable used to compute functional responses
    # microbes # ----
    resource_density<-array(0,n_trophic_group) # total quantity of detritus (for FOM and faeces uptake)
    for (i in which(!sapply(detritus_fraction_available, is.null))){
      for (j in 1:FOM_index_detritus){ # faeces + FOM (shared functional response)
        resource_density[i]=resource_density[i] + density_vector[n_tot_species+j]*detritus_fraction_available[[i]][j]
      }
    }
    for (i in microbes_index){
      if (density_vector[i]>0){ # check if species i is not extinct
        # faeces + FOM (shared functional response)
        for (j in 1:FOM_index_detritus){ # faeces + FOM (shared functional response)
          instant_consumption_matrix[i,n_tot_species+j]=Rcomputation_uptake_instant(resource_density[species_traits$trophic_group_index[i]], # resource density
                                                                                    density_vector[i]*detritus_fraction_available[[species_traits$trophic_group_index[i]]][j], # consumer density and available fraction of detritus
                                                                                    active_vector[i],
                                                                                    max_uptake_matrix[i,j],
                                                                                    half_saturation_matrix[i,j])
          sum_consumption_array[n_tot_species+j]=sum_consumption_array[n_tot_species+j]+instant_consumption_matrix[i,n_tot_species+j] # total consumption of prey (sum of F_ij)
        }
        # SOM
        instant_consumption_matrix[i,SOM_index]=Rcomputation_uptake_instant(density_vector[SOM_index]*detritus_fraction_available[[species_traits$trophic_group_index[i]]][SOM_index_detritus], # resource density
                                                                            density_vector[i]*detritus_fraction_available[[species_traits$trophic_group_index[i]]][SOM_index_detritus], # consumer density
                                                                            active_vector[i],
                                                                            max_uptake_matrix[i,SOM_index_detritus],
                                                                            half_saturation_matrix[i,SOM_index_detritus])
        sum_consumption_array[SOM_index]=sum_consumption_array[SOM_index]+instant_consumption_matrix[i,SOM_index] # total consumption of prey (sum of F_ij)
        # DOC + nutrients (dissolved and fully available)
        for (j in DOC_index_detritus:(n_tot_detritus+n_nutrients)){ # DOC + mineral N (dissolved and fully available)
          instant_consumption_matrix[i,n_tot_species+j]=Rcomputation_uptake_instant(density_vector[n_tot_species+j], # resource density
                                                                                    density_vector[i], # consumer density
                                                                                    active_vector[i],
                                                                                    max_uptake_matrix[i,j],
                                                                                    half_saturation_matrix[i,j])
          sum_consumption_array[n_tot_species+j]=sum_consumption_array[n_tot_species+j]+instant_consumption_matrix[i,n_tot_species+j] # total consumption of prey (sum of F_ij)
        }
      }
    }
    # detritivores # ----
    for (i in detritivores_index){
      if (density_vector[i]>0){ # check if species i is not extinct
        denominator=1
        # interference
        denominator=denominator+interference_rate_vector[i]*density_vector[i]
        # update of the preference for detritus (detritus C:N changes at each time step)
        preference_CN_vector=Rmake_preference_CN_detritus_vector(density_vector,CN_vector,CN_vector[i],light_interaction_list[[i]],trophic_group_parameters$pref_CN_sigma[species_traits$trophic_group_index[i]],pref_threshold)
        k=0 # counter for preference_CN_vector
        for (j in light_interaction_list[[i]]){ # index of the resources consumed by the predator
          k=k+1
          denominator=denominator+(density_vector[j]*preference_CN_vector[k]*attack_rate_matrix[i,j]*handling_time_matrix[i,j])
        }
        k=0 # resetting of k
        # computation of F_i,j(t) and sum_predation_matrix [j] = sum_k F_kj
        for (j in light_interaction_list[[i]]){ # index of the resources consumed by the predator
          k=k+1
          instant_consumption_matrix[i,j]=active_vector[i]*density_vector[i]*preference_CN_vector[k]*attack_rate_matrix[i,j]/denominator # instantaneous consumption
          sum_consumption_array[j]=sum_consumption_array[j]+instant_consumption_matrix[i,j] # total consumption of prey (sum of F_ij)
        }
      }
    }
    # carnivores and microbivores # ----
    for (i in c(microbivores_index,carnivores_index)){
      if (density_vector[i]>0){ # check if species i is not extinct
        denominator=1
        # interference
        denominator=denominator+interference_rate_vector[i]*density_vector[i]
        for (j in light_interaction_list[[i]]){ # index of the resources consumed by the predator
          denominator=denominator+((active_vector[j]*density_vector[j])^(1+trophic_group_parameters$hill[species_traits$trophic_group_index[i]])*attack_rate_matrix[i,j]*handling_time_matrix[i,j])
        }
        # computation of F_i,j(t) and sum_predation_matrix [j] = sum_k F_kj
        for (j in light_interaction_list[[i]]){ # index of the resources consumed by the predator
          instant_consumption_matrix[i,j]=active_vector[i]*density_vector[i]*(active_vector[j]*density_vector[j])^(trophic_group_parameters$hill[species_traits$trophic_group_index[i]])*attack_rate_matrix[i,j]/denominator # instantaneous consumption
          sum_consumption_array[j]=sum_consumption_array[j]+instant_consumption_matrix[i,j] # total consumption of prey (sum of F_ij)
        }
      }
    }
    # output # ----
    return(list(instant_consumption_matrix=instant_consumption_matrix,
                sum_consumption_array=sum_consumption_array))
  })
}
# resource consumed
Rcomputation_consumption_matrix<-function(density_vector,
                                          CN_vector,
                                          active_vector,
                                          instant_consumption_matrix,
                                          sum_consumption_array,
                                          light_interaction_list,
                                          respiration,
                                          detritus_surface,
                                          deltaT,
                                          species_traits,
                                          foodweb_inputs){
  with(foodweb_inputs,{
    consumed_prey_density=array(0,n_tot) # quantity of consumed resource over a time step deltaT
    consumption_matrix=matrix(0,n_tot_species,n_tot) # consumption of each resource by each consumer over a time step deltaT
    ### general consumption (needs to be updated for particular trophic groups) ###
    # organisms
    for (j in 1:n_tot_species){
      consumed_prey_density[j]=active_vector[j]*density_vector[j]*Rcomputation_real_rate(deltaT,sum_consumption_array[j]) # loss of resource j due to consumption over a time step deltaT
    }
    # detritus
    for (j in (n_tot_species+1):n_tot){
      consumed_prey_density[j]=density_vector[j]*Rcomputation_real_rate(deltaT,sum_consumption_array[j]) # loss of resource j due to consumption over a time step deltaT
    }
    for (i in species_index){
      if (density_vector[i]>0){
        for (j in light_interaction_list[[i]]){ # index of resources consumed by trophic species i
          if (sum_consumption_array[j]>0){
            consumption_matrix[i,j]=consumed_prey_density[j]*instant_consumption_matrix[i,j]/sum_consumption_array[j] # quantity of resource j eaten by consumer i
          }
        }
      }
    }
    ### microbes ###
    growth=array(0,n_tot_species) # growth after stoichiometric homeostasis
    limitation=array("",n_tot_species) # qualitative limitation (C or N limitation)
    limitation_index=array(0,n_tot_species) # quantitative limitation (how far is the species from the TER)
    C_uptake_lim=0 # maximum C uptake from detritus pools
    N_uptake_lim=0 # maximum N uptake from mineral nitrogen
    ND_uptake_lim=0 # N uptake from detritus
    CN_detritus_average=0 # average C:N of consumed detritus
    ratio=0 # reduction coefficient of detritus consumption
    deltaC=0 # C budget (C intake minus respiration)
    detlaN=0 # N budget (N from minerals and detritus)
    for (i in microbes_index){
      if (density_vector[i]>0){
        # determination of limitation
        C_uptake_lim=0
        ND_uptake_lim=0
        for (j in detritus_index){
          if (consumption_matrix[i,j]>0){
            C_uptake_lim=C_uptake_lim+consumption_matrix[i,j]
            ND_uptake_lim=ND_uptake_lim+consumption_matrix[i,j]/CN_vector[j]
          }
        }
        CN_detritus_average=Rcomputaion_CN_detritus_average(C_uptake_lim,ND_uptake_lim)
        N_uptake_lim=consumption_matrix[i,N_index]
        deltaC=C_uptake_lim-respiration[i] # C budget
        detlaN=N_uptake_lim+ND_uptake_lim # N budget
        stoichio=Rcomputation_limitation(deltaC,detlaN,CN_vector[i]) # computation of limitation
        growth[i]=stoichio$growth # net biomass production
        limitation[i]=stoichio$limitation # qualitative limitation (C or N limitation)
        limitation_index[i]=stoichio$limitation_index # quantitative limitation (how far is the species from the TER)
        # recalculation of C and N uptakes (update of the consumption matrix)
        if (limitation[i]=="C-limitation"){
          consumption_matrix[i,N_index]=Rcomputation_uptake_nutrient_C_lim(C_uptake_lim,respiration[i],CN_detritus_average,CN_vector[i])
        }else if (limitation[i]=="N-limitation"){
          ratio=Rcomputation_uptake_detritus_N_lim(N_uptake_lim,respiration[i],CN_detritus_average,CN_vector[i])/C_uptake_lim # reduction coefficient of detritus consumption
          consumption_matrix[i,detritus_index]=ratio*consumption_matrix[i,detritus_index]
        }
        # in case of co-limitation, consumption_matrix remains unchanged
      }
    }
    ### microbes consumed by detritivores ### 
    for (i in detritivores_index){
      if (density_vector[i]>0){
        detritus_fraction<-array(0,n_trophic_group) # fraction of eaten detritus
        for (j in which(!sapply(detritus_surface, is.null))){
          detritus_fraction[j]=sum(detritus_surface[[j]]*consumption_matrix[i,detritus_index]/density_vector[detritus_index],na.rm=TRUE) # total "surface" of detritus eaten by 
        }
        for (j in microbes_index){
          consumption_matrix[i,j]=epsilon_microbivory*density_vector[j]*detritus_fraction[species_traits$trophic_group_index[j]] # microbial biomass eaten with detritus
        }
      }
    }
    ### updates the vector of consumed resource biomass for microbes, detritus and nutrients
    for (j in c(microbes_index,detritus_index,nutrients_index)){
      consumed_prey_density[j]=sum(consumption_matrix[,j]) # loss of resource j due to consumption over a time step deltaT
    }
    ### N mineralised by microbes
    N_mineralised=0
    for (i in microbes_index){
      if (consumption_matrix[i,N_index]<0){
        N_mineralised = N_mineralised - consumption_matrix[i,N_index]
      }
    }
    return(list(consumption_matrix=consumption_matrix,
                consumed_prey_density=consumed_prey_density,
                growth=growth, # for microbes only
                limitation=limitation, # for microbes only
                limitation_index=limitation_index, # for microbes only
                N_mineralised=N_mineralised)) # for microbes only
  })
}
  
# MICROBE STOICHIOMETRY # ----
# instantaneous uptake (Michaelis-Menten functional response)
Rcomputation_uptake_instant<-function(resource_density,consumer_density,active,phi_max,K){
  return(phi_max*active*consumer_density/(K+resource_density))
}
# average C:N of ingested detritus
Rcomputaion_CN_detritus_average<-function(C_uptake_lim,ND_uptake_lim){
  return(C_uptake_lim/ND_uptake_lim)
}
# detritus uptake in the case of N-limitation
Rcomputation_uptake_detritus_N_lim<-function(N_uptake_lim,respiration,CN_detritus_average,species_CN_average){
  return(CN_detritus_average*(species_CN_average*N_uptake_lim+respiration)/(CN_detritus_average-species_CN_average))
}
# nitrogen uptake in the case of C-limitation
Rcomputation_uptake_nutrient_C_lim<-function(C_uptake_lim,respiration,CN_detritus_average,species_CN_average){
  return((C_uptake_lim*(CN_detritus_average-species_CN_average)-CN_detritus_average*respiration)/(CN_detritus_average*species_CN_average))
}

# STOICHIOMETRY # ----
# limitation (C, N or co-limitation)
Rcomputation_limitation<-function(deltaC,deltaN,species_CN){
  growth_lim=c(deltaC,species_CN*deltaN) # maximal growth from C and N perspectives
  growth=min(growth_lim) # Liebig's law
  growth_max=max(growth_lim)
  limitation=which(growth_lim==growth) # identification of limitation
  limitation_index=(growth_lim[1]-growth_lim[2])/growth_max # quantification of limitation (<0 C limitation ; >0 N limitation)
  if (growth_max<1e-8){
    limitation_index=0
    limitation="C-limitation" # only C intake can be equal to zero or being negative
  }else if (length(limitation)==2){
    limitation="co-limitation"
  }else if(limitation==1){
    limitation="C-limitation"
  }else if(limitation==2){
    limitation="N-limitation"
  }
  return(list(growth=growth,
              limitation=limitation,
              limitation_index=limitation_index))
}
# assimilated, unassimilated and excreted C and N (stoichiometric budget)
Rcomputation_budget<-function(density_vector,
                              CN_vector,
                              active_vector,
                              species_traits,
                              light_interaction_list,
                              consumption_matrix,
                              respiration,
                              dormancy_rate_vector,
                              dormancy_asymmetry_vector,
                              assimilation_efficicency_matrix,
                              deltaT,
                              foodweb_inputs){
  with(foodweb_inputs,{
    growth=array(0,n_tot_species) # growth after stoichiometric homeostasis
    limitation=array("",n_tot_species) # qualitative limitation (C or N limitation)
    limitation_index=array(0,n_tot_species) # quantitative limitation (how far is the species from the TER)
    deltaC=0 # carbon budget
    deltaN=0 # nitrogen budget
    mineral_N_waste=array(0,n_tot_species) # mineral nitrogen excreted by each species
    mineral_C_waste=array(0,n_tot_species) # DOC excreted by each species
    detritus_N_waste=array(0,n_tot_detritus) # nitrogen excreted in detritus
    detritus_C_waste=array(0,n_tot_detritus) # carbon excreted in detritus
    diet_matrix=matrix(0,n_tot_species,n_tot) # biomass of each resource that is assimilated
    # animal species
    for (i in animals_index){
      if (density_vector[i]>0){
        deltaC=0
        deltaN=0
        pos_f=species_traits$faeces_pos[i] # associated faeces pool index (in 1:n_faeces)
        for (j in light_interaction_list[[i]]){
          diet_matrix[i,j]=assimilation_efficicency_matrix[i,j]*consumption_matrix[i,j] # assimilated biomass
          if (diet_matrix[i,j]>0){
            deltaC=deltaC+diet_matrix[i,j] # total C intake
            deltaN=deltaN+diet_matrix[i,j]/CN_vector[j] # total N intake
            if (j%in%microbes_index){ # microbivores detritus excretion goes in the SOM compartment
              detritus_C_waste[SOM_index_detritus]=detritus_C_waste[SOM_index_detritus]+(1-assimilation_efficicency_matrix[i,j])*consumption_matrix[i,j] # unassimilated C
              detritus_N_waste[SOM_index_detritus]=detritus_N_waste[SOM_index_detritus]+(1-assimilation_efficicency_matrix[i,j])*consumption_matrix[i,j]/CN_vector[j] # unassimilated N
            }else{ # production of faeces
              detritus_C_waste[pos_f]=detritus_C_waste[pos_f]+(1-assimilation_efficicency_matrix[i,j])*consumption_matrix[i,j] # unassimilated C
              detritus_N_waste[pos_f]=detritus_N_waste[pos_f]+(1-assimilation_efficicency_matrix[i,j])*consumption_matrix[i,j]/CN_vector[j] # unassimilated N
            }
          }
        }
        deltaC=deltaC-respiration[i] # net carbon budget
        stoichio=Rcomputation_limitation(deltaC,deltaN,CN_vector[i])
        growth[i]=stoichio$growth # growth after stoichiometric homeostasis
        limitation[i]=stoichio$limitation # qualitative limitation (C or N limitation)
        limitation_index[i]=stoichio$limitation_index # quantitative limitation (how far is the species from the TER)
        mineral_C_waste[i]=(deltaC-growth[i])*species_traits$recycling[i] # direct C recycling
        mineral_N_waste[i]=(deltaN-growth[i]/CN_vector[i])*species_traits$recycling[i] # direct N recycling
        detritus_C_waste[pos_f]=detritus_C_waste[pos_f]+(deltaC-growth[i])*(1-species_traits$recycling[i]) # indirect C recycling
        detritus_N_waste[pos_f]=detritus_N_waste[pos_f]+(deltaN-growth[i]/CN_vector[i])*(1-species_traits$recycling[i]) # indirect N recycling
      }
    }
    # microbes (calculation already done in the consumption matrix)
    C_uptake=0 # total C uptake
    q_act=0 # activation function
    switch=array(0,n_tot_species)
    for (i in microbes_index){
      if (density_vector[i]>0){
        # activity of microbes
        diet_matrix[i,detritus_index]=consumption_matrix[i,detritus_index] # diet of microbes
        C_uptake=sum(consumption_matrix[i,detritus_index]) # total carbon intake of microbe species i
        q_act=Rcomputation_activation(C_uptake,respiration[i]) # resource-dependent factor modulating the activation of dormant microbes
        switch[i]=Rcomputation_switch(q_act,active_vector[i],dormancy_rate_vector[i],dormancy_asymmetry_vector[i],deltaT)# variation of the fraction of active microbes over a time step deltaT
        # mineralisation of nitrogen
        if(consumption_matrix[i,n_tot_species+n_tot_detritus]<0){
          mineral_N_waste[i]=consumption_matrix[i,n_tot_species+n_tot_detritus] # nitrogen mineralisation by microbes
        }
      }
    }
    return(list(growth=growth, # for animals only
                limitation=limitation, # for animals only
                limitation_index=limitation_index, # for animals only
                mineral_C_waste=mineral_C_waste,
                mineral_N_waste=mineral_N_waste,
                detritus_C_waste=detritus_C_waste,
                detritus_N_waste=detritus_N_waste,
                diet_matrix=diet_matrix,
                switch=switch))
  })
}

# GROWTH RATE AGE STRUCTURED POPULATIONS # ----
# update the growth vector to include the age-structured trophic groups (AFTER ALL OTHER DEMOGRAPHIC PROCESSES)
Rcomputation_growth_age_structured<-function(density_vector,new_density_vector,growth,species_traits,n_mass_list,n_CN_list,age_structure_list,age_factors_vector,foodweb_inputs){
  with(foodweb_inputs,{
    #growth=budget$growth # FOR TEST
    growth_age=array(0,n_tot_species) # new vector of growth
    growth_biomass=0 # biomass transiting from class i to class i+1
    # age_structure_list: matrix with the position of trophic species in the density vector (row=size, column=C:N)
    for (k in age_structured_index){
      for (j in 1:n_CN_list[k]){ # C:N in column (each column corresponds to a species)
        for (i in 1:(n_mass_list[k]-1)){ # body mass in row (each row corresponds to an age class)
          if (growth[age_structure_list[[k]][i,j]]>0){
            growth_biomass=min(1,age_factors_vector[k]*growth[age_structure_list[[k]][i,j]]/(growth[age_structure_list[[k]][i,j]]+density_vector[age_structure_list[[k]][i,j]]))* # fraction of biomass transiting from class i to class i+1
              new_density_vector[age_structure_list[[k]][i,j]] # biomass after calculating all other demographic processes
            growth_age[age_structure_list[[k]][i,j]] = growth_age[age_structure_list[[k]][i,j]] - growth_biomass # loss due to the growth of individuals to the next class
            growth_age[age_structure_list[[k]][i+1,j]] = growth_age[age_structure_list[[k]][i+1,j]] + growth_biomass # gain due to the growth of individuals from the previous class
          }
        }
        if (growth[age_structure_list[[k]][n_mass_list[k],j]]>0){
          # first age class getting off-springs from adults
          growth_age[age_structure_list[[k]][1,j]] = growth_age[age_structure_list[[k]][1,j]] + growth[age_structure_list[[k]][n_mass_list[k],j]] # off-springs from adult reproduction
          growth_age[age_structure_list[[k]][n_mass_list[k],j]] = growth_age[age_structure_list[[k]][n_mass_list[k],j]] - growth[age_structure_list[[k]][n_mass_list[k],j]] # all the net growth is invested to reproduction
        }
      }
    }
    return(growth_age)
  })
}

# MICROBE DORMANCY # ----
# resource-dependent factor modulating the activation of dormant microbes
Rcomputation_activation<-function(C_uptake,respiration){
  q_act=C_uptake/(C_uptake+respiration)
  if (is.na(q_act)){
    q_act=0
  }
  return(q_act)
}
# variation of the fraction of active microbes over a time step deltaT
Rcomputation_switch<-function(q_act,active,dormancy_rate,dormancy_asymmetry,deltaT){
  return((dormancy_asymmetry*q_act/(dormancy_asymmetry*q_act + 1-q_act)-active)*Rcomputation_real_rate(deltaT,dormancy_rate*(dormancy_asymmetry*q_act + 1-q_act)))
  #return((q_act-active)*Rcomputation_real_rate(deltaT,dormancy_rate))
}

# MORTALITY AND DENSITY-DENPENDENT REGULATION # ----
# mortality due to natural mortality and self-regulation
Rcomputation_mortality<-function(density_vector,
                                 active_vector,
                                 mortality_rate_vector,
                                 self_regulation_rate_vector,
                                 shared_self_regulation_rate_vector,
                                 trophic_group_index_list,
                                 species_traits,
                                 deltaT,
                                 foodweb_inputs){
  with(foodweb_inputs,{
    mortality=array(0,n_tot_species)
    mortality_instant=0
    sum_density_trophic_group=array(0,n_trophic_group) # vector with the sum of biomass of each trophic group
    for (i in 1:n_trophic_group){
      sum_density_trophic_group[i]=sum(density_vector[trophic_group_index_list[[i]]])
    }
    for (i in species_index){
      mortality_instant=mortality_rate_vector[i] + # natural linear mortality
        self_regulation_rate_vector[i]*density_vector[i] + # independent self-regulation
        shared_self_regulation_rate_vector[i]*sum_density_trophic_group[species_traits$trophic_group_index[i]] # shared self-regulation
      mortality[i]=active_vector*density_vector[i]*Rcomputation_real_rate(deltaT,mortality_instant) # mortality over a time step deltaT
    }
    return(mortality)
  })
}
# carbon and nutrients released in FOM and SOM due to mortality
Rcomputation_mortality_detritus<-function(CV_vector,mortality,foodweb_inputs){
  with(foodweb_inputs,{
    FOM_C=0
    FOM_N=0
    SOM_C=0
    SOM_N=0
    for (i in microbes_index){
      SOM_C=SOM_C+mortality[i]
      SOM_N=SOM_N+mortality[i]/CV_vector[i]
    }
    for (i in animals_index){
      FOM_C=FOM_C+mortality[i]
      FOM_N=FOM_N+mortality[i]/CV_vector[i]
    }
    return(list(FOM_C=FOM_C,
                FOM_N=FOM_N,
                SOM_C=SOM_C,
                SOM_N=SOM_N))
  })
}

# BIOMASS DYNAMICS # ----
# dynamics of species biomass, nutrients and detritus
Rcomputation_foodweb_step_dynamics<-function(density_vector,
                                             CN_vector,
                                             active_vector,
                                             TL_vector,
                                             species_traits,
                                             detritus_traits,
                                             biological_rates,
                                             deltaT,
                                             trophic_group_parameters,
                                             foodweb_inputs){
  with(foodweb_inputs,{
    # computation of matter flows # ----
    respiration=Rcomputation_respiration(density_vector,active_vector,biological_rates$metabolic_rate_vector,deltaT,foodweb_inputs) # respiration over a time step deltaT
    detritus_surface=Rcomputation_detritus_surface(density_vector[detritus_index],biological_rates$detritus_fraction_available,n_trophic_group) # surface of detritus hosting microbes
    instant_consumption=Rcomputation_instant_consumption_matrix(density_vector, # instantaneous resource consumption
                                                                CN_vector,
                                                                active_vector,
                                                                biological_rates$interference_rate_vector,
                                                                biological_rates$attack_rate_matrix,
                                                                biological_rates$handling_time_matrix,
                                                                biological_rates$light_interaction_list,
                                                                biological_rates$max_uptake_matrix,
                                                                biological_rates$half_saturation_matrix,
                                                                biological_rates$detritus_fraction_available,
                                                                species_traits,
                                                                trophic_group_parameters,
                                                                foodweb_inputs)
    #instant_consumption_matrix=instant_consumption$instant_consumption_matrix # FOR TEST
    #sum_consumption_array=instant_consumption$sum_consumption_array # FOR TEST
    #fraction_available=detritus_traits$fraction_available # FOR TEST
    consumption=Rcomputation_consumption_matrix(density_vector, # resource consumption over a time step deltaT
                                                CN_vector,
                                                active_vector,
                                                instant_consumption$instant_consumption_matrix,
                                                instant_consumption$sum_consumption_array,
                                                biological_rates$light_interaction_list,
                                                respiration,
                                                detritus_surface,
                                                deltaT,
                                                species_traits,
                                                foodweb_inputs)
    #consumption_matrix=consumption$consumption_matrix # FOR TEST
    budget=Rcomputation_budget(density_vector, # C and N budget according to stoichiometric constrains, computation of recycling
                               CN_vector,
                               active_vector,
                               species_traits,
                               biological_rates$light_interaction_list,
                               consumption$consumption_matrix,
                               respiration,
                               biological_rates$dormancy_rate_vector,
                               biological_rates$dormancy_asymmetry_vector,
                               biological_rates$assimilation_efficicency_matrix,
                               deltaT,
                               foodweb_inputs)
    mortality=Rcomputation_mortality(density_vector,
                                     active_vector,
                                     biological_rates$mortality_rate_vector,
                                     biological_rates$self_regulation_rate_vector,
                                     biological_rates$shared_self_regulation_rate_vector,
                                     trophic_group_index_list,
                                     species_traits,
                                     deltaT,
                                     foodweb_inputs)
    mortality_detritus=Rcomputation_mortality_detritus(CN_vector,mortality,foodweb_inputs)
    
    # update the density of each compartment # ----
    # begin with all the loss of detritus (ORDER IS IMPORTANT)
    # general #
    new_density_vector=density_vector-consumption$consumed_prey_density # loss due to consumption/predation
    new_detritus_N=array(0,n_tot_detritus) # N gain by each detritus compartment
    # loss from leaching
    C_burrowing=leaching_SOM*density_vector[SOM_index]*deltaT
    new_density_vector[FOM_index]=new_density_vector[FOM_index]-leaching_FOM*density_vector[FOM_index]*deltaT # FOM leaching
    new_density_vector[SOM_index]=new_density_vector[SOM_index]-C_burrowing # SOM leaching
    new_density_vector[DOC_index]=new_density_vector[DOC_index]-leaching_DOC*density_vector[DOC_index]*deltaT # DOC  leaching
    new_density_vector[N_index]=new_density_vector[N_index]-leaching_N*density_vector[N_index]*deltaT # nitrogen leaching
    pos=0 # position in the new_detritus_N vector
    for (i in detritus_index){
      pos=pos+1
      if (new_density_vector[i]>0){
        new_detritus_N[pos]=new_density_vector[i]/CN_vector[i] # N in the remaining detritus before gain
      }
    }
    # species #
    new_density_vector[species_index]=new_density_vector[species_index]+
      consumption$growth + budget$growth - # net biomass production (after respiration)
      mortality # natural mortality and self-regulation
    growth_age_structured=Rcomputation_growth_age_structured(density_vector,
                                                             new_density_vector,
                                                             budget$growth, # growth of age structured populations
                                                             species_traits,
                                                             biological_rates$n_mass_list,
                                                             biological_rates$n_CN_list,
                                                             biological_rates$age_structure_list,
                                                             biological_rates$age_factors_vector,
                                                             foodweb_inputs)
    new_density_vector[species_index]=new_density_vector[species_index]+growth_age_structured # growth between age classes
    # detritus #
    new_density_vector[FOM_index]=new_density_vector[FOM_index]+mortality_detritus$FOM_C # detritus from natural mortality
    new_density_vector[SOM_index]=new_density_vector[SOM_index]+mortality_detritus$SOM_C # detritus from natural mortality
    new_detritus_N[FOM_index_detritus]=new_detritus_N[FOM_index_detritus]+mortality_detritus$FOM_N # detritus N from natural mortality
    new_detritus_N[SOM_index_detritus]=new_detritus_N[SOM_index_detritus]+mortality_detritus$SOM_N # detritus N from natural mortality
    new_density_vector[detritus_index]=new_density_vector[detritus_index]+budget$detritus_C_waste # detritus excreted
    new_density_vector[DOC_index]=new_density_vector[DOC_index]+sum(budget$mineral_C_waste) # DOC excreted
    new_detritus_N=new_detritus_N+budget$detritus_N_waste # N excreted in detritus
    # nutrients
    new_density_vector[N_index]=new_density_vector[N_index]+sum(budget$mineral_N_waste) # mineral nutrients excreted
    # gain from external inputs
    new_density_vector[FOM_index]=new_density_vector[FOM_index]+input_FOM*deltaT # FOM input
    new_density_vector[SOM_index]=new_density_vector[SOM_index]+input_SOM*deltaT # SOM input
    new_density_vector[DOC_index]=new_density_vector[DOC_index]+input_DOC*deltaT # DOC input
    new_density_vector[N_index]=new_density_vector[N_index]+input_N*deltaT # nitrogen input
    new_detritus_N[c(FOM_index_detritus,SOM_index_detritus)]=new_detritus_N[c(FOM_index_detritus,SOM_index_detritus)]+c(input_FOM/CN_ini_FOM,input_SOM/CN_ini_SOM)*deltaT
    # extinction #
    extinction<-which(new_density_vector[species_index]<extinction_threshold & new_density_vector[species_index]>0) # species getting extinct at this time step
    #print(extinction)
    new_density_vector[new_density_vector<extinction_threshold]=0 # no negative biomass allowed for biotic and abiotic compartments
    for (j in extinction){
      for (i in biological_rates$reverse_interaction_list[[j]]){
        # removes the extinct species from consumer diets and recompute the relative preference for C:N ratios
        biological_rates$light_interaction_list[[i]]<-biological_rates$light_interaction_list[[i]][-which(biological_rates$light_interaction_list[[i]]==j)]
        biological_rates$preference_bodymass_matrix[i,j]=0
        biological_rates$preference_CN_matrix=Rmake_preference_CN_matrix(species_traits,trophic_group_parameters,biological_rates$preference_bodymass_matrix,foodweb_inputs) # preference for prey depending on C:N ratio
        biological_rates$attack_rate_matrix=biological_rates$attack_rate_base_matrix*biological_rates$preference_bodymass_matrix*biological_rates$preference_CN_matrix
      }
    }
    
    # update the other dynamical vectors # ----
    # detritus C:N ratio
    pos=0 # position in the new_detritus_N vector
    new_CN_vector=CN_vector
    for (i in detritus_index){
      pos=pos+1
      if (new_detritus_N[pos]>0){
        new_CN_vector[i]=new_density_vector[i]/new_detritus_N[pos]
      }else{
        new_CN_vector[i]=0
      }
    }
    new_CN_vector[DOC_index]=Inf
    # fraction of active microbes
    new_active_vector=active_vector+budget$switch
    # limitation
    limitation=budget$limitation # limitation of animals
    limitation[microbes_index]=consumption$limitation[microbes_index] # limitation of microbes
    limitation_index=budget$limitation_index # quantitative limitation of animals
    limitation_index[microbes_index]=consumption$limitation_index[microbes_index] # quantitative limitation of microbes
    # conversion of the flows into rates
    # biomass of each resource that is consumed
    consumption_matrix=consumption$consumption_matrix/deltaT
    # biomass of each resource that is assimilated
    diet_matrix=budget$diet_matrix/deltaT
    # decomposition of detritus
    decomposition=consumption$consumed_prey_density[detritus_index]/deltaT
    # respiration
    respiration=respiration/deltaT
    # mineralisation of nitrogen
    N_mineralisation=(budget$mineral_N_waste + consumption$N_mineralised)/deltaT
    # burrowing of carbon in deep soil
    C_burrowing=C_burrowing/deltaT
    # detritus production rate
    detritus_production=array(0,4) # production rate of faeces, FOM, SOM and DOC by organisms
    detritus_production[1]=sum(budget$detritus_C_waste[faeces_index_detritus])/deltaT # faeces
    detritus_production[2]=(mortality_detritus$FOM_C + budget$detritus_C_waste[FOM_index_detritus])/deltaT # FOM
    detritus_production[3]=(mortality_detritus$SOM_C + budget$detritus_C_waste[SOM_index_detritus])/deltaT # SOM
    detritus_production[4]=sum(budget$mineral_C_waste)/deltaT # DOC

    return(list(density_vector=new_density_vector,
                CN_vector=new_CN_vector,
                active_vector=new_active_vector,
                limitation=limitation,
                limitation_index=limitation_index,
                consumption_matrix=consumption_matrix,
                diet_matrix=diet_matrix,
                decomposition=decomposition,
                respiration=respiration,
                N_mineralisation=N_mineralisation,
                C_burrowing=C_burrowing,
                detritus_production=detritus_production))
  })
}

#########################################
### VERIFICATIONS OF PARAMETER FILES  ### ----
#########################################
# species presence
Rcheck_species_presence<-function(trophic_group_parameters,n_trophic_group){
  for (i in 1:n_trophic_group){
    if (trophic_group_parameters$n_mass[i]>0 & trophic_group_parameters$n_CN[i]==0){
      stop(paste("Trophic group ",i," ",trophic_group_parameters$trophic_group[i],": C:N distribution is not defined",sep=""))
    }
    if (trophic_group_parameters$n_mass[i]==0 & trophic_group_parameters$n_CN[i]>0){
      stop(paste("Trophic group ",i," ",trophic_group_parameters$trophic_group[i],": body mass distribution is not defined",sep=""))
    }
  }
}
# age structured species
Rcheck_age_structured<-function(trophic_group_parameters,n_trophic_group){
  for (i in 1:n_trophic_group){
    if (trophic_group_parameters$n_mass[i]==1 & trophic_group_parameters$population_structure[i]=="age_structured"){
      stop(paste("Trophic group ",i," ",trophic_group_parameters$trophic_group[i],": more than 1 mass class is required",sep=""))
    }
  }
}
# integration time definition
Rcheck_time_steps<-function(T_simu,T_trans,T_step,n_step_output){
  if (T_trans>=T_simu){
    stop("T_trans must be smaller than T_simu")
  }
  if (floor(T_simu/T_step)<1){
    stop("No time steps are simulated, check the values of T_trans and deltaT")
  }
  # if (floor((T_simu-T_trans-deltaT)/deltaT)<n_step_output){
  #   stop("n_step_output is to large, check the values of T_trans, T_simu and n_step_output")
  # }
}

##########################
### OUTPUT FORMATTING  ### ----
##########################
# computation of trophic levels
Rcomputation_TL<-function(TL_vector,diet,foodweb_iputs){
  with(foodweb_iputs,{
    TL_vector_new<-TL_vector
    # TL microbes = 1
    # TL microbivores = 2
    denominator=0
    for (i in c(detritivores_index,carnivores_index)){
      TL_vector_new[i]=0 # initialisation
      denominator=0 # initialisation
      for (j in 1:n_tot){
        TL_vector_new[i]=TL_vector_new[i]+TL_vector[j]*diet[i,j]
        denominator=denominator+diet[i,j]
      }
      if (denominator<1e-20){
        TL_vector_new[i]=0
      }else{
        TL_vector_new[i]=1+TL_vector_new[i]/denominator
      }
    }
    return(TL_vector_new)
  })
}
# time series formatting for one simulation
Rmake_time_series<-function(output_simulation,n_step_output,name_output,names){
  n=length(output_simulation[[1]][[name_output]])
  if (n==0){
    return(NULL)
  }
  time_series=matrix(0,n_step_output,n+1)
  for (i in 1:n_step_output){
    #print(i)
    time_series[i,1]=output_simulation[[i]]$time
    time_series[i,1+(1:n)]=output_simulation[[i]][[name_output]]
  }
  colnames(time_series)=c("time",names)
  return(time_series)
}
# aggregate the data from all simulation into a single table (ALL SIMULATIOSN MUST HAVE THE SAME NUMBER OF COMPARTMENTS)
Rmake_data_tables<-function(results,sim_param,n_simu){
  # initalisation of the tables
  density=as.data.frame(matrix(0,n_simu,length(results[[1]]$density_mean)))
  CN=as.data.frame(matrix(0,n_simu,length(results[[1]]$CN_mean)))
  TL=as.data.frame(matrix(0,n_simu,length(results[[1]]$TL_vector)))
  active=as.data.frame(matrix(0,n_simu,length(results[[1]]$active_mean)))
  limitation_index=as.data.frame(matrix(0,n_simu,length(results[[1]]$limitation_index_mean)))
  decomposition=as.data.frame(matrix(0,n_simu,length(results[[1]]$decomposition_mean)))
  respiration=as.data.frame(matrix(0,n_simu,length(results[[1]]$respiration_mean)))
  N_mineralisation=as.data.frame(matrix(0,n_simu,length(results[[1]]$N_mineralisation_mean)))
  C_burrowing=as.data.frame(matrix(0,n_simu,1))
  detritus_production=as.data.frame(matrix(0,n_simu,4))
  # complete the tables
  for (i in 1:n_simu){
    density[i,]=results[[i]]$density_mean
    CN[i,]=results[[i]]$CN_mean
    TL[i,]=results[[i]]$TL_vector
    active[i,]=results[[i]]$active_mean
    limitation_index[i,]=results[[i]]$limitation_index_mean
    decomposition[i,]=results[[i]]$decomposition_mean
    respiration[i,]=results[[i]]$respiration_mean
    N_mineralisation[i,]=results[[i]]$N_mineralisation_mean
    C_burrowing[i,]=results[[i]]$C_burrowing_mean
    detritus_production[i,]=results[[i]]$detritus_production_mean
  }
  # adds the names of each variable
  colnames(density)=(results[[1]]$info)$names_tot
  colnames(CN)=(results[[1]]$info)$names_tot
  colnames(TL)=(results[[1]]$info)$names_species
  colnames(active)=(results[[1]]$info)$names_species
  colnames(limitation_index)=(results[[1]]$info)$names_species
  colnames(decomposition)=c((results[[1]]$info)$names_faeces,"FOM","SOM","DOC")
  colnames(respiration)=(results[[1]]$info)$names_species
  colnames(N_mineralisation)=(results[[1]]$info)$names_species
  colnames(C_burrowing)="C_burrowing"
  colnames(detritus_production)=c("faeces","FOM","SOM","DOC")
  # pastes the value of the variables
  density=cbind(sim_param,density)
  CN=cbind(sim_param,CN)
  TL=cbind(sim_param,TL)
  active=cbind(sim_param,active)
  limitation_index=cbind(sim_param,limitation_index)
  decomposition=cbind(sim_param,decomposition)
  respiration=cbind(sim_param,respiration)
  N_mineralisation=cbind(sim_param,N_mineralisation)
  C_burrowing=cbind(sim_param,C_burrowing)
  detritus_production=cbind(sim_param,detritus_production)
  
  return(list(density=density,
              CN=CN,
              TL=TL,
              active=active,
              limitation_index=limitation_index,
              decomposition=decomposition,
              respiration=respiration,
              N_mineralisation=N_mineralisation,
              C_burrowing=C_burrowing,
              detritus_production=detritus_production))
}
# save the outputs of simulations for each food web
Rsave_data_tables<-function(results,sim_param,n_simu,path_results){
  for (i in 1:n_simu){
    # consumption matrix
    consumption<-as.data.frame(results[[i]]$consumption_mean)
    colnames(consumption)<-results[[i]]$info$names_tot
    write.table(consumption,paste(path_results,"consumption_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    # diet matrix
    diet<-as.data.frame(results[[i]]$diet_mean)
    colnames(diet)<-results[[i]]$info$names_tot
    write.table(diet,paste(path_results,"diet_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    # save the time series data if recorded
    if (sim_param$TS_record[i]=="yes"){
      TS_names<-grep("TS_", names(results[[i]]), value = TRUE) # find the time series data
      index<-which(names(results[[i]])%in%TS_names) # index of the time series data in each row of results
      for (j in 1:length(TS_names)){
        write.table(results[[i]][[index[j]]],paste(path_results,TS_names[j],"_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
      }
    }
  }
  # save the tables summarising the results of all simulations
  for (i in unique(sim_param$trophic_group_file)){ # one table for each trophic group file
    sub<-which(sim_param$trophic_group_file==i) # subsets the simulations related to the trophic group
    tables<-Rmake_data_tables(results[sub], # generates the tables
                              sim_param[sub,],
                              length(sub))
    for (j in 1:length(tables)){
      write.table(tables[[j]],paste(path_results,names(tables)[j],"_",i,".txt",sep=""),sep=",",row.names=F,col.names=T) # writes all the results tables
    }
    # write additional files with miscellaneous information
    write.table(results[[sub[1]]]$species_traits,paste(path_results,"species_traits_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$detritus_traits,paste(path_results,"detritus_traits_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$info$names_species,paste(path_results,"names_species_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$info$names_faeces,paste(path_results,"names_faeces_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$info$names_detritus,paste(path_results,"names_detritus_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$info$names_nutrients,paste(path_results,"names_nutrients_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(results[[sub[1]]]$info$names_tot,paste(path_results,"names_tot_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
    write.table(as.data.frame(results[[sub[1]]]$info[1:7]),paste(path_results,"info_",i,".txt",sep=""),sep=",",row.names=F,col.names=T)
  }
  return("Output data saved")
}
# Make a table ready to use for ggplot
table_for_plot<-function(data,melt.var,id.vars,variable.name,value.name){
  # data=data_TL
  # melt.var=names_species
  # id.vars=""
  # variable.name="species"
  # value.name="TL"
  data<-data[,which(names(data)%in%c(melt.var,id.vars))] # only keeps the relevant columns
  data<-melt(data,
             id.vars = id.vars,
             variable.name = variable.name,
             value.name = value.name)
  return(data)
}

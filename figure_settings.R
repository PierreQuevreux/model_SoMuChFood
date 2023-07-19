################
### PACKAGES ### ----
################
library(ggplot2)
library(viridis)
library(scales)
library(ggnewscale)

#################
### FUNCTIONS ### ----
#################

`%notin%` <- Negate(`%in%`)

# reformat the data sets
Rmake_table_for_plot<-function(data,id.vars,variable.name,value.name){
  if (length(data)>0){
    data<-melt(data,
             id.vars = id.vars,
             variable.name = variable.name,
             value.name = value.name)
  }else{data=NULL}
  return(data)
}

# load and format the data sets
Rload_data<-function(path_results,trophic_group_file){
  # load the data sets
  info<-read.table(paste(path_results,"info_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  active<-read.table(paste(path_results,"active_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  C_burrowing<-read.table(paste(path_results,"C_burrowing_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  CN<-read.table(paste(path_results,"CN_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  decomposition<-read.table(paste(path_results,"decomposition_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  if (is.na(file.size(paste(path_results,"detritus_production_",trophic_group_file,".txt",sep="")))){
    detritus_production=NULL
  }else{
    detritus_production<-read.table(paste(path_results,"detritus_production_",trophic_group_file,".txt",sep=""),sep=",",header=T) 
  }
  density<-read.table(paste(path_results,"density_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  detritus_traits<-read.table(paste(path_results,"detritus_traits_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  limitation_index<-read.table(paste(path_results,"limitation_index_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  N_mineralisation<-read.table(paste(path_results,"N_mineralisation_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  respiration<-read.table(paste(path_results,"respiration_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  species_traits<-read.table(paste(path_results,"species_traits_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  TL<-read.table(paste(path_results,"TL_",trophic_group_file,".txt",sep=""),sep=",",header=T)
  n_tot=info$n_tot
  n_tot_species=info$n_tot_species
  n_faeces=info$n_faeces
  n_detritus_pool=info$n_detritus_pool
  n_tot_detritus=info$n_tot_detritus
  n_nutrients=info$n_nutrients
  n_parameters=info$n_parameters
  names_faeces<-read.table(paste(path_results,"names_faeces_",trophic_group_file,".txt",sep=""),sep=",",header=T)$x
  names_species<-read.table(paste(path_results,"names_species_",trophic_group_file,".txt",sep=""),sep=",",header=T)$x
  names_tot<-read.table(paste(path_results,"names_tot_",trophic_group_file,".txt",sep=""),sep=",",header=T)$x
  names_parameters=names(density)[1:n_parameters]
  names_detritus<-read.table(paste(path_results,"names_detritus_",trophic_group_file,".txt",sep=""),sep=",",header=T)$x
  
  # big data set summarising the data of each species
  ## density
  data_species<-Rmake_table_for_plot(density[,1:(n_parameters+n_tot_species)],names_parameters,"names","biomass")
  ## activity
  data_temp<-Rmake_table_for_plot(active,names_parameters,"names","active")
  data_species$active<-data_temp$active
  ## CN
  data_temp<-Rmake_table_for_plot(CN[,1:(n_parameters+n_tot_species)],names_parameters,"names","CN")
  data_species$CN<-data_temp$CN
  ## limitation_index
  data_temp<-Rmake_table_for_plot(limitation_index,names_parameters,"names","limitation_index")
  data_species$limitation_index<-data_temp$limitation_index
  ## N_mineralisation
  data_temp<-Rmake_table_for_plot(N_mineralisation,names_parameters,"names","N_mineralisation")
  data_species$N_mineralisation<-data_temp$N_mineralisation
  ## respiration
  data_temp<-Rmake_table_for_plot(respiration,names_parameters,"names","respiration")
  data_species$respiration<-data_temp$respiration
  ## activity
  data_temp<-Rmake_table_for_plot(TL,names_parameters,"names","TL")
  data_species$TL<-data_temp$TL
  # final touch
  levels(data_species$names)<-names_species
  data_species<-merge(data_species,species_traits[,-which(names(species_traits)%in%c("CN"))],by=c("names"))
  data_species$trophic_group<-factor(data_species$trophic_group, levels=unique(species_traits$trophic_group))
  data_species$names<-factor(data_species$names, levels=unique(species_traits$names))
  # species ID (position for the food web plots)
  data_species$num<-data_species$names
  levels(data_species$num)<-c(1:length(levels(data_species$num)))
  data_species$num<-as.numeric(data_species$num)
  
  # big data set summarising the data of each detritus
  ## density
  data_detritus<-Rmake_table_for_plot(density[,which(names(density)%in%c(names_parameters,names_detritus))],names_parameters,"names","biomass")
  ## CN
  data_temp<-Rmake_table_for_plot(CN[,which(names(density)%in%c(names_parameters,names_detritus))],names_parameters,"names","CN")
  data_detritus$CN<-data_temp$CN
  ## decomposition
  data_temp<-Rmake_table_for_plot(decomposition,names_parameters,"compartment","decomposition")
  data_detritus$decomposition<-data_temp$decomposition
  # final touch
  levels(data_detritus$names)<-names_detritus
  data_detritus<-merge(data_detritus,detritus_traits,by=c("names"))
  data_detritus$type<-factor(data_detritus$type,levels=c("faeces","FOM","SOM","DOC"))
  
  # format the data
  data_density<-Rmake_table_for_plot(density,names_parameters,"names","biomass")
  data_CN<-Rmake_table_for_plot(CN,names_parameters,"names","CN")
  data_detritus_production<-Rmake_table_for_plot(detritus_production,names_parameters,"type","production")
  levels(data_density$names)<-names_tot
  levels(data_CN$names)<-names_tot
  
  return(list(data_species=data_species,
              data_detritus=data_detritus,
              data_density=data_density,
              data_CN=data_CN,
              data_detritus_production=data_detritus_production,
              data_C_burrowing=C_burrowing,
              species_traits=species_traits,
              detritus_traits=detritus_traits,
              n_tot=n_tot,
              n_tot_species=n_tot_species,
              n_faeces=n_faeces,
              n_detritus_pool=n_detritus_pool,
              n_tot_detritus=n_tot_detritus,
              n_nutrients=n_nutrients,
              n_parameters=n_parameters,
              names_faeces=names_faeces,
              names_species=names_species,
              names_tot=names_tot,
              names_parameters=names_parameters,
              names_detritus=names_detritus))
}

# set the correct  number of species for each trophic group
Rcheck_count_species<-function(count){
  if (length(count)==0){
    return(0)
  }else{
    return(count)
  }
}

# create the colour vector corresponding to trophic groups
Rmake_colour_trophic_groupes<-function(species_traits,colour_trophic_groups_base,trophic_groups){
  count<-as.data.frame(table(species_traits$trophic_group))
  count$Var1<-factor(count$Var1, levels=unique(species_traits$trophic_group))
  count<-count[order(levels(count$Var1)),]
  colour_trophic_groups<-colour_trophic_groups_base[which(trophic_groups%in%count$Var1)]
  return(colour_trophic_groups)
}

# create the colour vector corresponding to trophic groups plus detritus classes
Rmake_colour_trophic_groupes_total<-function(species_traits,colour_trophic_groups_base,trophic_groups){
  count<-as.data.frame(table(species_traits$trophic_group))
  count$Var1<-factor(count$Var1, levels=unique(species_traits$trophic_group))
  count<-count[order(levels(count$Var1)),]
  colour_trophic_groups<-colour_trophic_groups_base[c(which(trophic_groups%in%count$Var1),1:4+length(trophic_groups))]
  return(colour_trophic_groups)
}

# create a colour palette for species
Rmake_colour_palette<-function(species_traits,n_faeces){
  colour_microbes<-colorRampPalette(c("hotpink1","hotpink4"))
  colour_bacteria<-colorRampPalette(c("hotpink1","hotpink4"))
  colour_fungi<-colorRampPalette(c("gold1","gold4"))
  colour_micro_foodweb_microbivores<-colorRampPalette(c("seagreen1","seagreen4"))
  colour_micro_foodweb_carnivores<-colorRampPalette(c("lightblue1","lightblue4"))
  colour_macro_foodweb_detritivores<-colorRampPalette(c("purple1","purple4"))
  colour_macro_foodweb_carnivores<-colorRampPalette(c("red1","red4"))
  colour_trophic_whales<-colorRampPalette(c("dodgerblue1","dodgerblue4"))
  colour_faeces<-colorRampPalette(c("sienna1","sienna4"))
  colour_FOM="orange1"
  colour_SOM="orange4"
  colour_DOC="black"
  colour_N="chartreuse4"
  count<-as.data.frame(table(species_traits$trophic_group))
  palette<-c(colour_microbes(Rcheck_count_species(count$Freq[count$Var1=="microbes"])),
             colour_bacteria(Rcheck_count_species(count$Freq[count$Var1=="bacteria"])),
             colour_fungi(Rcheck_count_species(count$Freq[count$Var1=="fungi"])),
             colour_micro_foodweb_microbivores(Rcheck_count_species(count$Freq[count$Var1=="micro-food web microbivores"])),
             colour_micro_foodweb_carnivores(Rcheck_count_species(count$Freq[count$Var1=="micro-food web carnivores"])),
             colour_macro_foodweb_detritivores(Rcheck_count_species(count$Freq[count$Var1=="macro-food web detritivores"])),
             colour_macro_foodweb_carnivores(Rcheck_count_species(count$Freq[count$Var1=="macro-food web carnivores"])),
             colour_trophic_whales(Rcheck_count_species(count$Freq[count$Var1=="trophic whales"])),
             colour_faeces(n_faeces),
             colour_FOM,
             colour_SOM,
             colour_DOC,
             colour_N)
  return(palette)
}

# line type palette
Rmake_linetype_palette<-function(n_species,n_faeces){
  return(c(rep("solid",n_species),rep("22",n_faeces),rep("solid",3),"22"))
}

# create the colour and line scales depending on the community
Rmake_palettes<-function(species_traits,n_tot_species,n_faeces){
  colour_group_microbes="hotpink2"
  colour_group_bacteria="hotpink2"
  colour_group_fungi="gold2"
  colour_group_micro_foodweb_microbivores="seagreen2"
  colour_group_micro_foodweb_carnivores="lightblue2"
  colour_group_macro_foodweb_detritivores="purple2"
  colour_group_macro_foodweb_carnivores="red2"
  colour_group_trophic_whales="dodgerblue2"
  colour_group_trophic_faeces<-"sienna2"
  colour_FOM="orange1"
  colour_SOM="orange4"
  colour_DOC="black"
  colour_N="chartreuse4"
  colour_trophic_groups_base=c(colour_group_microbes,
                               colour_group_bacteria,
                               colour_group_fungi,
                               colour_group_micro_foodweb_microbivores,
                               colour_group_micro_foodweb_carnivores,
                               colour_group_macro_foodweb_detritivores,
                               colour_group_macro_foodweb_carnivores,
                               colour_group_trophic_whales,
                               colour_group_trophic_faeces,
                               colour_FOM,
                               colour_SOM,
                               colour_DOC,
                               colour_N)
  trophic_groups=c("microbes",
                   "bacteria",
                   "fungi",
                   "micro-food web microbivores",
                   "micro-food web carnivores",
                   "macro-food web detritivores",
                   "macro-food web carnivores",
                   "trophic whales")
  colours_total=Rmake_colour_palette(species_traits,n_faeces)
  colours_species=Rmake_colour_palette(species_traits,0)
  colours_predators=Rmake_colour_palette(species_traits[species_traits$trophic_type%in%c("microbivores","carnivores"),],0)
  colours_detritivores=Rmake_colour_palette(species_traits[species_traits$trophic_type%in%c("microbes","detritivores"),],0)
  lines_total=Rmake_linetype_palette(n_tot_species,n_faeces)
  lines_species=Rmake_linetype_palette(n_tot_species,0)
  
  scale_colour_line_total<-scale_color_manual(values=colours_total,
                                              name=NULL)
  scale_colour_line_species<-scale_color_manual(values=colours_species,
                                                name=NULL)
  scale_colour_line_predator<-scale_color_manual(values=colours_predators,
                                                 name=NULL)
  scale_colour_line_detritivores<-scale_color_manual(values=colours_detritivores,
                                                     name=NULL)
  scale_colour_bar_total<-scale_fill_manual(values=colours_total,
                                            name=NULL)
  scale_colour_bar_species<-scale_fill_manual(values=colours_species,
                                              name=NULL)
  scale_linetype_total<-scale_linetype_manual(values=lines_total,
                                              name=NULL)
  scale_linetype_species<-scale_linetype_manual(values=lines_species,
                                                name=NULL)
  # colours for detritus only
  colours_detritus=Rmake_colour_palette(NULL,n_faeces)
  scale_colour_line_detritus<-scale_color_manual(values=colours_detritus,
                                                 name=NULL)
  scale_colour_bar_detritus<-scale_fill_manual(values=colours_detritus,
                                               name=NULL)
  # colours for aggregated trophic groups
  colour_trophic_groups<-Rmake_colour_trophic_groupes(species_traits,colour_trophic_groups_base,trophic_groups)
  scale_colour_line_trophic_groups<-scale_color_manual(values=colour_trophic_groups,
                                                       name=NULL,
                                                       guide = guide_legend(reverse=TRUE))
  scale_colour_bar_trophic_groups<-scale_fill_manual(values=colour_trophic_groups,
                                                     name=NULL,
                                                     guide = guide_legend(reverse=TRUE))
  colour_trophic_groups_detritivores<-colour_trophic_groups[which(colour_trophic_groups%in%c(colour_group_microbes,
                                                                                             colour_group_macro_foodweb_detritivores,
                                                                                             colour_group_trophic_whales))]
  scale_colour_line_trophic_groups_detritivores<-scale_color_manual(values=colour_trophic_groups_detritivores,
                                                                    name=NULL)
  scale_colour_bar_trophic_groups_detritivores<-scale_fill_manual(values=colour_trophic_groups_detritivores,
                                                                  name=NULL,
                                                                  guide = guide_legend(reverse=TRUE))
  colour_trophic_groups_consumers<-colour_trophic_groups[which(colour_trophic_groups%in%c(colour_group_micro_foodweb_microbivores,
                                                                                          colour_group_micro_foodweb_carnivores,
                                                                                          colour_group_macro_foodweb_carnivores))]
  scale_colour_line_trophic_groups_consumers<-scale_color_manual(values=colour_trophic_groups_consumers,
                                                                 name=NULL)
  
  # colours for aggregated trophic groups plus detritus classes
  colour_trophic_groups_total<-Rmake_colour_trophic_groupes_total(species_traits,colour_trophic_groups_base,trophic_groups)
  scale_colour_line_trophic_groups_total<-scale_color_manual(values=colour_trophic_groups_total,
                                                             name=NULL)
  scale_colour_bar_trophic_groups_total<-scale_fill_manual(values=colour_trophic_groups_total,
                                                           name=NULL)
  # colours for aggregated detritus classes
  colour_detritus_classes<-Rmake_colour_trophic_groupes_total(NULL,colour_trophic_groups_base,trophic_groups)
  scale_colour_line_detritus_classes<-scale_color_manual(values=colour_detritus_classes,
                                                         name=NULL)
  scale_colour_bar_detritus_classes<-scale_fill_manual(values=colour_detritus_classes,
                                                       name=NULL)
  scale_colour_bar_detritus_production<-scale_fill_manual(values=colour_detritus_classes[c(1,3)],
                                                          name=NULL)
  
  return(list(scale_colour_line_total=scale_colour_line_total,
              scale_colour_bar_total=scale_colour_bar_total,
              scale_colour_line_species=scale_colour_line_species,
              scale_colour_bar_species=scale_colour_bar_species,
              scale_colour_line_predator=scale_colour_line_predator,
              scale_colour_line_detritus=scale_colour_line_detritus,
              scale_colour_bar_detritus=scale_colour_bar_detritus,
              scale_colour_line_trophic_groups=scale_colour_line_trophic_groups,
              scale_colour_bar_trophic_groups=scale_colour_bar_trophic_groups,
              scale_colour_line_trophic_groups_detritivores=scale_colour_line_trophic_groups_detritivores,
              scale_colour_bar_trophic_groups_detritivores=scale_colour_bar_trophic_groups_detritivores,
              scale_colour_line_trophic_groups_consumers=scale_colour_line_trophic_groups_consumers,
              scale_colour_line_trophic_groups_total=scale_colour_line_trophic_groups_total,
              scale_colour_bar_trophic_groups_total=scale_colour_bar_trophic_groups_total,
              scale_colour_line_detritus_classes=scale_colour_line_detritus_classes,
              scale_colour_bar_detritus_classes=scale_colour_bar_detritus_classes,
              scale_colour_bar_detritus_production=scale_colour_bar_detritus_production,
              scale_linetype_total=scale_linetype_total,
              scale_linetype_species=scale_linetype_species))
}

# diplay the graphs if they are not saved
Rdisplay<-function(save,graph){
  if(!save){
    print(graph)
  }
}

######################
### GRAPH SETTINGS ### ----
######################
# align legend title # ----
align_legend <- function(p, hjust = 0.5)
{
  # extract legend
  g <- cowplot::plot_to_gtable(p)
  grobs <- g$grobs
  legend_index <- which(sapply(grobs, function(x) x$name) == "guide-box")
  legend <- grobs[[legend_index]]
  
  # extract guides table
  guides_index <- which(sapply(legend$grobs, function(x) x$name) == "layout")
  
  # there can be multiple guides within one legend box  
  for (gi in guides_index) {
    guides <- legend$grobs[[gi]]
    
    # add extra column for spacing
    # guides$width[5] is the extra spacing from the end of the legend text
    # to the end of the legend title. If we instead distribute it by `hjust:(1-hjust)` on
    # both sides, we get an aligned legend
    spacing <- guides$width[5]
    guides <- gtable::gtable_add_cols(guides, hjust*spacing, 1)
    guides$widths[6] <- (1-hjust)*spacing
    title_index <- guides$layout$name == "title"
    guides$layout$l[title_index] <- 2
    
    # reconstruct guides and write back
    legend$grobs[[gi]] <- guides
  }
  
  # reconstruct legend and write back
  g$grobs[[legend_index]] <- legend
  g
}
# themes # ----
theme<-theme_bw()+
  theme(panel.background = element_blank(),
        panel.border = element_blank(),
        text = element_text(family="serif",size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        axis.line.x.top = element_blank(),
        axis.line.y.right = element_blank(),
        strip.background = element_blank(),
        legend.key = element_blank(),
        plot.title = element_text(hjust = 0.5))

theme_raster<-theme_classic()+
  theme(panel.background = element_blank(),
        text = element_text(family="serif",size=20),
        axis.text = element_text(size=20),
        axis.line = element_line(),
        legend.key=element_blank(),
        legend.text.align=0,
        strip.background = element_blank(),
        plot.title = element_text(hjust = 0.5))

theme_vbar<-theme+
  theme(axis.title.x=element_blank(),
        panel.grid.minor.x=element_blank(),
        panel.grid.major.x=element_blank())

# axis # ----
x_axis_log10_short<-scale_x_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10_short<-scale_y_log10(labels = scales::trans_format("log10", scales::math_format(10^.x)))
x_axis_log10<-scale_x_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))
y_axis_log10<-scale_y_log10(breaks = trans_breaks("log10", function(x) 10^x),
                            labels = scales::trans_format("log10", scales::math_format(10^.x)))

# labels # ----
label_time="Time (days)"
label_time_short="Time (x100 days)"
label_biomass=expression("Biomass (mgC m"^-2*")")
label_biomass_potapov=expression(atop("Fresh biomass density","(mg gC"[soil]^-1*")"))
label_biomass_fraction="Fraction of the total biomass"
label_biomass_fraction_2lines="Fraction of the\ntotal biomass"
label_biomass_fraction_active="Fraction of the\ntotal active biomass"
label_bodymass_approx="Increasing body mass"
label_bodymass="Fresh body mass (mg)"
label_stock=expression("Stock (mgC m"^-2*")")
label_CN="C:N ratio"
label_active="Fraction of active biomass"
label_TL="Trophic level"
label_limitation_index="Limitation index"
label_decomposition=expression("Decomposition (mgC m"^-2*" d"^-1*")")
label_decomposition_rate=expression("Decomposition rate (d"^-1*")")
label_production=expression("Production (mgC m"^-2*" d"^-1*")")
label_respiration=expression("Respiration (mgC m"^-2*" d"^-1*")")
label_respiration_fraction=expression("Fraction of the total respiration")
label_respiration_fraction_2lines=expression("Fraction of the\ntotal respiration")
label_N_mineralisation=expression("N mineralisation (mgN m"^-2*" d"^-1*")")
label_diet="Predator diet composition"
label_DOC_input=expression("DOC input (mgC m"^-2*" d"^-1*")")
label_FOM_input=expression("FOM input (mgC m"^-2*" d"^-1*")")
label_c0_prey=expression(italic(c[0*",prey"]))
label_c0_pred=expression(italic(c[0*",pred"]))
label_multichannel="Multichannel"
label_size_structured="Size-stuctured"
label_detritivores_spectrum="Full size spectrum"
label_a_phi_DOC=expression("DOC decomposition factor "*italic(a[phi])[DOC])
label_a_phi_FOM=expression("FOM decomposition factor "*italic(a[phi])[FOM])
label_a_phi_SOM=expression("SOM decomposition factor "*italic(a[phi])[SOM])
label_a_phi_DOC_short=expression("DOC decomp. factor "*italic(a[phi])[DOC])
label_a_phi_FOM_short=expression("FOM decomp. factor "*italic(a[phi])[FOM])
label_a_phi_SOM_short=expression("SOM decomp. factor "*italic(a[phi])[SOM])
label_a_K_DOC=expression("DOC half-saturation factor "*italic(a[K])[DOC])
label_a_K_FOM=expression("FOM half-saturation factor "*italic(a[K])[FOM])
label_a_K_SOM=expression("SOM half-saturation factor "*italic(a[K])[SOM])
label_a_K_DOC_short=expression("DOC half-s. factor "*italic(a[K])[DOC])
label_a_K_FOM_short=expression("FOM half-s. factor "*italic(a[K])[FOM])
label_a_K_SOM_short=expression("SOM half-s. factor "*italic(a[K])[SOM])
label_Q_asy=expression("Asymmetry of dormancy rate "*italic(q[asy]))
label_radius_FOM=expression("Radius of FOM particles "*italic(r)[FOM]*" (\u03BCm)")
# palettes # ----
scale_colour_data_type<-scale_colour_manual(values=c("red","deepskyblue2"))

# model_SoMuChFood
Model used in the article "Consequences of the multichannel structure of soil food webs on ecosystem functioning" 

# Code of the model
main.R: code running parallelised simulations  
internal_functions_R.R: functions encoding the model  
parameter_generator.R: code generating the parameter tables feeding the model  
A folder named "results" must be created to save the output files.

# Input files
parameters_simulation.csv: general parameters  
parameters_trophic_groups_1.csv: parameters specific to trophic groups for the multichannel model  
parameters_trophic_groups_2.csv: parameters specific to trophic groups for the size spectrum model  
parameters_trophic_groups_3.csv: parameters specific to trophic groups for the size-structured model  
These files must be placed in a folder named "input_data".

# Statistical analysis and figures
parameter_calculation.R: calculation of parameters from data available online  
figure_analysis.R: ggplot code of the figures  
figure_settings.R: ggplot options used in figures (e.g. labels, colour scales and themes)

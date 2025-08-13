setwd('/Users/verjim/laptop_D_17.01.2022/Schmitz_lab/code/R/human_slice/hp_stats/')

# Slice incubation only
# won't save anything new - just warn that the files already exist
source("slice_incubation_only.R")

# AIS and MEA
source("AIS_staining_MEA.R")

# firing properties
source("firing_main_final_combined.R")

# intrinsic properties
source("intrinsic_main_final_combined.R")

# slice hrs_after_OP
source("slice_all_CTR_hrs_inc.R")

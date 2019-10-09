import numpy as np
import pandas as pd
import math
import glob
import sys
import os

sys.path.append('scripts/')
from analysis_trajectories import *

#################################
########### MAIN ################
#################################


nproc = 2
top_dir =  sys.argv[1] 
cluster = str(sys.argv[2])
analys_dir = './'
print(top_dir)

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')
print('out', out_dirs)

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence

XLs_cutoffs = {'BS3':30.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          analysis_dir = analys_dir,
                          nproc=nproc)

# Create dir
gsms_A_dir = analys_dir+'GSMs_cl'+cluster+'/sample_A'
gsms_B_dir = analys_dir+'GSMs_cl'+cluster+'/sample_B'

AT.create_gsms_dir(gsms_A_dir)
AT.create_gsms_dir(gsms_B_dir)

HA = AT.get_models_to_extract('selected_models_A_cluster-1_detailed.csv')
HB = AT.get_models_to_extract('selected_models_B_cluster-1_detailed.csv')
AT.do_extract_models(HA, 'h1', gsms_A_dir,'../output')
AT.do_extract_models(HB, 'h2', gsms_B_dir,'../output')

exit()



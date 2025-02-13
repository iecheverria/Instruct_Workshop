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

nproc = 1
top_dir =  sys.argv[1] 
analys_dir = './'

# Check if analysis dir exists
if not os.path.isdir(analys_dir):
    os.makedirs(analys_dir)

# How are the trajectories dir names
dir_head = 'run_'
out_dirs = glob.glob(top_dir+'/'+dir_head+'*/output/')

################################
# Get and organize fields for
# analysis
################################
# Read the total score, plot
# and check for score convengence
XLs_cutoffs = {'BS3':30.0}

# Load module
AT = AnalysisTrajectories(out_dirs,
                          dir_name = 'run_',
                          analysis_dir = analys_dir,
                          nproc=nproc,
                          nskip=10)

# Define restraints to analyze
AT.set_analyze_Connectivity_restraint()
AT.set_analyze_Excluded_volume_restraint()
AT.set_analyze_XLs_restraint(XLs_cutoffs = XLs_cutoffs)

# Read stat files
AT.read_stat_files()
AT.write_models_info()

AT.hdbscan_clustering(['EV_sum','CR_sum','XLs_sum'],
                      min_cluster_size=1000,
                      min_samples=1,
                      skip=10)

AT.summarize_XLs_info()

exit()



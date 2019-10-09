###############################
# Modeling of Rpb1/Rpb2
# dimer using XLs data
#
# Salilab - UCSF
###############################
import IMP
import RMF
import IMP.atom
import IMP.rmf
import IMP.pmi
import IMP.pmi.topology
import IMP.pmi.dof
import IMP.pmi.macros
import IMP.pmi.restraints
import IMP.pmi.restraints.stereochemistry
import IMP.pmi.restraints.crosslinking
from IMP.pmi.io.crosslink import CrossLinkDataBaseKeywordsConverter
import numpy as np
import sys

##############################
# Sampling paramemters
##############################   

num_frames = 50000
if '--test' in sys.argv: num_frames=100

##############################
# SYSTEM SETUP
##############################   

mdl = IMP.Model()

topo = 'topology_rnap.dat'
reader_sys = IMP.pmi.topology.TopologyReader(topo,
                                             pdb_dir = '../data/',
                                             fasta_dir = '../data/')

bs_sys = IMP.pmi.macros.BuildSystem(mdl,
                                      resolutions=[1])
bs_sys.add_state(reader_sys)

root_hier,  dof = bs_sys.execute_macro(max_rb_trans=2.0,
                                       max_rb_rot=0.2)
mols = bs_sys.get_molecules()[0]
##############################
# Connectivity
##############################
output_objects = [] # keep a list of functions that need to be reported
sample_objects = []
rmf_restraints = []

crs = []
for molname in mols:
    for mol in mols[molname]:
        copy_n = IMP.atom.Copy(mol.get_hierarchy()).get_copy_index()
        cr = IMP.pmi.restraints.stereochemistry.ConnectivityRestraint(mol)
        cr.set_label(mol.get_name()+'.'+str(copy_n))
        cr.add_to_model()
        output_objects.append(cr)
        crs.append(cr)

##############################
# Excluded Volume
##############################
evr1 = IMP.pmi.restraints.stereochemistry.ExcludedVolumeSphere(included_objects=mols.values(),
                                                               resolution=1)
evr1.add_to_model()
evr1.set_weight(1.0)
output_objects.append(evr1)


##############################
# XLs restraint
##############################    

cldbkc=IMP.pmi.io.crosslink.CrossLinkDataBaseKeywordsConverter()
cldbkc.set_protein1_key("Protein1")
cldbkc.set_protein2_key("Protein2")
cldbkc.set_residue1_key("AbsPos1")
cldbkc.set_residue2_key("AbsPos2")
cldbkc.set_unique_id_key("Id")
cldbkc.set_psi_key("Score")

# XLs RESTRAINT
cldb=IMP.pmi.io.crosslink.CrossLinkDataBase(cldbkc)
cldb.create_set_from_file("../data/data_xls_180215.dat")

xl1 = IMP.pmi.restraints.crosslinking.CrossLinkingMassSpectrometryRestraint(root_hier=root_hier,
                                                                            CrossLinkDataBase=cldb,
                                                                            resolution=1.0,
                                                                            length=21.0,
                                                                            slope=0.01)
xl1.add_to_model()
xl1.set_weight(1.0)

rmf_restraints.append(xl1)
output_objects.append(xl1)
dof.get_nuisances_from_restraint(xl1)

##############################    
# Write coordinates
##############################
out = IMP.pmi.output.Output()
out.init_rmf("all_ini.rmf3", [root_hier])
out.write_rmf("all_ini.rmf3")
out.close_rmf("all_ini.rmf3")

##############################
# Shuffle configurations
##############################    
IMP.pmi.tools.shuffle_configuration(root_hier,
                                    bounding_box=((-80,-80,-80),(80,80,80)))

dof.optimize_flexible_beads(500)

##############################
# MC sampling
##############################    

rex=IMP.pmi.macros.ReplicaExchange0(mdl,
                                    root_hier=root_hier,                          
                                    crosslink_restraints=rmf_restraints,          
                                    monte_carlo_sample_objects=dof.get_movers(),
                                    replica_exchange_minimum_temperature=1.0,
                                    replica_exchange_maximum_temperature=4.0,
                                    global_output_directory="output/",
                                    output_objects=output_objects,
                                    monte_carlo_steps=10,
                                    number_of_frames=num_frames,
                                    number_of_best_scoring_models=10)

rex.execute_macro()
exit()

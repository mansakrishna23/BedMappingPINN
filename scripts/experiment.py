# PINN experiment
# ========================================
# author: Mansa Krishna
# ========================================
# import relevant packages
from PINNICLE import pinnicle as pinn
import os
import numpy as np
from datetime import datetime
import deepxde as dde
import mat73
import json
import scipy.io as sio
# ========================================
# experiment definition
glacier = 'Deception' # glacier name: either Upernavik, Narssap, or Deception
random_seed = np.random.randint(1, 9999)
ptype = 'mc sb vel' # mc = mass conservation, sb = stress balance, vel = include log velocity
output_subfolder = glacier+'-MCSB-NBC'
#output_subfolder = 'tests'
iters = 700000 
# ========================================
dde.config.set_default_float('float64')
dde.config.disable_xla_jit()
dde.config.set_random_seed(random_seed)
# =======================================
pinnkeys = ptype.split(" ")
pname = ""
for i in range(len(pinnkeys)):
    pname += pinnkeys[i].upper()+"_"
modname = pname+glacier+'_'+str(random_seed)
print(modname)
# --------------------------------------
datestr = datetime.now().strftime("%Y%m%d_%H%M%S")
daystr = datetime.now().strftime("%d-%m-%Y")

rootpath = "./data/"+glacier
# path to input data and domain
issm_file = glacier.lower()+"_md.mat" # ISSM model loaded with all training data
radar_file = "ProcessedTracks.mat" # ice thickness data along flight tracks
exp_file = glacier+".exp" # PINN domain

# path to files
issm_path = os.path.join(rootpath, issm_file)
radar_path = os.path.join(rootpath, radar_file)
domain_path = os.path.join(rootpath, exp_file)
# path for saving models
output_filename = modname
model_folder = "./saved_models/"+output_subfolder+'/'+output_filename+"_"+datestr+"/"

# general parameters
hp = {}
hp["epochs"] = iters
hp["learning_rate"] = 0.0001
hp["loss_function"] = "MSE"
hp["save_path"] = model_folder
hp["is_save"] = True
hp["is_plot"] = True

# neural network
hp["activation"] = "tanh"
hp["initializer"] = "Glorot uniform"
hp["num_neurons"] = 128
hp["num_layers"] = 6

# data
print(pinnkeys)

issm = {}
issm["data_size"] = {'u':4000, 'v':4000}
issm["data_path"] = issm_path
if 'mc' in pinnkeys:
    issm["data_size"]["a"] = 4000
if 'sb' in pinnkeys:
    issm["data_size"]["s"] = 4000
if 'vel' in pinnkeys:
    issm["data_size"]["vel"] = 4000

mat = {}
mat["data_size"] = {"H":4000}
mat["data_path"] = radar_path
mat["name_map"] = {"H":"thickness"}
mat["source"] = "mat"

hp["data"] = {"ISSM":issm, "MAT":mat}

# domain
hp["shapefile"] = domain_path
hp["num_collocation_points"] = 9000

# additional loss functions
if 'vel' in pinnkeys:
    hp["additional_loss"] = {}

# vel log
vel_loss = {}
vel_loss['name'] = "vel MAPE"
vel_loss['function'] = "MAPE"
vel_loss['weight'] = 1.0e-5

if 'vel' in pinnkeys:
    hp["additional_loss"]["vel"] = vel_loss

# physics
MC = {}
MC['pde_weights'] = [1e12]

SSA = {}
SSA["pde_weights"] = [1e-12, 1e-12]

hp["equations"] = {}
if 'mc' in pinnkeys:
    hp["equations"]["MC"] = MC
if 'sb' in pinnkeys:
    hp["equations"]["SSA"] = SSA

# create experiment
md = pinn.PINN(hp)
print(md.params)
md.compile()

# train
md.train()


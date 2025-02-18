import os
import glob
import shutil
import h5py as h5
import numpy as np

#
# Revert MRPs to DM types
#

# Parameter file
fn = "DD0000/output_0000"

# Copy directory to save data (no overwriting)
origdir = fn.split('/')[0]
savdir = origdir + ".sav"
if not os.path.exists(savdir):
    shutil.copytree(origdir, savdir)

# Change all MRPs (type == 4) to DM (type == 1)
files = glob.glob(f"{fn}.cpu*")
for f in files:
    with h5.File(f, "a") as fp:
        for name, group in fp.items():
            if not name.startswith('Grid'): continue
            mrp = group['particle_type'][:] == 4
            group['particle_type'][mrp] = 1
            
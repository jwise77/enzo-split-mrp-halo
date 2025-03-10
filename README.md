# Enzo resimulation of a halo with particle splitting
_John H. Wise, 18 February 2025_

## Locating the particles to split

Depending on your scientific goals, select the halo of interest.  You will need its center and mass at the final output.  As an example, I have provided the notebook titled `FindResimHalo.ipynb`, showing me locating a halo "by eye" and identifying it in the halo finder list.  Given this halo, we want to determine its Lagrangian volume at the initial redshift z0, i.e. where this matter existed at this time.  The python code `get_halo_initial_extent.py` performs this operation and outputs an HDF5 file (starting with `initial_particle`) with the particle positions and IDs at z=z0.  To run this code, you will need to give

* The `my_halo` dictionary with
  * An ID (not used but for user convenience for bookkeeping)
  * Halo mass
  * Halo center
  * Note: `r_units` is not used
* Appropriate arguments in the `get_center_and_extent` routine
  * Initial output parameter filename
  * Final output parameter file
  * `round_size`, not used in the calculation but should equal the resolution of the base AMR level
  * `radius_factor`, used when defining the final Lagrangian volume

This routine will create a sphere centered on the given halo with a radius of `radius_factor` times the virial radius.  It will record the particle IDs within this sphere and then search for the same particles in the initial output.  It will output a file(s) in either plain text or HDF5 format with the initial particle positions and IDs, which Enzo can use to restrict AMR to its Lagrangian volume.

## Determining the particle splitting region

The Lagrangian region of a halo can be quite irregular.  Sometimes there are disconnected small volumes within this Lagrangian volume.  Usually these disconnected are small in fractional volume and can be ignored.  You will need to identify a bounding box for the Lagrangian volume.  Currently, Enzo only splits particles within a cube.  I have approached this in the notebook `Split-ICs.ipynb` by calculating some statistics and plotting some projections and histograms of the particles.  Here I have used the center of mass (all particles have equal mass) and manually determined the extent of the bounding box "by eye."  This Lagrangian volume center and width will be used in Enzo to split particles within this cube.

**Update (10 March 2025):** I have updated Enzo to consider non-cubic particle splitting regions.  I've made some associated changes to `Split-ICs.ipynb` to remove outliers (outside of 3-sigma) from the Lagrangian region and to make it a tight 3D bounding box.  Testing on the `SG64-L2` simulation, I see good results, applying AMR to only the halo progenitors.

## Calculating the Lyman-Werner (LW) intensity at the halo

* _Required:_ Star particle data in the directory `star-data/`
* _Required:_ Merger tree data readable by ytree.  Provided in `full_arbor.tgz`.

When resimulating a halo from a "full physics" simulation, we should include the effects of nearby stellar radiation.  We can calculate the LW intensity, using the star particle data, at the halo as a function of time and consider this intensity as a background in the resimulation.  This calculation is performed by the code `dcbh-lwb-file.py`, which assumes that the star metadata are stored in the directory `star-data` with one HDF5 file per simulation output.  It will also output plots of the halo merger tree, the LWB, contributions from star particles, and an input table (HDF5) for use with Grackle.

At the end of the code, you will need to give the halo unique ID from the halo catalog and the redshift range that you want the LWB calculated.  Usually, you'll want to select the final redshift and when the halo becomes resolved by the simulation.  The initial redshift can be determined by looking at the earliest redshift in the merger tree, which is produced by the code, and rerunning the code.

## Altering the simulation parameters

I recommend making a copy of the original data directory before modifying it.  The resimulation uses so-called "must refine particles" that restricts AMR to where these particles exist to increase computational efficiency, focusing our high resolution cells on the Lagrangian volume only.  If the previous simulation used must refine particles, we must reset their particle types to their original dark matter type so we follow a new Lagrangian volume.  I have provided a script called `reset-mrp.py` to perform this operation.  The SG64-L2 simulations requires this, but the Renaissance Simulations do not.

I have generated a diff file `changed_params.diff` between the original and modified parameter files.  The modified parameters descriptions are listed below.

* **StopFirstTimeAtLevel = 0 -> 14** :: This will stop the simulation when it reaches this AMR level, which will happen when a halo collapses.
* **grackle_data_file = filename** :: Use the HDF5 file produced by `dcbh-lwb-file.py` that provides the LW background.  The previous cooling table most likely provided additional background sources, i.e. Haardt-Madau background, metal cooling, etc.
* **UVbackground = 0 -> 1** :: Turns on the UV background defined in `grackle_data_file`
* **MetalCooling = 1 -> 0** :: Turns off metal cooling because we don't have any metals because there won't be any star formation/feedback in this simulation.
* **RadiativeTransfer = 1 -> 0** :: Turns off radiative transfer because there will be no stars
* **RefineByJeansLengthSafetyFactor = 4 -> 64** :: Resolves the local Jeans length by 64+ cells.
  * ⚠️ *NOTE: For testing purposes, leave this at 4 so that the simulations run faster* ⚠️
* **ParticleSplitterIterations = 0 -> 2** :: Each iteration will split a particle into 13 daughter particles.  Here will have 2 iterations, so 169 split particles.
* **ParticleSplitterCenter = x y z** :: Center of the bounding cube containing the particles that will be split
* **ParticleSplitterCenterRegion[dim] = width1 width2** :: Width of the bounding cube in the first and second iteration in dimension `dim`.  They must be of decreasing width because they are nested.  I usually make it ~10% wider with each iteration.
* **ParticleSplitterMustRefine = 0 -> 1** :: Enables marking must-refine particles from a file in the following parameter
* **ParticleSplitterMustRefineIDFile = filename** :: HDF5 filename with the positions and ID of the must refine particles. Use the file produced by `get_halo_initial_extent.py`
* **StarParticleCreation = 40 -> 0** :: Turn off star formation
* **StarParticleFeedback = 40 -> 0** :: Turn off stellar feedback
* **MustRefineParticlesRefineToLevel = 2 -> 3** :: To restrict the AMR grids to the volume encompassed by the must-refine particles, we must increase this parameter to one above the innermost "zoom-in" level

## Restarting

Now you're ready to restart the simulation from its initial output.  With these changes, Enzo will split the particles within the specified cube and mark a subset of those as must-refine particles.  Be sure that the two HDF5 files containing the UV/LW background and must refine particle IDs are located in the simulation top-directory (one level up from DD0000, for example).  The following command is an example with 8 cores, restarting from `DD0000/output_0000` and piping the stdout and stderr to the file `estd.out.0`.  If using an HPC system, change `mpirun` to the necessary MPI wrapper.
```
mpirun -n 8 /path/to/enzo.exe -d -r DD0000/output_0000 >& estd.out.0
```
The `-d` and `-r` flags respectively denote debug mode for extra (useful) output and a restart from an output.
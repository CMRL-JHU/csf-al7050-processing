# Cold Spray Formed AL7050 Pre/Post-processing #

### What is this repository for? ###

This repository is one part of two which together can create a 3D statistically equivalent virtual microstructure (SEVM) for cold spray formed (CSF) AL7050. The other repository is responsible for generating 3D bimodal masks representing the coarse grain (CG) and ultra-fine grain (UFG) regions via a modified version of [SliceGAN](https://github.com/stke9/SliceGAN) and can be found [here](https://github.com/CMRL-JHU/sliceGAN). This repository contains all the necessary scripts to:
1) Preprocess EBSD imagery of CSF AL7050 to obtain 2D bimodal masks as input for SliceGAN, and morphological and crystallographic characteristics of grains for 3D grain reconstruction in the SEVM.
2) Reconstruct the coarse grains within the CG region of the 3D bimodal mask from SliceGAN to complete the SEVM.

### How do I get set up? ###

* #### Dependencies ####
    | Name and Link | Version used for development
    |---------------|-----------------------------
    | Windows, Linux, MacOS | Ubuntu 22.04.3 LTS
    | [Dream3D](https://dream3d.bluequartz.net/Download/) 6 | Dream3D 6.5.141
    | [Python](https://www.python.org/downloads/) >= 3.10 and Pip or [Anaconda](https://www.anaconda.com/download) | Anaconda 24.1.2,  Python 3.10.14
    | [Matlab](https://www.mathworks.com/downloads/) | Matlab R2021a
    | [MTEX](https://mtex-toolbox.github.io/download) | MTEX 5.10.2
    Python packages and required versions contained in repository's requirements.txt file.
    
* #### Installation ####
    * Git clone the repository to a directory of your choosing. Repository zip is ~16 MB, so typical download time should be negligible.
    * Install Matlab as well as Python or Anaconda, whichever is preferred, according to their own instructions.
    * Download and uncompress Dream3D and MTEX, then copy the path to the MTEX root directory (e.g. ~/Downloads/mtex-5.10.2) to:
        * `grain_packing/pipeline_16_GrainSwapFast.m:path_mtex`
        * `grain_packing/pipeline_15_ReplaceCrystallography.m:path_mtex`
    * Install Python packages using the the included requirements file:
        * Using Anaconda:
            ```sh
            conda create -n sevm --file requirements.txt
            ```
        * Using Pip: 
            ```sh
            pip install -r requirements.txt
            ```
        * If the former did not work, it can also be manually created with Anaconda using the following commands:
            ```sh
            conda create -n sevm python=3.10
            conda install -n sevm h5py numpy scipy pandas seaborn opencv
            ```

### Contents ###

* *Preprocessing* - Contains scripts required to preprocess EBSD imagery of CSF AL7050 to obtain 2D bimodal masks as input for SliceGAN, and morphological and crystallographic characteristics of grains for 3D grain reconstruction in the SEVM.
    * *CG* - Contains scripts required to obtain characteristics of CGs.
        * *XY* - Contains scripts required to obtain prior particle scale information (volume fraction, bimodal mask) from low resolution EBSD scans in the plane comprised of the raster and deposition directions.
        * *YZ* - Contains scripts required to obtain prior particle scale information (volume fraction, bimodal mask) from low resolution EBSD scans in the plane comprised of the transverse and deposition directions.
        * *YZ Small* - Contains scripts required to obtain CG crystallographic and morphological data from high resolution EBSD scans in the plane comprised of the transverse and deposition directions.
    * *UFG* - Contains scripts required to opbtain characteristics of UFGs.
        * *YZ Small* - Contains scripts required to obtain UFG crystallographic and morphological data from high resolution EBSD scans in the plane comprised of the transverse and deposition directions.

* *Grain Packing* - Contains scripts required to reconstruct the coarse grains within the CG region of the 3D bimodal mask from SliceGAN to complete the SEVM.
   
### Running ###

If the installation succeeded, it should now be possible to run the scripts using the included example inputs.
* #### Preprocessing: ####
    Here we need only concern ourselves with the contents of the "CG" directory.
    The "xy" and "yz" directories contain scripts necessary to obtain information about the prior particles. In the "pipeline\_input" directory of each can be found example low resolution, high area EBSD scans.
    The "yz\_small" directory contains scripts necessary to obtain information about the CGs. In its "pipeline\_input" directory can be found an example high resolution, low area EBSD scan.
    The scripts in each of these directories will need to be run in numerical order. Explanations of, and directions for each script can be found in the first lines when opened in a code editor or in the first annotation filter in dream3d. *.py must be executed with Python, and *.json files must be executed with Dream3D. Scripts without a numerical identifier are support files and should not be run directly.
    The table below shows the intended inputs and the intended destination for the outputs of each directory.
    | Input                     | Directory   | Output     |
    |---------------------------|-------------|------------|
    | High Area (> 40 micrometer per side) EBSD file (ang, ctf, etc..) | CG/xy       | *.dream3d file containing 2D bimodal mask for SliceGAN       |
    | Orthogonal High Area (> 40 micrometer) EBSD file (ang, ctf, etc..)    | CG/yz       | *.dream3d file containing 2D bimodal mask for SliceGAN       |
    | High Resolution (< 1/3 micrometer/voxel edge) EBSD file (ang, ctf, etc..)         | CG/yz_small | *.dream3d file containing segmented grains for Processing |
    Typical runtime for each of the of the 7 scripts in "xy" and "yz" as well as the 6 scripts in "yz_small" should measure in single digit seconds using the provided example inputs.
* #### SliceGAN ####
    SliceGAN is not contained in this repository (see above) so it will not be covered in depth here, but below is a table of its intended inputs and the intended destination for its outputs.
    | Input                     | Output     |
    |---------------------------|------------|
    | *.dream3d file containing 2D bimodal mask from Preprocessing (both) | *.dream3d file containing 3D bimodal mask for Processing |
    Training time for SliceGAN should take approximately 16 hours an a single Nvidia Tesla A100 using the provided sample inputs. Generation time should measure in the double digit seconds.
* #### Grain Packing ####
    The grain packing directory contains all the scripts necessary to create a SEVM for CSF microstructures. It takes in the bimodal CG/UFG domain mask from SliceGAN and packs it with grains having statistically equivalent morphological and crystallographic characteristics.
    As in Preprocessing, the scripts in each of these directories will need to be run in numerical order. Explanations and instructions are once again included in at the top of every script.  *.py and *.pyw must be executed with Python, *.json files must be executed with Dream3D, and *.m files must be executed with matlab. Scripts without a numerical identifier are support files and should not be run directly.
    Note that if you would like to run a small example script, it would be very beneficial to sample a small subdomain with pipeline component #8:
    `python pipeline_08_ChooseSubdomain.py --x 0 20 --y 0 20 --z 0 20`
    Note also that due to the non-determinate nature of pipeline components #10, #16, and #17, no two SEVMs created by this pipeline will be exactly alike.
    The table below shows the intended inputs and the intended destination for the output.
    | Input                                 | Output                |
    |---------------------------------------|-----------------------|
    | *.dream3d file containing segmented grains from Preprocessing, *.dream3d file containing 3D bimodal mask from SliceGAN | *.dream3d file containing 3D SEVM for Meshing, Fem Software |
    Typical run time for each of the scripts varies on the size of the of the chosen subdomain. Each of the 20 scripts should measure in the single or double digit seconds provided the above [0, 20] voxel size was used. The voxel size of [0, 317] used in the accompanying manuscript can take upwards of several days on an Intel© Core™ i7-9750H CPU @ 2.60GHz × 6 with 16 GBs of non-ECC, Synchronous DDR4 @ 2667 MHz.

### Who do I talk to? ###

* Joshua Stickel: jsticke3@jh.edu
* Brayan Murgas: bmurgas2@jh.edu

### License ###

MIT
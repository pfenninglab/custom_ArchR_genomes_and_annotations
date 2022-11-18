INSTALLING ARCHR/SEURAT IN CONDA ON LANE  

This is a quick guide on how to set up ArchR on Lane. First, you should understand why these steps are necessary. Conda is a package management system that allows users to create and enter virtual environments. Virtual environments are designed to make working on projects simpler by mitigating many of the issues which commonly come up when trying to run code a cluttered environment. For example, if one of my projects requires r-base 4.1.3, but another requires r-base 4.2.0, then installing both on my local machine or cluster will cause a lot of headaches when trying to switch versions. How would the computer know which version to use? Which package library should it reference? Where to install new packages? These and a litany of other issues can cause problems when trying to install everything you need for all projects or packages into you native environment.

Conda allows you to create virtual environments with remarkable ease, taking care of all the conflicts, installations, and organization on its own. It also allows users to recreate identical Conda environemnts through the use of .yaml or .yml files. Conda installation instructions below are taken from: https://github.com/pfenninglab/mouse_sst#setup


INSTALLING CONDA
"Download the latest Miniconda installer. This is the correct installer for lane and bridges:

-cd /tmp
-curl -O https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh

If you're not on lane or bridges, check your system's architecture and download the correct installer from Latest Miniconda Installer Links.
Run the installer:

-bash Miniconda3-latest-Linux-x86_64.sh

You will need to answer a few questions when prompted. The default values for each question should work.
Cleanup and exit the interactive node:

-rm Miniconda3-latest-Linux-x86_64.sh
-exit

This will return you to the head node on the cluster."


CREATING A VIRUAL ENVIRONMENT FROM A .YAML FILE
Using an SCP transfer program add the files 'ArchR_Env.yaml' and  'install_ArchR_dependencies.R' to a directory on lane. Run the following commands.

Be warned! Conda has gone through many versions, so some functions are out of date. The proper code is NOT:
conda create --file ArchR_Env.yaml

Rather, it is:
-conda env create --file ArchR_Env.yaml

In this file are R packages that are difficult to install thorugh R, Python packages which ArchR relies upon i.e. macs2, the Jupyter Kernel so you can run Jupyter notebooks, etc. 
Once you have created the virtual environment from the ArchR_Env.yaml, you can run 

-conda list -n ArchR_Env

to look at what packages are installed in your virtual environment.


INSTALLING R PACKAGES WITHIN YOUR CONDA ENVIRONMENT

Ideally, we'd be able to install all the R packages we needed through our .yaml, and that would be the end of it. Indeed, we've already installed plenty of essential R packages like r-xml, r-ggplot2, r-devtools, etc. However, not all of the packages we want are available through Conda's channels. This is why we need to install them through R, which when done in our Conda environment, will create an R package library unique to that environment. You can see this for yourself by running library() in R when in your Conda environment's R kernel.

Before you run install_ArchR_dependencies.R, make sure you request a compute node. Instructions on how to do that can be found on the lane_cluster tutorial. 

Activate your newly created conda environment:
-conda activate ArchR_Env

Run:
-echo on
-Rscript install_ArchR_dependencies.R

Once this completes, congrats! You should be able to launch an R session and see ArchR, Seurat, and all the necessary dependencies installed into your environment. 


BONUS TIPS: JUPYTER NOTEBOOKS

The lane_cluster tutorial has instructions on how to install and run jupyter notebooks on lane. In ArchR_Env.yaml, the R kernel required for Jupyter is already provided, so all you need to do is switch the kernel after launching your notebook. To make accessible 'ArchR_Env' from your jupyter notebook, make sure to set 'source activate Archr_Env' in your jupyter_script.sbatch

In the top right of a running notebook, you should see something like Python[conda:env something something]. In the menu bar, go to KERNEL -> CHANGE KERNEL -> R[conda env:ArchR_Env]. Dont't worry about activating your conda environment in lane first, you will have access to all your different kernels and virutal environments in Jupyter. You should be good to go!

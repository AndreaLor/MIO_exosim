# ExoSim_N
 
ExoSim_N (Exoplanet Observation Simulator - New) is a generic time-domain simulator of exoplanet transit spectroscopy.
It is a modification of the original ExoSim, now written in Python 3 and re-structured for improved operability.

Installation
------
We recommend setting up a virtual environment  for ExoSim_N to avoid package conflicts.  The environment should install Python=3.8.5, matplotlib=3.3.1-0, setuptools=49.6.0, numpy=1.19.1, numba=0.50.1.  Using conda the following command line instruction will do all this for you:.

    conda create -n exosim python=3.8.5 matplotlib=3.3.1-0 setuptools=49.6.0 numpy=1.19.1 numba=0.50.1

Then activate this environment. Depending on the system the activation command may be any one of the following:

    source activate exosim
    
or    

    conda activate exosim
    
or    
    
    activate exosim


### GitHub

Next, clone the repository from github.

    git clone https://github.com/subisarkar/ExoSim_N.git


### Databases

Next, download the following databases.  

[Phoenix BT-Settl database](https://phoenix.ens-lyon.fr/Grids/BT-Settl/CIFIST2011_2015/FITS/BT-Settl_M-0.0a+0.0.tar) (Allard F., Homeier D., Freytag B., 2012, Philos. Trans. Royal Soc. A, 370, 2765).  
These are the stellar spectrum models.  Unzip the folder to give the folder labelled `BT-Settl_M-0.0a+0.0`.  (There is no need to unzip the enclosed fits files.)  Then moved the folder into  `ExoSim_N/archive/` .

[LDC](https://drive.google.com/file/d/1lWRdqW_wI3y31ugqq2HfyyekGyOSteL_/view?usp=sharing)  
These contain pre-calculated limb darkening coefficients obtained using ExoTETHyS (Morello, G. et al. (2020). AJ, 159,  75).  Unzip the folders and move them into `ExoSim_N/archive/` . 

From the [NASA Exoplanet Archive](https://exoplanetarchive.ipac.caltech.edu/cgi-bin/TblView/nph-tblView?app=ExoTbls&config=PS&constraint=default_flag=1), download the 'Planetary Systems' table as follows: under 'Download Table' choose 'CSV format', 'Download all columns', 'Download all rows', and then hit 'Download Table'.  Move the .csv file into the `ExoSim_N/archive/` .


### Set up
Next navigate to inside the `ExoSim_N` folder and run `setup.py`.  This will setup the remaining package dependencies for ExoSim-N.

    cd ExoSim_N
    python setup.py install


### Output folder
By default ExoSim_N will place the results from the simulations into the `ExoSim_N/output/` folder. However you may choose a different location to store these results.  To do so find the file `exosim_paths.txt` in the `ExoSim_N/exosim_n/input_files/` folder, and next to `output_directory` (leaving at least one white space gap),  enter the full path to the location including folder name, e.g. `/Users/UserA/Desktop/ExoSim_Results`.  The folder will be automatically generated. If this is blank, the default location will be chosen.

### Your system
ExoSim_N is a generic shell.  You will need a 'system' folder that contains all instrument-specific data to make it work.  The folder should be named with the name of the system: e.g. 'ProjectX', and placed in the `ExoSim_N/systems/` folder.  The folder must contain an ExoSim-compatible XML input configuration file with the name `ProjectX.XML`  (replace 'ProjectX with name of the system).  The folder also needs to have a folder named `PSF` in which specific PSFs are contained, if you wish to use pre-generated PSFs.  Alternately the code can use Airy functions to produce PSFs.  All transmission files, QE files, etc. should be contained in this system folder.    To obtain the system folders for Ariel and other missions please contact: subhajit.sarkar@astro.cf.ac.uk.

Running a simulation
------
Navigate to inside the `ExoSim_N` folder, and run the `run_exosim.py` file with an input parameter file (e.g. `exosim_input_params_ex1.txt`) as the argument.  Some example input parameter files are provided (see below).

      cd ExoSim_N
      python run_exosim.py exosim_input_params_ex1.txt
      
Alternately if using an IDE (e.g. Spyder), you can open the file `ExoSim_N/exosim)n/run_files/run_exosim.py` and run from within the environment using the parameter file name (e.g. `exosim_input_params_ex1.txt`) as the argument for the function `run`.  Results will be packaged as a .pickle file in the designated output folder.  The code will also display results once completed.

Results
------
The results from any given simulation are placed in the output folder.  To identify the simulation, the file name includes the type of simulation, the planet, instrument channel and end time of the simulation.  A file with the suffix 'TEMP' is a temporary file for an ongoing multi-realisation simulation, but can be accessed in the same way as a completed simulation file.  The data is in the form of a dictionary.   To display results from a results file (e.g. `xxxx.pickle`):

    cd ExoSim_N
    python results.py xxxx.pickle

Make sure the file is in the designated output directory for this to work.
In addition a .txt with the same name as the dictionary file will be generated containing derived simulation values (e.g. number of NDRs, integration time etc.) and the input parameter file values used.


Use of code
------

Please contact: subhajit.sarkar@astro.cf.ac.uk
 
Citing
------

If you use ExoSim in your research, please cite:
Sarkar, S.et al. 2020. ExoSim: the exoplanet observation simulator. Experimental Astronomy.

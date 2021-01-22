# ExoSim_N tutorial
 
We provide three example input files.  These can be used to run any system, however the default channels in these files correspond to Ariel channels. 

Example 1 : Allan deviation analysis
------
Navigate to inside the `ExoSim_N` folder. 

      cd ExoSim_N
      
Enter the following.
      
      python run_exosim.py exosim_n_input_params_ex1.txt
      
Alternately if using an IDE (e.g. Spyder), you can open the file `ExoSim_N/exosim)n/run_files/run_exosim.py` and run from within the environment using the parameter file name (e.g. `exosim_n_input_params_ex1.txt`) as the argument for the function `run`.  

This example will run an out-of-transit simulation, followed by Allan deviation analysis.  It will automatically run  `results.py` to display the final results.  The results can also be displayed by entering the following:
 
    python results.py xxxx.pickle

where  `xxxx.pickle`  is the output file name in the output directory.  The results will show the signal, noise, fractional noise at T14, and the predicted noise on the transit depth.


Example 2 : Full transit simulation with Monte Carlo method
------
Navigate to inside the `ExoSim_N` folder. 

cd ExoSim_N

Enter the following.

python run_exosim.py exosim_n_input_params_ex2.txt

Alternately if using an IDE (e.g. Spyder), you can open the file `ExoSim_N/exosim)n/run_files/run_exosim.py` and run from within the environment using the parameter file name (e.g. `exosim_n_input_params_ex2.txt`) as the argument for the function `run`.  

This example will run a Monte Carlo full transit simulation with 25 realizations.  It will automatically run  `results.py` to display the final results.  The results can also be displayed by entering the following:

python results.py xxxx.pickle

where  `xxxx.pickle`  is the output file name in the output directory.  The results will show the the predicted noise on the transit depth and spectra with error bars.



Example 3 : Noise budget
------
Navigate to inside the `ExoSim_N` folder. 

cd ExoSim_N

Enter the following.

python run_exosim.py exosim_n_input_params_ex3.txt

Alternately if using an IDE (e.g. Spyder), you can open the file `ExoSim_N/exosim)n/run_files/run_exosim.py` and run from within the environment using the parameter file name (e.g. `exosim_n_input_params_ex3.txt`) as the argument for the function `run`.  

This example will run an out of transit simulation, cycling over all noise sources.  It will automatically run  `results.py` to display the final results.  The results can also be displayed by entering the following:

python results.py xxxx.pickle

where  `xxxx.pickle`  is the output file name in the output directory.  The results will show the the signal and noise per spectral bin, per noise source.


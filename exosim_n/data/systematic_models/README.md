# Insert your systematic grids here

Each folder has a name, e.g. 'Systematics1'.  This path needs to added to the input parameters file if you are simulating the systematic.  In the folder, three files are needed in either .txt, or .npy format.  Name the wavelength file 'wl.txt' or 'wl.npy'.  Name the time file 'time.txt' or 'time.npy'.  Name the systematic file 'syst.txt', or 'syst.npy'.   Wavelength should be in microns, and time in seconds.  The systematic should be the fractional change on the baseline with time and wavelength, e.g. an increase of 10% is a value of 1.1.



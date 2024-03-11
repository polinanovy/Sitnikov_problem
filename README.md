# Computer solution for the Sitnikov's problem
Solves and plots the Poincare map for the Sitnikov's problem

The program implements [Sitnikov's problem](https://en.wikipedia.org/wiki/Sitnikov_problem) based on the theory from the workbook [Celestial Dynamics: Chaoticity and Dynamics of Celestial Systems](https://onlinelibrary.wiley.com/doi/book/10.1002/9783527651856)

## Files
- **INPUT.json** - file with input parameters
- **mod.py** - module-file with all needed functions for the main program
- **run.py** - main program
- **e_0.01_example.png**, **e_0.01_zoom_example.png**, **e_0.11_example.png** - example plots for different eccentricity and scales

To run the program, you can use the command ```python3 run.py```

## Required python packages
- json
- os
- subprocess
- numpy
- scipy
- matplotlib
- pathlib
- tqdm

You can install them using ```pip``` or ```pip3```

## INPUT-file includes
- eccentricity $0 \leqslant e < 1$
- initial velocity *v* for the falling body
- absolute tolerance *atol*. The solver keeps the local error estimates less than atol + rtol * abs(y)
- relative tolerance *rtol*
- name of the method for integration *method*
- the height *z1* to which the body falls
- initial height *z2* for the falling body
- *step* for heights from z1 to z2
- number of rotations *n_rot*
- name of the directory for saving the results *save_dir*. If you don't have this directory, the program will ask you to create it and do it for you



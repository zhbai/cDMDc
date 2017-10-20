Compressive DMDc (cDMDc)

Copyright 2017, All Rights Reserved
Code by Zhe Bai
For Paper, "Dynamic mode decomposition for compressive system identification"
by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton.

The algorithms are in the “utils” directory:
	
	func_DMD   —  (paper) Algorithm 1 Exact DMD.__
	func_DMDc  -  (paper) Algorithm 2 DMD with control.__
	func_cDMD  -  (paper) Algorithm 3 Compressive DMD.__	
	func_cDMDc -  (paper) Algorithm 4 Compressive DMD with control.__
	checkModes -  reorder the modes and return the index of the mode with zero eigenvalues.__
	normalize  -   normalize the DMD modes for calibration.__
	turn_sign  -   regulate the DMD modes for plots and error report.__
	cosamp 	   -   iterative signal recovery from incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.__
	
Code to generate the model, figures, and output error for the example of “stochastically forced linear system” is in matlab file “EXAMPLE” in the main directory. (Example 4.1 in the paper)
    - Reconstruction errors are reported in ./output/
    - Figures are generated in ./figures/

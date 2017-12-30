# Compressive Dynamic Mode Decomposition with Control (cDMDc)
cDMDc [1] is a novel framework for compressive system identification unifying two recent innovations that extend DMD to systems with actuation [2] and systems with heavily subsampled measurements [3]. Using cDMDc it is possible to identify a low-order model from limited input–output data and reconstruct the associated full-state dynamic modes with compressed sensing, adding interpretability to the state of the reduced-order model. Moreover, when full-state data is available, it is possible to dramatically accelerate downstream computations by first compressing the data. 

The publication "Dynamic mode decomposition for compressive system identification"
by Z. Bai, E, Kaiser, J. L. Proctor, J. N. Kutz and S. L. Brunton. is available on [arXiv](https://arxiv.org/abs/1710.07737).

## Installation

1. Clone this repository to your desktop.
2. Add path to `cDMDc/utils` folder to Matlab search path using `addpath('<path to mds>/cDMDc/utils')`.

## Dependencies
For the reconstruction of the modes, cDMDc uses the CoSaMP algorithm developed by Deanna Needell and Joel Tropp [4].

An implementation written by David Mary and Bob L. Sturm is added in the folder `cDMDc/utils` , which is released under the the Commons Creative Licence with Attribution Non-commercial Share Alike (by-nc-sa).


## Getting Started

See `Ex01_SFLS.m` for demonstrating the approach on a stochastically forced linear system with known low-rank dynamics and an artificially inflated state dimension (Example 4.1 in [1]). Just execute this file in MatLab and it will generate the model, figures, and compute the reconstruction error.

## Organization

The algorithms are in the `cDMDc/utils` directory:
	
	func_DMD   —  Exact DMD, Algorithm 1 in [1].
	func_DMDc  -  DMD with control, Algorithm 2 in [1].
	func_cDMD  -  Compressive DMD, Algorithm 3 in [1].	
	func_cDMDc -  Compressive DMD with control, Algorithm 4 in [1].
	checkModes -  reorder the modes and return the index of the mode with zero eigenvalues.
	normalize  -   normalize the DMD modes for calibration.
	turn_sign  -   regulate the DMD modes for plots and error report.
	cosamp 	   -   iterative signal recovery from incomplete and inaccurate samples" by Deanna Needell and Joel Tropp.
    
## License (MIT license)

See the [LICENSE file](LICENSE) for details.

## References

[1] Zhe Bai, Eurika Kaiser, Joshua L. Proctor, J. Nathan Kutz, and Steven L. Brunton. arXiv Preprint, 2017.__
[2] J. L. Proctor, S. L. Brunton, and J. N. Kutz. Dynamic mode decomposition with control. SIAM J. Applied Dynamical Systems, 15(1):142–161, 2016.__
[3] Steven L. Brunton, Joshua L. Proctor, Jonathan H. Tu, and J. Nathan Kutz. Compressed sensing and dynamic mode decomposition. Journal of Computational Dynamics, 2(2):165–191, 2015.__
[4] D. Needell and J. A. Tropp. CoSaMP: iterative signal recovery from incomplete and inaccurate samples. Communications of the ACM, 53(12):93–100, 2010.

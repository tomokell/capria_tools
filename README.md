# Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA) Tools
[![DOI](https://zenodo.org/badge/513092852.svg)](https://zenodo.org/badge/latestdoi/513092852)

Matlab code for CAPRIA signal simulation (e.g. for flip angle schedule optimisation) and image reconstruction.
## CAPRIA Signal Simulation (*models*)
Included here are some simple models for both the angiographic and perfusion signal, including effects from a train of excitation pulses with constant or variable flip angles (assuming a spoiled gradient echo was used). Currently, dispersion is ignored, but we plan to include this soon. Calculations use standard angiographic and perfusion signal models (referenced below) modified to account for the signal attenuation and changes in the amount of excited transverse magnetisation present with a variable flip angle readout. 

Angiographic or perfusion signals can be calculated using calls to `CAPRIASignalSimple.m`. Some examples are given in `ExampleCAPRIASignalSimulations.m`. 

## CAPRIA image reconstruction (*recon*)
This code allows a variety of techniques to be used to reconstruct CAPRIA angiograms and/or perfusion images from raw k-space data in Siemens twix (meas.dat) format. This includes:
- Extracting the relevant parts of the raw data
- Phase correction
- Construction of an "anatomical" image, using the mean label/control data
- Coil sensitivity estimation (with adaptive combine)
- Coil compression
- A variety of reconstruction options:
    - adjoint non-uniform fast Fourier transform (NUFFT), similar to regridding
    - regularised conjugate gradient SENSE (cgSENSE)
    - spatially sparse with L2 temporal smoothness
    - xf-sparse (i.e. sparse in the spatial-temporal frequency domain)
    - locally low rank (LLR)
    - spatial total variation (xTV)
    - xTV + L2 temporal smoothness
    - low rank plus sparse (L+S)

In our experience, cgSENSE, spatial sparsity with L2 temporal smoothness and LLR work most robustly on CAPRIA data, although a comprehensive comparison has not been performed.

The main call is to `CAPRIARecon.m`, which can accept various options to turn on or off sections of the reconstruction. Due to the relatively long reconstruction times, the code is structured to save out intermediate files (e.g. coil sensitivity maps), so subsequent calls will avoid re-running the same stages that have already produced results. Some example angiographic and perfusion-like reconstruction calls are given in `Example_CAPRIA_Angio_Recon.m` and `Example_CAPRIA_Perf_Recon.m`.

## Utilities (*utils*)
This folder contains various helper functions that are used by `models` and `recon` and should be included in the Matlab path. 

## External dependencies
Tested using MATLAB 2017b and 2019a. Some of the recon scripts may rely on the following toolboxes:
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

In addition, the recon scripts require:
- `mapVBVD` by Philipp Ehses and colleagues, for Siemens raw data reading: [[GitHub Link]](https://github.com/CIC-methods/FID-A/blob/master/inputOutput/mapVBVD)
- `NUFFT` operators within the Image Reconconstruction Toolbox (IRT) by Jeff Fessler and colleagues: [[Download here]](http://web.eecs.umich.edu/~fessler/irt) [[GitHub link]](https://github.com/JeffFessler/mirt)
- `FSL` Matlab tools and utilities for saving results to Nifti format: [[Download here]](https://fsl.fmrib.ox.ac.uk)

Note that we make heavy use of a subset of Mark Chiew's Matlab image reconstruction code, which is included here for simplicity, but the full set can be found on his [website](https://users.fmrib.ox.ac.uk/~mchiew/Tools.html) or on [GitHub](https://github.com/mchiew/recon-tools-matlab).

For some of the plotting, we have also included the handy utility `hatchfill2` and associated license, but the original version can be found [here](https://uk.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2).

## References
A manuscript describing an implementation of 4D CAPRIA and optimisation using variable flip angle schedules and an LLR reconstruction has been accepted for publication in Magnetic Resonance in Medicine:

- *Okell, T. W. & Chiew, M. Optimization of 4D Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling. Magnetic Resonance in Medicine. (2022) doi:10.1002/mrm.29558. [[DOI link](https://doi.org/10.1002/mrm.29558)]*

A 2D implementation of CAPRIA was previously described in:

- *Okell, T. W. Combined angiography and perfusion using radial imaging and arterial spin labeling. Magnetic Resonance in Medicine 81, 182–194 (2019).*

Code and data to reproduce the statistical analyses and figures in the preprint can be found here:

- *Okell, Thomas W & Chiew, Mark. 4D Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA): Data and Code to Reproduce Statistics and Figures. (Zenodo, 2022). doi:10.5281/zenodo.7390441.* [[Zenodo link](https://zenodo.org/record/7390441)]

The angiographic signal model uses a modified version of the kinetic model described in:

- *Okell, T. W., Chappell, M. A., Schulz, U. G. & Jezzard, P. A kinetic model for vessel-encoded dynamic angiography with arterial spin labeling. Magnetic Resonance in Medicine 68, 969–979 (2012).*

The base perfusion model is the (P)CASL model described in:

- *Buxton, R. B. et al. A general kinetic model for quantitative perfusion imaging with arterial spin labeling. Magnetic Resonance in Medicine 40, 383–396 (1998).*

The 3D extension of the golden ratio for MRI is described in:

- *Chan, R. W., Ramsay, E. A., Cunningham, C. H. & Plewes, D. B. Temporal stability of adaptive 3D radial MRI using multidimensional golden means. Magn Reson Med 61, 354–363 (2009).*

The phase correction step is described in:

- *Schauman, S. S., Chiew, M. & Okell, T. W. Highly accelerated vessel‐selective arterial spin labeling angiography using sparsity and smoothness constraints. Magnetic Resonance in Medicine 83, 892–905 (2020).*

Coil sensitivity estimation uses the adaptive combine algorithm:

- *Walsh, D. O., Gmitro, A. F. & Marcellin, M. W. Adaptive reconstruction of phased array MR imagery. Magnetic resonance in medicine 43, 682–90 (2000).*

Coil compression is described in:

- *Buehrer, M., Pruessmann, K. P., Boesiger, P. & Kozerke, S. Array compression for MRI with large coil arrays. Magn. Reson. Med. 57, 1131–1139 (2007).*

Conjugate gradient SENSE is described in:

- *Pruessmann, K. P., Weiger, M., Börnert, P. & Boesiger, P. Advances in sensitivity encoding with arbitrary k -space trajectories: SENSE With Arbitrary k -Space Trajectories. Magn. Reson. Med. 46, 638–651 (2001).*

xf-sparse and xTV are forms of compressed sensing described in:

- *Lustig, M., Donoho, D. & Pauly, J. M. Sparse MRI: The application of compressed sensing for rapid MR imaging. Magnetic Resonance in Medicine 58, 1182–1195 (2007).*

The xTV algorithms use Constrained Fast Gradient Projection:

- *Beck, A. & Teboulle, M. Fast Gradient-Based Algorithms for Constrained Total Variation Image Denoising and Deblurring Problems. IEEE Trans. on Image Process. 18, 2419–2434 (2009).*

x-sparse algorithms use FISTA:

- *Beck, A. & Teboulle, M. A Fast Iterative Shrinkage-Thresholding Algorithm for Linear Inverse Problems. SIAM Journal on Imaging Sciences (2009).*

Locally low rank (LLR) methods are described in:

- *Trzasko J, Manduca A, Borisch E. Local versus global low-rank promotion in dynamic MRI series reconstruction. In: Proceedings 19th Scientific Meeting, ISMRM. Montreal; 2011. p. 4371.*

Our LLR implementation uses cycle spinning, described in:

- *Coifman RR, Donoho DL. Translation-Invariant De-Noising. In: Antoniadis A, Oppenheim G, editors. Wavelets and Statistics. Vol. 103. Lecture Notes in Statistics. New York, NY: Springer New York; 1995. pp. 125–150.*

This also makes use of the proximal optimized gradient method:

- *Taylor AB, Hendrickx JM, Glineur F. Exact Worst-Case Performance of First-Order Methods for Composite Convex Optimization. SIAM J. Optim. 2017;27:1283–1313.*

Low rank plus sparse (L+S) methods are described in:

- *Lin, C. Y. & Fessler, J. A. Efficient Dynamic Parallel MRI Reconstruction for the Low-Rank Plus Sparse Model. IEEE Trans. Comput. Imaging 5, 17–26 (2019).*

## Acknowledgements

Thanks to Philipp Ehses, Jeff Fessler and the FSL team for releasing their excellent toolboxes that we rely on here.

## Contributors
Thomas W. Okell: https://orcid.org/0000-0001-8258-0659

Mark Chiew: https://orcid.org/0000-0001-6272-8783

## How to cite

Okell, Thomas W, & Chiew, Mark. (2022). Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA) Tools: Initial release (v1.0.0). Zenodo. https://doi.org/10.5281/zenodo.6821643

## Copyright

Copyright © 2022, University of Oxford


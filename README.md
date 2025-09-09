# Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA) Tools
[![DOI](https://zenodo.org/badge/513092852.svg)](https://zenodo.org/badge/latestdoi/513092852)

Matlab code for CAPRIA/CAPRIA+S signal simulation (e.g. for flip angle schedule optimisation), image reconstruction and signal quantification.

## Brief Explanation
Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA) is an MRI method that allows the simultaneous and non-invasive acquisition of time-resolved 3D angiograms and perfusion images. This has recently been extended to allow the reconstruction of T1-weighted structural images (CAPRIA+S) from the same raw k-space dataset. See references below for more details.

## CAPRIA+S Signal Simulation (*models*)
Included here are some models for angiographic, perfusion and structural signals, including effects from a train of excitation pulses with constant or variable flip angles (assuming a spoiled gradient echo was used). Models with and without dispersion are now included. Calculations use standard angiographic, perfusion and static tissue relaxation signal models (referenced below) modified to account for the signal attenuation and changes in the amount of excited transverse magnetisation present with a variable flip angle readout. 

Angiographic or perfusion signals can be calculated using calls to `CAPRIASignalSimple.m` (no dispersion) or `CAPRIASignalDisp.m` (with dispersion). Some examples are given in `ExampleCAPRIASignalSimulations.m`. The static tissue signal can be simulated at the acquired timepoints using `CAPRIAStaticSignal.m`, or a more complete (but time-consuming) simulation of the full evolution of the magnetisation (blood and static tissue) can be run using `CompareCAPRIAandCAPRIAplusSSignals.m`, which relies on a simple Bloch equation simulator (which can be found [here](https://github.com/tomokell/bloch_sim)).

## CAPRIA+S image reconstruction (*recon*)
This code allows a variety of techniques to be used to reconstruct CAPRIA(+S) angiograms, perfusion images and/or structural images from raw k-space data in Siemens twix (meas.dat) format. This includes:
- Extracting the relevant parts of the raw data
- Phase correction
- Construction of an "anatomical" image, using the mean label/control data
- Coil sensitivity estimation (with adaptive combine)
- Coil compression
- A variety of reconstruction options:
    - adjoint non-uniform fast Fourier transform (NUFFT), similar to regridding
    - **locally low rank (LLR) [our preferred approach]**
    - regularised conjugate gradient SENSE (cgSENSE)
    - spatially sparse with L2 temporal smoothness
    - xf-sparse (i.e. sparse in the spatial-temporal frequency domain)
    - spatial total variation (xTV)
    - xTV + L2 temporal smoothness
    - low rank plus sparse (L+S)

In our experience, LLR works most robustly on CAPRIA data, with cgSENSE or spatial sparsity with L2 temporal smoothness also working well, although a comprehensive comparison has not been performed.

The main call is to `CAPRIARecon.m`, which can accept various options to turn on or off sections of the reconstruction. Due to the relatively long reconstruction times, the code is structured to save out intermediate files (e.g. coil sensitivity maps), so subsequent calls will avoid re-running the same stages that have already produced results. Some example angiographic and perfusion-like reconstruction calls used with our original 4D CAPRIA implementation are given in `Example_CAPRIA_Angio_Recon.m` and `Example_CAPRIA_Perf_Recon.m`. Updated scripts used with the CAPRIA+S pulse sequence (which includes an inversion pulse just prior to the readout, slab-selective excitation and modified golden ratio looping) are given in `Example_CAPRIAplusS_[Angio/Perf/Perf_Calib/Struct]_Recon.m`. Note that the `Perf_Calib` reconstruction results in time-resolved tissue (structural) images similar to those from the structural approach, but with identical image resolution and signal scaling to the perfusion reconstruction so they can be used for signal calibration.

## Quantification (*quant*)
In this folder are some scripts for quantifying physiological parameters from CAPRIA(+S) angiograms, perfusion images or structural images. Note that the angiography and perfusion fitting algorithms rely on FABBER (a Bayesian model fitting tool, part of FSL). The scripts make a call to `fabber_asl`, which must be compiled from a modified version that incorporates the CAPRIA flip angle modulation and dispersion and can be found at the `pcasldisp` branch [here](https://github.com/tomokell/fabber_models_asl/tree/pcasldisp).

`FitCAPRIAAngioData.m` runs the angiographic model fitting after first constructing a vessel mask to run the fitting within. 

`CAPRIACalibration.m` gives an example of how to calibrate the CAPRIA(+S) perfusion signal, making some assumptions around other directories that are available, but this could be adapted to different scenarios as needed. Note that this Matlab script calls some shell scripts in the `quant/shellscripts` folder, which must be on the `PATH` within the terminal used to launch Matlab.

`ExampleCAPRIAPerfAnalysis.sh` is a shell script which gives an example call to the `oxasl` tool (which can be found [here](https://github.com/physimals/oxasl)), taking in calibrated data from the `CAPRIACalibration.m` script to perform perfusion quantification (incorporating the required options file for FABBER).

`CAPRIAStaticTissueSignalFit.m` phase corrects the data and fits a modified inversion recovery model appropriate for CAPRIA(+S) to reconstructed structural images, yielding maps of M0, T1, inversion pulse efficiency and (optionally) relative B1+.

## Utilities (*utils*)
This folder contains various helper functions that are used by `models`, `recon` and `quant` and should be included in the Matlab path. 

## External dependencies
Tested using MATLAB 2017b, 2019a and 2024b. Some scripts may rely on the following toolboxes:
- Optimization Toolbox
- Signal Processing Toolbox
- Image Processing Toolbox
- Statistics and Machine Learning Toolbox

In addition, some of the simulation scripts require:
- a simple Bloch equation simulator (which can be found [here](https://github.com/tomokell/bloch_sim))

The recon scripts require:
- `mapVBVD` by Philipp Ehses and colleagues, for Siemens raw data reading: [[GitHub Link]](https://github.com/CIC-methods/FID-A/blob/master/inputOutput/mapVBVD)
- `NUFFT` operators within the Image Reconconstruction Toolbox (IRT) by Jeff Fessler and colleagues: [[Download here]](http://web.eecs.umich.edu/~fessler/irt) [[GitHub link]](https://github.com/JeffFessler/mirt)
- `FSL` Matlab tools and utilities for saving results to Nifti format: [[Download here]](https://fsl.fmrib.ox.ac.uk)

Note that we make heavy use of a subset of Mark Chiew's Matlab image reconstruction code, which is included here for simplicity, but the full set can be found on his [website](https://mchiew.github.io/Tools.html) or on [GitHub](https://github.com/mchiew/recon-tools-matlab).

The quantification scripts rely on:
- `FSL` for various utilities: [[Download here]](https://fsl.fmrib.ox.ac.uk)
- `oxasl` for perfusion quantification: [[Download here]](https://github.com/physimals/oxasl)
- A custom version of `fabber`: see the `pcasldisp` branch [[here]](https://github.com/tomokell/fabber_models_asl/tree/pcasldisp).

For some of the plotting, we have also included the handy utility `hatchfill2` and associated license, but the original version can be found [here](https://uk.mathworks.com/matlabcentral/fileexchange/53593-hatchfill2).

## References
A pre-print describing CAPRIA+S is available here:

- *Okell, T.W., Woods J.G. and Chiew M. Combined Angiographic, Structural and Perfusion Radial Imaging using Arterial Spin Labeling. bioRxiv. 2025. [[DOI link](https://doi.org/10.1101/2025.04.29.651006)]*

4D CAPRIA and optimisation using variable flip angle schedules and an LLR reconstruction is published here:

- *Okell, T. W. & Chiew, M. Optimization of 4D Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling. Magnetic Resonance in Medicine. (2022) doi:10.1002/mrm.29558. [[DOI link](https://doi.org/10.1002/mrm.29558)]*

A 2D implementation of CAPRIA was previously described in:

- *Okell, T. W. Combined angiography and perfusion using radial imaging and arterial spin labeling. Magnetic Resonance in Medicine 81, 182–194 (2019).*

Code and data to reproduce the statistical analyses and figures in the 4D CAPRIA paper can be found here:

- *Okell, Thomas W & Chiew, Mark. 4D Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA): Data and Code to Reproduce Statistics and Figures. (Zenodo, 2022). doi:10.5281/zenodo.7390441.* [[Zenodo link](https://zenodo.org/record/7390441)]

The angiographic signal model uses a modified version of the kinetic model described in:

- *Okell, T. W., Chappell, M. A., Schulz, U. G. & Jezzard, P. A kinetic model for vessel-encoded dynamic angiography with arterial spin labeling. Magnetic Resonance in Medicine 68, 969–979 (2012).*

The base perfusion model is the (P)CASL model described in:

- *Buxton, R. B. et al. A general kinetic model for quantitative perfusion imaging with arterial spin labeling. Magnetic Resonance in Medicine 40, 383–396 (1998).*

Perfusion quantification with `oxasl` relies on the tools described here:

- *Chappell, M.A., Kirk, T.F., Craig, M.S., McConnell, F.A.K., Zhao, M.Y., MacIntosh, B.J., Okell, T.W., Woolrich, M.W. BASIL: A toolbox for perfusion quantification using arterial spin labelling. Imaging Neuroscience. 2023; 1: 1–16. [[DOI link](https://doi.org/10.1162/imag_a_00041)]*

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

Thanks to Philipp Ehses, Jeff Fessler, Michael Chappell and the FSL team for releasing their excellent toolboxes that we rely on here.

## Contributors
Thomas W. Okell: https://orcid.org/0000-0001-8258-0659

Joseph G. Woods: https://orcid.org/0000-0002-0329-824X

Mark Chiew: https://orcid.org/0000-0001-6272-8783

## How to cite

Okell, Thomas W, Woods, Joseph G & Chiew, Mark. (2025). Combined Angiography and Perfusion using Radial Imaging and Arterial Spin Labeling (CAPRIA) Tools (v2.0.0). Zenodo. https://doi.org/10.5281/zenodo.6821642

## Copyright

Copyright © 2025, University of Oxford


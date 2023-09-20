# Designing Phase Masks for Under-Display Cameras

This repository provides code for our paper "Designing Phase Masks for Under-Display Cameras", ICCV 2023, Anqi Yang, Eunhee Kang, Hyong-Euk Lee, and Aswin C. Sankaranarayanan. If you find our data and code useful, please cite our work.

<img src="./images/teaser.png" width="400">

## Setup configuration
The setup specifications for the proposed TOLED w/ phase mask is saved in `TOLED_phase_mask_package/PhasePlate_Params.mat`. Each configuration includes the display pixel pitch, focal length of microlens arrays, parameter T that controls the height of phase mask, and whether to use optimized / uniform sampled / fixed height map for each lenslet. You can specify the index of the configuration and simualte its PSFs, captured image, and evaluate its performance.

## Download pre-computed PSFs
We pre-compute the spatially-varying PSFs for all setup configurations. You can find the matfiles [here](https://drive.google.com/file/d/1fN_86QKuJ1nvzxaJcqOJcNIWFHv60Jt3/view?usp=sharing). Please make an `output` folder and extract `densePSFs` inside it.
```
./output/densePSFs/SETUP_NAMES/PSFs.mat
```
Each `PSFs.mat` includs PSFs for TOLED with phase masks `PSFs_y`  and PSFs for TOLED without phase masks `PSFs_y_lensOnly`. `PSFs_y` is an array of 1024 x 1024 x 3. Each row is PSF for one incidient direction and there all a total of 1024 directions.  Each one-dimensional PSF is 1 x 1024.

## Simulating spatially-varying PSFs

If you wish to compute the spatially-varying PSFs yourself, we provide MATLAB script to computes them. Please run the following script:
```matlab
main_compute_PSFs(ids) 
% ids the indices of setup configurations
```
For each setup, the script simulates and saves PSFs for 400 wavelengths between 300 nm and 700 nm and for 1024 directions of the incident light within the field of view.



## Simulate captured image

We provide a script to simulate captured images under the proposed setup, denoise and deblur the captured images, and evaluate the image quality with PSNR and SSIM.

```matlab
main_quantitative_evalution(ids)
% ids is the indices of setup configuration
```

The denoising and deblurring algorithm requires several hyperparameters. We include the tuned hyperparameters for each setup in the pre-computed folders `output/densePSF/\SETUP\NAME/tuned_parameters.mat`.

## Optimize the height maps of phase masks

We provide the optimization scripts and optimized results in this folder:
```matlab
cd optimize_phase_masks
```

**Pre-computed heights** We pre-compute the optimal lenslet heights and save them at `optimize_phase_mask/lambda0s_pitch*_m_*.mat`.

If you wish to optimize the heights by yourself, there are two steps. First, you need to compute the invertibility matrix $V$ where $V_{j,k}$ represents the system invertibility for height $h_k$ and wavelength $\lambda_k$. Second, you will compute the optimal heights by solving the target function in Equation 11-12. 

**Compute invertibility matrix** The code for computing the invertibility matrix is provided here:
```matlab
compute_invertible_scores.m
```

**Optimize for heights** You can solve for the optimal heights by running:
```matlab
optimize_varying_lambda0s.m
```



### Acknowledgement
We use BM3D code from http://www.cs.tut.fi/~foi/GCF-BM3D/index.html#ref_software for image denoising. This work is/was supported by Global Research Outreach program of Samsung Advanced Institute of Technology and the NSF CAREER award CCF-1652569.
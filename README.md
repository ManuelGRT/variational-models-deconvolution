# Variational Models for Image Deconvolution (MATLAB)

This MATLAB project implements three variational models for image deconvolution — the process of restoring an image that has been blurred and potentially corrupted by noise. It also includes an analytical reconstruction method using Fourier transform inverse.

## Models Implemented

### 1. **Deconvolution Model (Tikhonov Fidelity Term)**
Located in: `variational_models/Deconvolution_model.m`

This model assumes only blurring has affected the image. It uses Tikhonov fidelity term, minimizing a cost function with a quadratic penalty to favor smooth solutions.

### 2. **Deconvolution and Denoising Model (Tikhonov Total Variation Regularization)**
Located in: `variational_models/Deconvolution_Denoising_model.m`

This model addresses both blur and noise. It employs total variation (TV, p=1) regularization to preserve edges while reducing noise, which results in sharper reconstructions than using a linear regularization (p=2).

### 3. **Blind Deconvolution Model**
Located in: `variational_models/Blind_Deconvolution_model.m`

Inspired by the model of Perrone and Favaro (2014), this approach simultaneously estimates the original image and the blur kernel from a single degraded observation. It's suitable for real-world scenarios where the blur is unknown.

## Analytical Reconstruction

Script: `utils/analytic_reconstruction.m`

This script computes the inverse of a known blurring operator (in the frequency domain) to analytically reconstruct the original image, serving as a baseline for comparison with the variational approaches.

## Project Structure
``` 
variational-models-deconvolution/ 
├── commons/ # Gradient, divergence, and noise estimation functions
├── utils/ # Blurring kernels, convolution routines, and analytical reconstruction
├── variational_models/ # Implementations of the three main variational models and execution scripts
├── images/ # Example images and results (if available)
├── README.md # Project documentation 
```

## Dependencies

This project requires:

- MATLAB R2019a or later (tested)
- **Image Processing Toolbox** – for image I/O, filtering, and quality metrics like PSNR
- **Signal Processing Toolbox**


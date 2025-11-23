# MATLAB wrapper for BM3D denoising - from Tampere with love

MATLAB wrapper for BM3D for stationary correlated noise (including white noise) for color,
grayscale and multichannel images and deblurring.

BM3D is an algorithm for attenuation of additive spatially correlated
stationary (aka colored) Gaussian noise. This package provides a wrapper
for the BM3D binaries for use for grayscale, color and other multichannel images
for denoising and deblurring.

This implementation is based on
- Y. Mäkinen, L. Azzari, A. Foi, 2020, "Collaborative Filtering of Correlated Noise: Exact Transform-Domain Variance 
for Improved Shrinkage and Patch Matching", in IEEE Transactions on Image Processing, vol. 29, pp. 8339-8354.
- K. Dabov, A. Foi, V. Katkovnik, K. Egiazarian, 2007, "Image Denoising by Sparse 3-D Transform-Domain Collaborative
Filtering", in IEEE Transactions on Image Processing, vol. 16, pp. 2080-2095.

The package depends on the "BM4D" matlab package (included in the `/bm3d` folder). Please see the "BM4D" package for 
supported platforms (available from https://webpages.tuni.fi/foi/GCF-BM3D/#ref_software).

The package is available for non-commercial use only. For details, see LICENSE.

Authors:
	Ymir Mäkinen   <ymir.makinen@tuni.fi>
	Lucio Azzari
	Alessandro Foi




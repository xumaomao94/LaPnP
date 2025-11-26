# LaPnP

This is the source code for our paper "Radio Map Estimation via Latent Domain Plug-and-Play Denoising". You can start by running the demo code.

Remember to cite our paper [1] if you use the code.

## Update (2025.11.26)

- Add demo code for RT data under 256x256 resolution.
- Now all three denoisers support running under both Windows and Mac with M-series chips (thanks to the updates in BM3D).

## Requirements

LaPnP plugs in exsisting denoisers as an update step. Following denoisers are used.

### BM3D

We have included the BM3D codes (matlab) in this repo, which downloads from https://webpages.tuni.fi/foi/GCF-BM3D/index.html. Please refer to [2] for BM3D denoiser.

### DRUnet

We have included the DRUnet codes in this repo, which downloads from https://github.com/cszn/DPIR. There are a few modifications to the original codes, as we need to call functions from Matlab. Please refer to [3] for more details.

- To run DRUnet, download the model "drunet_gray.pth" according to https://github.com/cszn/DPIR/blob/master/model_zoo/README.md, and put it in "Denoisers/DPIR-master/model_zoo/".

### DSG-NLM

It is implemented by ourselves in matlab. Please refer to [4] for details.

## References

[1] L. Xu, L, Cheng, J, Chen, W, Pu, and F. Xiao, "Radio Map Estimation via Latent Domain Plug-and-Play Denoising." submitted for possible publications

[2] K. Dabov, A. Foi, V. Katkovnik, and K. Egiazarian, “Image denoising by sparse 3D transform-domain collaborative filtering,” IEEE Trans. Image Process., vol. 16, no. 8, pp. 2080-2095, August 2007.

[3] K. Zhang, Y. Li, W. Zuo, L. Zhang, L. Van Gool, and R. Timofte, "Plug-and-play image restoration with deep denoiser prior." IEEE Trans. Pattern Analysis and Machine Intelligence, 2021.

[4] S. Sreehari, S. V. Venkatakrishnan, B. Wohlberg, G. T. Buzzard, L. F. Drummy, J. P. Simmons, and C. A. Bouman, "Plug-and-play priors for bright field electron tomography and sparse interpolation." IEEE Trans. Computational Imaging, 2016.

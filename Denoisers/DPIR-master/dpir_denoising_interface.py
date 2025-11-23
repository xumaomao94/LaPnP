import sys
import os.path
import logging

import numpy as np
from collections import OrderedDict

import torch

dir_path = os.path.dirname(os.path.realpath(__file__))
sys.path.append(os.path.join(dir_path, 'utils'))
sys.path.append(os.path.join(dir_path, 'models'))
from utils1 import utils_logger
from utils1 import utils_model
from utils1 import utils_image as util

def dpir(input,NoiseSigma):
    # ----------------------------------------
    # Preparation
    # ----------------------------------------
    # print('test')
    noise_level_img = NoiseSigma         # set AWGN noise level for noisy image
    noise_level_model = noise_level_img  # set noise level for model
    model_name = 'drunet_gray'           # set denoiser model, 'drunet_gray' | 'drunet_color'
    x8 = False                           # default: False, x8 to boost performance
    border = 0                           # shave boader to calculate PSNR and SSIM

    if 'color' in model_name:
        n_channels = 3                   # 3 for color image
    else:
        n_channels = 1                   # 1 for grayscale image

    # dir_path = 'c:\\Users\\xul2\\Codes\\pnp-denoiser\\DPIR-master'
    dir_path = os.path.join(os.getcwd(), 'Denoisers', 'DPIR-master')
    model_pool = 'model_zoo'             # fixed
    task_current = 'dn'                  # 'dn' for denoising

    model_path = os.path.join(dir_path, model_pool, model_name+'.pth')
    device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
    torch.cuda.empty_cache()

    # ----------------------------------------
    # load model
    # ----------------------------------------

    from models.network_unet import UNetRes as net
    model = net(in_nc=n_channels+1, out_nc=n_channels, nc=[64, 128, 256, 512], nb=4, act_mode='R', downsample_mode="strideconv", upsample_mode="convtranspose")
    model.load_state_dict(torch.load(model_path), strict=True)
    model.eval()
    for k, v in model.named_parameters():
        v.requires_grad = False
    model = model.to(device)

    # ------------------------------------
    # (1) img_L
    # ------------------------------------

    # img_name, ext = os.path.splitext(os.path.basename(img))
    # logger.info('{:->4d}--> {:>10s}'.format(idx+1, img_name+ext))
    input = np.reshape(input, input.shape + (1,))
    img_H = input
    # img_L = util.uint2single(img_H)

    img_L = util.single2tensor4(img_H)
    img_L = torch.cat((img_L, torch.FloatTensor([noise_level_model/255.]).repeat(1, 1, img_L.shape[2], img_L.shape[3])), dim=1)
    img_L = img_L.to(device)

    # ------------------------------------
    # (2) img_E
    # ------------------------------------

    if not x8 and img_L.size(2)//8==0 and img_L.size(3)//8==0:
        img_E = model(img_L)
    elif not x8 and (img_L.size(2)//8!=0 or img_L.size(3)//8!=0):
        img_E = utils_model.test_mode(model, img_L, refield=64, mode=5)
    elif x8:
        img_E = utils_model.test_mode(model, img_L, mode=3)

    img_E = util.tensor2uint(img_E)

    return img_E.copy()

if __name__ == '__main__':
    pass
    # img_H = util.imread_uint("C:\\Users\\xul2\\Codes\\pnp-denoiser\\DPIR-master\\testsets\\set12\\01.png", n_channels=1)
    # img_L = util.uint2single(img_H)
    # img_L = img_L[0:51,0:51]
    # img_L = util.uint2single(img_L)
    # img_o = dpir(img_L,15)
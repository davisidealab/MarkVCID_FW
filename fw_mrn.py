#!/usr/bin/python
# -*- coding: utf-8 -*-



import matplotlib.pyplot as plt
import os, sys
import nibabel as nib
import numpy as np
from dipy.io import read_bvals_bvecs
from dipy.core.gradients import gradient_table

import dipy.reconst.dti as dti

from dipy.reconst.dti import (design_matrix, decompose_tensor,
                           from_lower_triangular, lower_triangular)
from dipy.reconst.vec_val_sum import vec_val_vect
from dipy.core.ndindex import ndindex
import scipy.optimize as opt
from scipy.ndimage.filters import gaussian_filter

# Constant declarations
def wls_fit_tensor_fw(W, data, md_data,S0, Diso=3e-3, mask=None, 
                min_signal=1.0e-6, piterations=2, mdreg=2.0e-3, MDm = 0.0006):

    fw_params = np.zeros(data.shape[:-1] + (9,))


    # Prepare mask
    if mask is None:
        mask = np.ones(data.shape[:-1], dtype=bool)
    else:
        if mask.shape != data.shape[:-1]:
            raise ValueError("Mask is not the same shape as data.")
        mask = np.array(mask, dtype=bool, copy=False)


    index = ndindex(mask.shape)
    for v in index:
        if mask[v]:
            params = wls_iter_fw(W, data, md_data, S0, v, min_signal=min_signal,
                        Diso=3e-3, piterations=piterations, mdreg=mdreg, MDm = MDm)
            fw_params[v] = params
               
    return fw_params
    
    
def wls_iter_fw(W, data, md_data, S0 , v, Diso=3e-3, mdreg=2.0e-3,
             min_signal=1.0e-6, piterations=2, MDm = 0.0006):

    MDm1 = MDm

    sig = data[v]
    MD = md_data[v]
    dmatrix = W

    if (MD < mdreg):
        
        SS = S0[v]
        
        W = dmatrix
    
        # Define weights
        S2 = np.diag(sig**2)
    
        # Defining matrix to solve fwDTI wls solution
        WTS2 = np.dot(W.T, S2)
        inv_WT_S2_W = np.linalg.pinv(np.dot(WTS2, W))
        invWTS2W_WTS2 = np.dot(inv_WT_S2_W, WTS2)
    
        # Process voxel if it has significant signal from tissue
        if np.mean(sig) > min_signal and SS > min_signal:
 
            fwsig = np.exp(np.dot(dmatrix,
                                  np.array([Diso, 0, Diso, 0, 0, Diso, 0])))
    
            df = 1  # initialize precision
            flow = 0  # lower f evaluated
            fhig = 1  # higher f evaluated
            ns = 21  # initial number of samples per iteration
#            p = 0
            for p in range(piterations):
                df = df * 0.1
                fs = np.linspace(flow, fhig, num=ns)  # sampling f
                fs[ns-1] = 0.98
    
                SFW = np.array([fwsig, ]*ns)  # repeat contributions for all values
                FS, SI = np.meshgrid(fs, sig)
                SA = SI - FS*SS*SFW.T

                SA[SA <= 0] = min_signal
                y = np.log(SA / (1-FS))
                all_new_params = np.dot(invWTS2W_WTS2, y)
                # Select params for lower F2
                SIpred = (1-FS)*np.exp(np.dot(W, all_new_params)) + FS*SS*SFW.T
                F2 = np.sum(np.square(SI - SIpred), axis=0)
                evals, evecs =decompose_tensor(from_lower_triangular(all_new_params.T))
 
                MD2 = dti.mean_diffusivity(evals)
                FA2 = dti.fractional_anisotropy(evals)
                Mind0 = np.argmin(F2)
                if ( p == 0):
                    MDa = MD2[0]

                    
                if (MD > MDm1):
                    Mind1 = np.argmin(np.abs(MD2 - MDm1))
                    Mind2 = np.argmin(np.abs(FA2 - 3.0*FA2[0]))              
                    Mind1 = np.min([Mind1,Mind2])
                else:
                    MDm2 = 0.00042*MDa/MDm1
                    Mind1 = np.argmin(np.abs(MD2 - MDm2))
                    Mind2 = np.argmin(np.abs(FA2 - 3.0*FA2[0]))              
                    Mind1 = np.min([Mind1,Mind2])                    
                    
                    
                F2S1 =  F2[0] - F2[Mind1] 

                params1 = all_new_params[:, Mind1]                    
                f = fs[Mind1]  # Updated f
                flow = max([f - df,0])  # refining precision
                fhig = min([f + df, 0.98])
                



            fw_params = np.concatenate((params1,np.array([f]), np.array([F2S1])), axis=0)
  

        else:
            fw_params = np.zeros(9)
    else:
        fw_params = np.zeros(9)
        fw_params[7] = 1.0
        
    return fw_params

    
    


def wls_fit_dti(W, data, mask=None, min_signal=1.0e-6):

    fw_params = np.zeros(data.shape[:-1] + (9,))


    # Prepare mask
    if mask is None:
        mask = np.ones(data.shape[:-1], dtype=bool)
    else:
        if mask.shape != data.shape[:-1]:
            raise ValueError("Mask is not the same shape as data.")
        mask = np.array(mask, dtype=bool, copy=False)


    index = ndindex(mask.shape)
    for v in index:
        if mask[v]:
            params = wls_iter_dti(W, data, v, min_signal=min_signal)
            fw_params[v] = params
                
    return fw_params


def wls_iter_dti(W, data , v, min_signal=1.0e-6):
    


    sig = data[v]
    
    # Define weights
    S2 = np.diag(sig**2)
    SI = sig.copy()
 
    # solve fwDTI wls solution
    WTS2 = np.dot(W.T, S2)
    inv_WT_S2_W = np.linalg.pinv(np.dot(WTS2, W))
    invWTS2W_WTS2 = np.dot(inv_WT_S2_W, WTS2)
    
    SI[SI <= 0] = min_signal
    y = np.log(SI)
    params = np.dot(invWTS2W_WTS2, y)
    SIpred = np.exp(np.dot(W, params))
    
    F2 = np.sum(np.square(SI - SIpred), axis=0)  
    dti_params = np.concatenate((params,np.array([0]), np.array([F2])), axis=0)
    

    return dti_params
    



###########################################################################
# Main driver: process all visits listed in study_dirs    

def main():
    
    start_dir = os.getcwd()
    
#    Diso=3e-3
#    min_signal=1.0e-6
#    piterations=2 
    mdreg=2.0e-3

                    
    ## specify data   
    data_dir = sys.argv[1] 
    cur_out_dir = data_dir 
    fn_data = 'data.nii.gz'
    fn_mask = 'brain_mask.nii.gz'
    fn_bval = 'file.bval'
    fn_bvec = 'file.bvec'

       
    print("Read Data")
    
    if not os.path.exists(cur_out_dir):
        os.makedirs(cur_out_dir)
        
        # Change into current directory
    os.chdir(data_dir)
            
        # Load the dti data
    nii = nib.load(fn_data)
    niim = nib.load(fn_mask)
    img_data = nii.get_data()
    brain_mask = niim.get_data()
    bvals, bvecs = read_bvals_bvecs(fn_bval, fn_bvec)
    
    
    # Extract the data corresponding to b = 0 or b =tmd = MD_s    
    
    
    bvals_eq_0_1000 = (bvals == 0) | (bvals == 1000)
    new_bvals = bvals[bvals_eq_0_1000].copy()
    new_bvecs = bvecs[bvals_eq_0_1000].copy()
    data = img_data[:,:,:,bvals_eq_0_1000].copy()
    
    mask1 = brain_mask.copy()
       
    # Construct the gradient table
    gtab = gradient_table(new_bvals, new_bvecs)
    
    W = design_matrix(gtab)
    
    
    ## smooth the data

    fwhm = 1.25
    gauss_std = fwhm / np.sqrt(8 * np.log(2))  # converting fwhm to Gaussian std
    data_smooth = np.zeros(img_data.shape)
    for v in range(img_data.shape[-1]):
        data_smooth[..., v] = gaussian_filter(img_data[..., v], sigma=gauss_std)
        
#    data_img = nib.Nifti1Image(data_smooth.astype(np.float32), nii.affine)    
#    nib.save(data_img, data_dir + '/fw_mrn/'+ fn_smooth)
    data = data_smooth[:,:,:,bvals_eq_0_1000].copy()
    
    
    # start free water processing
    print('Weighted Least Squares')
                    
    S0 = np.mean(data[..., gtab.b0s_mask], axis=-1)
    
    # weighted least squares fit with no free water
    
    dti_params = wls_fit_dti(W, data, mask=mask1, min_signal=1.0e-6)
    
    evals, evecs = decompose_tensor(from_lower_triangular(dti_params))
    
     
    #evals = dti_params[..., :3]
    FA0 = dti.fractional_anisotropy(evals)   
    MD0 = dti.mean_diffusivity(evals)
    
    FA_img = nib.Nifti1Image(FA0.astype(np.float32), nii.affine)    
    nib.save(FA_img, cur_out_dir + '/wls_dti_FA.nii.gz')
    
#    MD1 = MD0.copy()
#    MD1 = MD1*1000.0
#    MD_img = nib.Nifti1Image(MD1.astype(np.float32), nii.affine)    
#    nib.save(MD_img, cur_out_dir + '/MD_wls_dti.nii.gz')
#    
#    RS_img = nib.Nifti1Image(RS0.astype(np.float32), nii.affine)    
#    nib.save(RS_img, cur_out_dir + '/RS_wls_dti.nii.gz')
#    
#    S0_img = nib.Nifti1Image(S00.astype(np.float32), nii.affine)    
#    nib.save(S0_img, cur_out_dir + '/S0_wls_dti.nii.gz')
    
    
    # start of free water
    pCSF = (MD0 > 0.002)
    mCSF = np.mean(MD0[pCSF])    
#    pCSF = [MD0 > 0.002]
#    mCSF = np.mean(MD0[pCSF])
    mdreg1 = 0.002*mCSF/0.0025
    mdreg = np.min([mdreg,mdreg1])
    MDm = 0.0006
    
    print('Free water')
    
    dti_params1 = wls_fit_tensor_fw(W, data, MD0, S0, Diso=3e-3, mask=mask1, 
                            min_signal=1.0e-6, piterations=2, mdreg=mdreg, MDm = MDm)
    
    
    evals, evecs = decompose_tensor(from_lower_triangular(dti_params1))

    FA1 = dti.fractional_anisotropy(evals)   
    MD1 = dti.mean_diffusivity(evals)
#    RS1 = dti_params1[...,8]
#    S01 = dti_params1[...,6]
#    S01 = np.exp(-S01)
    
    FW1 = dti_params1[..., 7]
    FW1_img = nib.Nifti1Image(FW1.astype(np.float32), nii.affine)    
    nib.save(FW1_img, cur_out_dir + '/wls_dti_FW.nii.gz')
    
              
    FA_img = nib.Nifti1Image(FA1.astype(np.float32), nii.affine)    
    nib.save(FA_img, cur_out_dir + '/fwc_wls_dti_FA.nii.gz')
    
    MD1 = MD1*1000.0
    MD_img = nib.Nifti1Image(MD1.astype(np.float32), nii.affine)    
    nib.save(MD_img, cur_out_dir + '/fwc_wls_dti_MD.nii.gz')
    
    
    
#    RS_img = nib.Nifti1Image(RS1.astype(np.float32), nii.affine)    
#    nib.save(RS_img, cur_out_dir + '/RS1_1shell_fwfit.nii.gz')
    
    
    # Return to start directory
    os.chdir(start_dir)
    
    
main()



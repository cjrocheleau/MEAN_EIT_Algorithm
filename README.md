# MEan Atlas Noser-based (MEAN) Algorithm
### C. Rocheleau, J.L. Mueller 
### Colorado State University
MEan Atlas Noser-based (MEAN) algorithm for Electical Impedance Tomography reconstructions using a one-step Gauss-Newton solver using the mean admittivity distribution of an anatomical atlas of infants as an initial condition. Additional use of the anatomical atlas may optionally be used in a Schur Complement-based post-processing method.  
## MEAN_Reconstruction.m 
  MEAN reconstruction algorithm implementation
  ### INPUTS
    V: (16x15 Matrix) Measured voltages from subject in mV representing each of 16 electrodes and 15 trigonimetric current patterns 
    Schur {OPTIONAL} (Struct): Struct containing fields A and b corresponding to the multiplicative A matrix and additive b vector in Schur complement correction processing
    ReconComps {OPTIONAL} (Struct): Struct containing components used for mean baby processing, to be read in as a variable to cut down on time cost of loading. Contents include
      - sigma_old (496x1 vector): Mean admittivity at each point on the joshua tree mesh
      - U (16x15 matrix): Estimated received voltages for each electrode and current pattern for the mean baby
      - D (16x15x496 tensor): containing gradients of  voltage for each current pattern and electrode at each point on Joshua tree mesh
      -voxel_vol_mL (scalar): Volume of each Joshua tree element in mL
    Reg_Params {Optional} (2x1 vector): Vector containing regularization parameters [alpha beta] where alpha is Tikhonov regularization parameter beta is NOSER-style regularization parameter
    joshmap {Optional} (nJ x nJ Matrix): Matrix used for reconstruction images, contains information about the Joshua Tree element corresponding to each pixel 
    
  ### OUTPUTS
    Recon: (nJ x nJ Matrix) Reconstruction results on the Joshua tree mesh 
    sigma_est: (1 x 496 vector) Vector of estimated admittivity values for each element of Joshua tree

 ## run_MEAN_example.m   
   Example usage of MEAN algorithm. Repository contains data for two test cases

Copyright 2025 Colorado State University Research Foundartion

[![DOI](https://zenodo.org/badge/1025197630.svg)](https://doi.org/10.5281/zenodo.16415830)

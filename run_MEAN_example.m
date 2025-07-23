%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% run_MEAN_example
%--------------------------------------------------------------------------
% C Rocheleau, Colorado State University
% 07/17/2025
%--------------------------------------------------------------------------
% This script serves as an example for running the MEAN algorithm function
% MEAN_Reconstruction.m
%--------------------------------------------------------------------------
% Copyright 2025 Colorado State University Research Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

close all; clear all; clc;

% Initialization ----------------------------------------------------------
% Option to automatically plot results
plot_flag = true; % false; % 

% Load in MEAN initial condition data
load('Amatrix_18_50_22_34.mat');

% Load in Schur Complement Matrices
load('Schur_Matrices_kappa_10e-3.mat');
Schur = Schur.Vent; % Assume this is a Ventilation Image

% Load in Test data
load('MEAN_TestCase1.mat');
% load('MEAN_TestCase2.mat');

if plot_flag
    figure('color','w','Outerposition',[50 50 600 900]);
    subplot(3,2,1)
    pcolor(real(Phantom)); shading flat; 
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Conductivity (S/m)';
    axis equal; axis off; 
    text(-30,32.5,'PHANTOM')
    subplot(3,2,2)
    pcolor(imag(Phantom)); shading flat; 
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Susceptivity (S/m)';
    axis equal; axis off; 
end

% Run MEAN without Schur Complement ---------------------------------------
Reg_Params = [2.5e-4 20];
MEAN = MEAN_Reconstruction(V, ReconComps, Reg_Params);

if plot_flag
    subplot(3,2,3)
    pcolor(real(MEAN)); shading flat; 
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Conductivity (S/m)';
    axis equal; axis off; 
    text(-25,32.5,'MEAN')
    subplot(3,2,4)
    pcolor(imag(MEAN)); shading flat; 
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Susceptivity (S/m)';
    axis equal; axis off; 
end

% Run MEAN with Schur Complement ------------------------------------------
Reg_Params = [2.5e-4 20];
MEAN_Schur = MEAN_Reconstruction(V, ReconComps, Reg_Params, Schur);

if plot_flag
    subplot(3,2,5)
    pcolor(real(MEAN_Schur)); shading flat; 
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Conductivity (S/m)';
    axis equal; axis off; 
    text(-25,32.5,'Schur')
    subplot(3,2,6)
    pcolor(imag(MEAN_Schur)); shading flat;
    temp = colorbar; clim([0 0.5]); temp.Label.String = 'Susceptivity (S/m)';
    axis equal; axis off; 
end
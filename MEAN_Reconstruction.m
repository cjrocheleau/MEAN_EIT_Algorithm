%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MEAN_Reconstruction.m
%--------------------------------------------------------------------------
% C Rocheleau, Colorado State University
% 7/14/2023
%--------------------------------------------------------------------------
% This function will reconstruct and EIT image using a least squares fit
% similar to the NOSER algorithm on the Joshua Tree mesh. The difference in
% this case is that the initial guess and gradients used in the
% computations are numerical approximations based on the Anatomical Atlas
% collected by the S2023 MATH633 class at CSU.
%--------------------------------------------------------------------------
% INPUTS
%       V: (16x15 Matrix) Measured voltages from subject in mV
%           representing each of 16 electrodes and 15 trigonimetric current patterns 
%       Schur {OPTIONAL} (Struct): Struct containing fields A and b
%           corresponding to the multiplicative A matrix and additive b
%           vector in Schur complement correction processing
%       ReconComps {OPTIONAL} (Struct): Struct containing components used
%           for mean baby processing, to be read in as a variable to cut down
%           on time cost of loading. Contents include
%           - sigma_old (496x1 vector): Mean admittivity at each point on 
%               the joshua tree mesh
%           - U (16x15 matrix): Estimated received voltages for each
%               electrode and current pattern for the mean baby
%           - D (16x15x496 tensor): containing gradients of  voltage for 
%               each current pattern and electrode at each point on 
%               Joshua tree mesh
%           -voxel_vol_mL (scalar): Volume of each Joshua tree element in
%               mL
%       Reg_Params {Optional} (2x1 vector): Vector containing
%           regularization parameters [alpha beta] where
%           alpha is Tikhonov regularization parameter
%           beta is NOSER-style regularization parameter
%       joshmap {Optional} (nJ x nJ Matrix): Matrix used for
%           reconstruction images, contains information about the Joshua 
%           Tree element corresponding to each pixel 
%--------------------------------------------------------------------------
% OUTPUTS
%       Recon: (nJ x nJ Matrix) Reconstruction results on the Joshua tree
%           mesh 
%       sigma_est: (1 x 496 vector) Vector of estimated admittivity values 
%           for each element of Joshua tree
%--------------------------------------------------------------------------
% REVISION HISTORY
%   10/23: C Rocheleau, Added regularization parameters and different sized 
%       Joshua Tree maps as optional inputs.
%   08/24: C Rocheleau, Updated to allow for complex-valued Schur
%       complement processing.
%   12/24: C Rocheleau, Renamed from MeanBaby to MEAN.
%--------------------------------------------------------------------------
% Copyright 2025 Colorado State University Research Foundation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Recon, sigma_est] = MEAN_Reconstruction(V, varargin) 
    % Bookkeeping to make V a full matrix
    if issparse(V)
        V = full(V);
    end

    % Optional Input handling
    for iIn = 1:length(varargin)
        % Check if this is Schur struct
        if isfield(varargin{iIn}, 'A')
            Schur = varargin{iIn};
        elseif isfield(varargin{iIn}, 'U')
            ReconComps = varargin{iIn};
        elseif isnumeric(varargin{iIn}) && length(varargin{iIn}) == 2
            alpha = varargin{iIn}(1);
            beta = varargin{iIn}(2);
        elseif isnumeric(varargin{iIn}) && ~isempty(varargin{iIn})
            joshmap = varargin{iIn};
        end
    end

    if ~exist('ReconComps','var')
        load('Amatrix_18_50_22_34.mat');
    end

    % Load in Data to translate Meanbaby to Joshua Tree
    if ~exist('joshmap','var')
        load('joshua_map_real_shape_64_baby_ACT5.mat','joshmap')
    end
     
    % Compute F and F' references in equation 19 of NOSER paper
    % F = squeeze(-2*sum(sum((V-Recon_Comps.U).*Recon_Comps.D)));
    tempD = reshape(ReconComps.D, 16*15, 496); 
    F = -2*tempD.'*(V(:) - ReconComps.U(:));
    
% Using Tikhonov and NOSER-style normalization
    if ~exist('alpha','var') || ~exist('beta','var') 
        alpha = 2.5e-4; 
        beta = 20;
    end
    

    Fprime = 2*(tempD.'*tempD);
    Fprime = Fprime + beta*diag(diag(Fprime));
    alpha  = alpha*eye(size(Fprime));
    Fprime = Fprime + alpha*max(diag(Fprime));
    x = Fprime\F;

    sigma_est = ReconComps.sigma_old  - x;

    Recon = (NaN + NaN*1i)*joshmap; 
    % Recon(joshmap > 0) = real(sigma_est(joshmap(joshmap > 0)));
    Recon(joshmap > 0) = sigma_est(joshmap(joshmap > 0));

    if nargout == 0
        figure()
        pcolor(real(Recon)); shading flat; colorbar;
        clim([0 0.5]);
        title('Mean Initial Condition Recon');
    end

    if exist('Schur','var')
        A = Schur.A; b = Schur.b;
        if isstruct(A)
            Rtemp = reshape(A.Cond*real(sigma_est) + b.Cond, 64, 64);
            Itemp = reshape(A.Susc*imag(sigma_est) + b.Susc, 64, 64);
            Recon = Rtemp + 1i*Itemp;
        elseif size(A,2) == 496
            Recon = reshape(A*sigma_est, 64, 64) + b;
        elseif size(A,2) == 992
            temp  = A*[real(sigma_est); imag(sigma_est)] + b;
            Recon = temp(1:4096) + 1i*temp(4097:end);
            Recon = reshape(Recon, 64, 64);
        end
        Recon(abs(Recon) == 0) = NaN + NaN*1i;
        
        if nargout == 0
            figure()
            pcolor(real(Recon)); shading flat; colorbar; axis ij;
            clim([0 0.5]);
            title('Mean Initial Condition Recon w Schur');
        end
    end
end

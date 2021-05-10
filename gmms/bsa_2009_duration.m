function [ median, sigma, tau, phi ] = bsa_2009_duration(rup,site,dur_type)
%
% Created by Jack Baker, October 19, 2015
% Updated by Emily Mongold, 4/12/2021
%
% Source Model:
% Bommer, J. J., Stafford, P. J., and Alarcón, J. E. (2009). "Empirical 
% Equations for the Prediction of the Significant, Bracketed, and Uniform 
% Duration of Earthquake Ground Motion." Bulletin of the Seismological 
% Society of America, 99(6), 3217-3233.
%
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   rup           = rupture object input containing the following
%                   variables:
%       M             = earthquake magnitude
%       Rrup          = closest distance to rupture (km)
%       Ztor          = depth to top of rupture (km)
%   site          = site object input containing the following
%                   variable:
%       Vs30      = average shear velocity over the first 30m of soil
%   dur_type      = 1 for 5-75% horizontal significant duration
%                 = 2 for 5-75% vertical significant duration
%                 = 3 for 5-95% horizontal significant duration
%                 = 4 for 5-95% vertical significant duration
%                 Locally, Def   = 1 for Ds5-75, =2 for DS5-95
% OUTPUT
%   median          = median significant duration prediction
%   sigma           = log standard deviation of significant duration
%   tau             = within-event log standard deviation
%   phi             = between-event log standard deviation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% check for valid "Def" input
if dur_type == 1 || dur_type == 2
    Def = 1;
elseif dur_type == 3 || dur_type == 4
    Def = 2;
else 
    fprintf('Error--invalid value for input paramter ''Def'' \n');
    return
end

% check for other valid inputs
if isempty(rup.M)||isempty(rup.Rrup)||isempty(site.Vs30)||isempty(rup.Ztor) % any input is empty
    median = nan;
    sigma = nan;
    tau = nan;
    phi = nan;
    return
end

% coefficients (from Table 2)
c0 = [-5.6298 -2.2393];
m1 = [1.2619   0.9368];
r1 = [2.0063   1.5686];
r2 = [-0.252  -0.1953];
h1 = [-2.3316  2.5];
v1 = [-0.29   -0.3478];
z1 = [-0.0522 -0.0365];
tauCoeff  = [0.3527 0.3252];
phiCoeff  = [0.4304 0.3460]; % (named 'sigma' in Table 2, but this is intra-event std and so is renamed to phi here) 
sigma_c   = [0.1729 0.1114];
sigma_Tgm = [0.5289 0.4616];

% log-median (equation 5)
i=Def; % rename variable for convenience below
lnDur = c0(i) + m1(i)*rup.M + (r1(i) + r2(i)*rup.M)*log(sqrt(rup.Rrup.^2+ h1(i).^2)) + v1(i)*log(site.Vs30) + z1(i)*rup.Ztor;
median = exp(lnDur);

% return proper values of standard deviation terms
sigma = sigma_Tgm(i);
tau = tauCoeff(i);
phi = phiCoeff(i);


end


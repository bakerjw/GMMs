function [rho] = gh_2008_spatial_corr(T, h)

% Created by Jack Baker, 10/21/2013
%
% Compute the spatial correlation of epsilons. Documentation is provided 
% in the following document:
%
% Goda, K., and Hong, H. P. (2008). "Spatial Correlation of Peak Ground 
% Motions and Response Spectra." Bulletin of the Seismological Society of 
% America, 98(1), 354-365.
%
% INPUT
%
%   T               = The period of interest. 
%
%   h               = The separation distance between two sites (units of km)
%
% OUTPUT
%
%   rho             = The predicted correlation coefficient



% compute parameters USING CALIFORNIA DATA (equation 11)
alpha = -0.16*log(T) + 0.62; 
beta = 0.5; 
    
% compute parameters USING UNCERTAIN PARAMETERS (equation 12)
% alpha = -0.16*log(T) + 0.68; 
% beta = 0.44; 
    
% compute correlation (equation 10)
rho = exp(-alpha.* (h.^beta) );

end
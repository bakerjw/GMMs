function [rho] = spatial_correlations(T, h, model_name)
% Function to call spectral acceleration correlation functions
% Created by Jack Baker
% 5/10/2021
%
% Note that this function does not evaluate cross-period correlations, even 
% though the lb_2013_spatial_corr and mcb_2018_spatial_corr functions have 
% that capability. To compute cross-period spatial correlations, call the
% functions directly
%   
% INPUT
%
%   T           = The period of interest (scalar). 
%
%   h           = Vector of separation distances between two sites 
%                 (units of km)
%
%   model_name  = the correlation model of interest. One of the following
%                 options:
%                       'gh_2008_spatial_corr'
%                       'hm_2019_spatial_corr'
%                       'jb_2009_spatial_corr'
%                       'lb_2013_spatial_corr'
%                       'mcb_2018_spatial_corr'                 
%
% OUTPUT
%
%   rho         = The predicted correlation coefficient
% 

switch model_name

    case 'gh_2008_spatial_corr'
        rho = gh_2008_spatial_corr(T, h);
    
    case 'hm_2019_spatial_corr'
        rho = hm_2019_spatial_corr(T, h);
    
    case 'jb_2009_spatial_corr'
        rho = jb_2009_spatial_corr(T, h);
    
    case 'lb_2013_spatial_corr'
        rho = lb_2013_spatial_corr(T, T, h);
    
    case 'mcb_2018_spatial_corr'
        rho = mcb_2018_spatial_corr(T, T, h);
            
end


end


function [rho] = sa_correlations(T1, T2, model_name)
% Function to call spectral acceleration correlation functions
% Created by Jack Baker
% 5/7/2021
%   
% INPUT
%
%   T1, T2      = The two periods of interest. The periods may be equal,
%                 and there is no restriction on which one is larger. T1
%                 and T2 may also be vectors
%
%   model_name  = the correlation model of interest. One of the following
%                 options:
%                       'a_2011_corr'
%                       'asa_2014_corr'
%                       'bb_2017_corr'
%                       'bc_2006_corr'
%                       'bj_2008_corr'
%                       'ga_2009_corr'
%
% OUTPUT
%
%   rho         = The predicted correlation coefficient
% 

switch model_name
    
    case 'a_2011_corr'
        rho = a_2011_corr(T1, T2);
    
    case 'asa_2014_corr'
        rho = asa_2014_corr(T1, T2);
    
    case 'bb_2017_corr'
        rho = bb_2017_corr(T1, T2);
    
    case 'bc_2006_corr'
        rho = bc_2006_corr(T1, T2);
    
    case 'bj_2008_corr'
        rho = bj_2008_corr(T1, T2);
        
    case 'ga_2009_corr'
        rho = ga_2009_corr(T1, T2);
        
end


end


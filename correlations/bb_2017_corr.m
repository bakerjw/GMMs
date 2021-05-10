function [rho] = bb_2017_corr(T1, T2)
%
% Created by Jack Baker, 5/7/2021
%
% Compute the correlation of epsilons for spectral accelerations. 
% Calculations use the positive definite version of the tabulated 
% correlations, as documented in:
%
% Baker, J. W., and Bradley, B. A. (2017). “Intensity measure correlations 
% observed in the NGA-West2 database, and dependence of correlations on 
% rupture and site parameters.” Earthquake Spectra, 33(1), 145–156.
%
% The function is fitted over the range the range 0.01s <= T1, T2 <= 10s
%
% INPUT
%
%   T1, T2      = Vectors of the two sets of SA periods of interest. 
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(T1)=n and
%                   length(T2)=m, then size(rho) = [n m]


load bb_2017_data

%% compute correlations
rho = nan*ones(length(T1), length(T2)); % initialize nans

% find indices of input periods that are within range
t1Idx = find(T1>=0.01 & T1<=10);
t2Idx = find(T2>=0.01 & T2<=10);

[X, Y] = meshgrid(log(T2(t2Idx)), log(T1(t1Idx)));

% compute correlations at periods that are within range
rho(t1Idx, t2Idx) = interp2(log(Periods), log(Periods), rhoData, X, Y);
            
end
        



function [rho] = ga_2009_corr(T1, T2)
%
% Created by Jack Baker, May 18, 2016
% Compute the correlation of epsilons using equation 8 of
% 
% Goda, K., and Atkinson, G. M. (2009). "Probabilistic Characterization of 
% Spatially Correlated Response Spectra for Earthquakes in Japan." Bulletin 
% of the Seismological Society of America, 99(5), 3003-3020.
%
% INPUT
%
%   T1, T2      = Vectors of the two sets of SA periods of interest. 
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(T1)=n and
%                   length(T2)=m, then size(rho) = [n m]

% model coefficients from the bottom of page 3009
theta1 = 1.374;
theta2 = 5.586;
theta3 = 0.728;


for i = 1:length(T1)
    for j = 1:length(T2)
        
        T_min = min(T1(i), T2(j));
        T_max = max(T1(i), T2(j));
        
        if T_max > 5 || T_min < 0.1 % limits based on plot limits of figure 6
            rho(i,j) = nan;
        else
            rho(i,j) = 1/3 * (1-cos(pi/2 - (theta1 + theta2 * (T_min < 0.25) ...
                           * (T_min/T_max)^theta3 * log10(T_min/0.25)) * log10(T_max/T_min))) ...
                           +1/3 * (1+cos(-1.5*log10(T_max/T_min)));                       
            if rho(i,j) > 1
                rho(i,j) = 1; % prevent too high of prediction (per 5/18/2016 JWB discussion with KG, this was unintended)
            end
        end
    end
end


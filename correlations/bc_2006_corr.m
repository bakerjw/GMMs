function [rho] = bc_2006_corr(T1, T2, opp)
%
% Created by Jack Baker, 2/2/05, 
% updated 5/7/2021 to allow vector period inputs 
%
% Compute the correlation of epsilons at two different periods. The
% correlation can be for epsilons in the same axis, or perpendicular axes.
%
% The function is strictly empirical, and should not be extrapolated beyond
% the range 0.05s <= T1, T2 <= 5s
%
% Documentation of this model is provided in the following:
%
% Baker J.W. and Cornell C.A. (2006). "Correlation of Response Spectral Values 
% for Multi-Component Ground Motions," Bulletin of the Seismological Society of 
% America, 96 (1), 215-227. 
%
%
% INPUT
%
%   T1, T2      = Vectors of the two sets of SA periods of interest. 
%
%   opp (optional)  = 0 for correlation in the same axis (default)
%                   = 1 for correlation on perpendicular axes
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(T1)=n and
%                   length(T2)=m, then size(rho) = [n m]


% if `opp' is not specified, compute correlations for opp=0
% if 
%     opp=0; 
% end

    
for i = 1:length(T1)
    for j = 1:length(T2)
                
        if (min([T1(i) T2(j)])<0.05 || max([T1(i) T2(j)])>5)
            % periods are out of allowable range
            rho(i,j) = nan;
            
            
        else % periods are in range--compute correlations
            
            T_min = min(T1(i), T2(j));
            T_max = max(T1(i), T2(j));
            
            X = 1 - log(T_max/T_min)*(0.32 - 0.133*exp( -(log(T_max)+1.51)^2 -5*(log(T_min)+2.59)^2));
            D = 9 - log(T_max/T_min^4);
            rho(i,j) = max([X, X-0.17*D, 0.15]);
            
            % adjust for opposite-direction correlations, if desired
            if nargin > 2 && opp == 1
                reduction = 0.78* sqrt( (1-0.05*log(T1)) * (1-0.05*log(T2)));
                rho(i,j) = rho(i,j) * reduction;
            end
        end
    end
end
        
        
        
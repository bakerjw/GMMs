function [median, sigma, period1] = c_1997_vert(T,rup,site,D)

% Created by Jack Baker, 2/1/05, bakerjw@stanford.edu
% Updated by Emily Mongold, 11/27/20
%
% Compute the Campbell attenuation prediction for vertical motions 
% 
% Source Model: 
% Campbell, K. W. (1997). “Empirical Near-Source Attenuation Relationships 
% for Horizontal and Vertical Components of Peak Ground Acceleration, Peak 
% Ground Velocity, and Pseudo-Absolute Acceleration Response Spectra.” 
% Seismological Research Letters, 68(1), 154–179.
%
% Note that this function makes a call to the external function
% "campbell_atten" to compute the horizontal attenuation prediction. Make
% sure that you have this function located in a directory accessible to
% Matlab.
% *EDITED to call external function 'c_1997_active.m' which has a different
% sigma calculation than 'campbell_atten' had, but is the same in every
% other way
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = period of vibration
%   rup           = rupture object input containing the following
%                   variables:
%       M             = earthquake magnitude
%       Rrup          = closest distance to fault rupture
%       lambda        = rake angle, used to set Fault_Type:
%                         = 0 for strike-slip fault
%                         = 1 for reverse, thrust, reverse-oblique, and thrust-oblique fault 
%   site          = site object input containing the following
%                   variable:
%       is_soil       = 0 for soil
%                     = 1 for soft rock
%                     = 2 for hard rock
%   D             = depth to basement bedrock (km)
% OUTPUT   
%   median          = median spectral acceleration prediction
%   sigma           = logarithmic standard deviation of spectral acceleration
%                     prediction
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Setting fault type based on rake angle
if rup.lambda <= 30 && rup.lambda >= -30 || rup.lambda <= 180 && rup.lambda >= 150 || rup.lambda <= -150 && rup.lambda >= -180  
    Fault_Type = 0; %Strike-slip
else
    Fault_Type = 1; %Other/unspecified
end   
% given period
period = [0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00];
% use the full period vector if T = 1000
if length (T) == 1 && T == 1000
    median=zeros(1,length(period));
    sigma=zeros(1,length(period));
    period1=period;
else
    median=zeros(1,length(T));
    sigma=zeros(1,length(T));
    period1=T;
end

% iterate in the case that T is a vector
for n = 1:length(period1)

% interpolate between periods if neccesary    

if (isempty(find(period == period1(n))))
    index_low = sum(period<period1(n));
    T_low = period(index_low);
    T_hi = period(index_low+1);
    
    [sa_low, sigma_low,~] = c_1997_vert(T_low,rup,site,D);
    [sa_hi, sigma_hi,~] = c_1997_vert(T_hi,rup,site,D);
    
    x = [log(T_low) log(T_hi)];
    Y_sa = [log(sa_low) log(sa_hi)];
    Y_sigma = [sigma_low sigma_hi];
    median(n) = exp(interp1(x,Y_sa,log(period1(n))));
    sigma(n) = interp1(x,Y_sigma,log(period1(n)));
    
else % no interpolation needed

    % compute Sa horizontal
    [sa_H, sigma_H,~] = c_1997_active(period1(n),rup,site,D);
    
    % get coefficients
    [c1, c2, c3, c4, c5] = get_coefs(period1(n));
    
    % compute Sa_vertical (Equation 13)
    ln_Sa = log(sa_H) + c1 - 0.1*rup.M + c2*tanh(0.71*(rup.M-4.7)) + c3*tanh(0.66*(rup.M-4.7)) ...
        - 1.5*log(rup.Rrup+0.071*exp(0.661*rup.M)) + 1.89*log(rup.Rrup+0.361*exp(0.576*rup.M)) ...
        - 0.11*Fault_Type + c4*tanh(0.51*D) + c5*tanh(0.57*D);
    median(n) = exp(ln_Sa);
    sigma(n) = sqrt(sigma_H^2 + 0.39^2);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    local function
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [c1, c2, c3, c4, c5] = get_coefs(T)
	% horizontal Sa coefficients from Table 6
	
    period = [0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00];
	c1 = [-1.32, -1.21, -1.29, -1.57, -1.73, -1.98, -2.03, -1.79, -1.82, -1.81, -1.65, -1.31, -1.35];
    c2 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.46, 0.67, 1.13, 1.52, 1.65, 1.28, 1.15];
    c3 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -0.74, -1.23, -1.59, -1.98, -2.23, -2.39, -2.03];
    c4 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.18, 0.57, 0.61, 1.07, 1.26];
    c5 = [0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00, -0.18, -0.49, -0.63, -0.84, -1.17];
    
    index = find(period == T);
    c1 = c1(index);
	c2 = c2(index);
	c3 = c3(index);
	c4 = c4(index);
	c5 = c5(index);
    
    
    
function [median, sigma,period1] = c_1997_active(T,rup,site,D)

% Created by Jack Baker, 2/1/05, bakerjw@stanford.edu
% Updated by Emily Mongold, 11/25/2020
%
% Source Model:
% Campbell, K. W. (1997). “Empirical Near-Source Attenuation Relationships 
% for Horizontal and Vertical Components of Peak Ground Acceleration, Peak 
% Ground Velocity, and Pseudo-Absolute Acceleration Response Spectra.” 
% Seismological Research Letters, 68(1), 154–179.
%
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
%   median        = median spectral acceleration prediction
%   sigma         = logarithmic standard deviation of spectral acceleration
%                   prediction
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

% interpolate between periods if neccesary    
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


if (isempty(find(period == period1(n), 1)))
    index_low = sum(period<period1(n));
    T_low = period(index_low);
    T_hi = period(index_low+1);
    
    [sa_low, sigma_low,~] = c_1997_active(T_low,rup,site,D);
    [sa_hi, sigma_hi,~] = c_1997_active(T_hi,rup,site,D);
    
    x = [log(T_low) log(T_hi)];
    Y_sa = [log(sa_low) log(sa_hi)];
    Y_sigma = [sigma_low sigma_hi];
    median(n) = exp(interp1(x,Y_sa,log(period1(n))));
    sigma(n) = interp1(x,Y_sigma,log(period1(n)));
    
else % no interpolation needed

    % convert soil parameter
    S_SR = 0;
    S_HR = 0;
    if(site.is_soil == 1)
        S_SR = 1;
    elseif(site.is_soil == 2)
        S_HR = 1;
    end
    
    % Compute PGA
    [Ah, sigma_lnAH] = get_PGA(rup.M,rup.Rrup,Fault_Type,S_SR,S_HR);
    
    % get coefficients
    [c1, c2, c3, c4, c5, c6, c7, c8] = get_coefs(period1(n));
    
    % compute Sa (Equation 8)
    [depth_term] = f_SA(D, c6, S_HR, S_SR);
    ln_sa = log(Ah) + c1 + c2*tanh(c3*(rup.M-4.7)) + (c4+c5*rup.M)*rup.Rrup + 0.5*c6*S_SR + c6*S_HR + c7*tanh(c8*D)*(1-S_HR) + depth_term;
    median(n) = exp(ln_sa);
    
    % compute sigma (Equations 5 and 10)
    if rup.M < 7.4
        sig = 0.889-0.0691.*rup.M;
    elseif rup.M >= 7.4
        sig = 0.38;
    end
    sigma(n) = sqrt(sig^2 + 0.27^2);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [Ah, sigma_lnAH] = get_PGA(M, R, F, S_SR, S_HR)
	% from Equation 3
	
    ln_Ah = -3.512 + 0.904*M - 1.328*log( sqrt(R^2 + (.149*exp(0.647*M))^2)) + (1.125-0.112*log(R)-0.0957*M)*F + (.44-.171*log(R))*S_SR + (.405-.222*log(R))*S_HR;
	Ah = exp(ln_Ah);
	
	if (Ah < 0.068)
        sigma_lnAH = 0.55;
	elseif (Ah <= 0.21)
        sigma_lnAH = 0.173 - 0.14 * ln_Ah;
	else
        sigma_lnAH = 0.39;
	end

function [c1, c2, c3, c4, c5, c6, c7, c8] = get_coefs(T)
	% horizontal Sa coefficients from Table 5
	
    period = [0.05, 0.075, 0.10, 0.15, 0.20, 0.30, 0.50, 0.75, 1.00, 1.50, 2.00, 3.00, 4.00];
	c1 = [.05 .27 .48 .72 .79 .77 -.28 -1.08 -1.79 -2.65 -3.28 -4.07 -4.26];
	c2 = [0.00	0.00 0.00 0.00 0.00 0.00 0.74 1.23 1.59 1.98 2.23 2.39 2.03];
	c3 = [0 0 0 0 0 0 0.66 0.66 0.66 0.66 0.66 0.66 0.66];
	c4 = [-0.0011 -0.0024 -0.0024 -0.001 0.0011 0.0035 0.0068 0.0077 0.0085 0.0094 0.01 0.0108 0.0112];
	c5 = [0.000055 0.000095 0.000007 -0.00027 -0.00053 -0.00072 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001 -0.001];
	c6 = [0.2 0.22 0.14 -0.02 -0.18 -0.4 -0.42 -0.44 -0.38 -0.32 -0.36 -0.22 -0.3];
	c7 = [0 0 0 0 0 0 0.25 0.37 0.57 0.72 0.83 0.86 1.05];
	c8 = [0 0 0 0 0 0 0.62 0.62 0.62 0.62 0.62 0.62 0.62];
    
    index = find(period == T);
	c1 = c1(index);
	c2 = c2(index);
	c3 = c3(index);
	c4 = c4(index);
	c5 = c5(index);
	c6 = c6(index);
	c7 = c7(index);
	c8 = c8(index);

function [depth_term] = f_SA(D, c6, S_HR, S_SR)
	% from equations below Equation 8
	
    if (D >= 1)
        depth_term = 0;
	else
        depth_term = c6*(1-S_SR)*(1-D) + 0.5*c6*(1-D)*S_SR;
	end

function [median, sigma, period1] = as_1997_vert(T,rup,site)

% Created by Jack Baker, 2/1/05, bakerjw@stanford.edu
% Updated by Emily Mongold, 11/25/20
%
% Purpose: Compute the Abrahamson and Silva attenuation prediction for 
% spectral accelerations of the vertical component
%
% Source Model:
% Abrahamson, N. A., and Silva, W. J. (1997). "Empirical Response Spectral
% Attenuation Relations for Shallow Crustal Earthquakes."Seismological 
% Research Letters, 68(1), 94-127. 
%
% This script has been modified to correct an error in the function
% f_3 that occured in the printed version of the attenuation.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = period of vibration
%   rup           = rupture object input containing the following
%                   variables:
%       M               = earthquake magnitude
%       R               = closest distance to fault rupture
%       lambda          = rake angle, used to set FaultType:
%                           = 1 for Reverse
%                           = 0.5 for reverse/oblique
%                           = 0 otherwise
%       HW              = 1 for Hanging Wall sites
%                       = 0 otherwise
%   site          = site object input containing the following
%                   variables:
%       is_soil         = 0 (soil),1(soft rock),2(hard rock), locally, S:
%                       = 1 for soil prediction
%                       = 0 for rock
% OUTPUT   
%   median          = median spectral acceleration prediction
%   sigma           = logarithmic standard deviation of spectral acceleration
%                     prediction
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if rup.lambda >= 60 && rup.lambda <= 120  % Converting lambda input to model specific fault type
    FaultType = 1;
elseif rup.lambda >= 30 && rup.lambda <= 60 || rup.lambda >= 120 && rup.lambda <= 150
    FaultType = 0.5;
else
    FaultType = 0;
end


% for the given period T, get the index for the constants
period = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.12, 0.15, 0.17, 0.2, 0.24, 0.3, 0.36, 0.4, 0.46, 0.5, 0.6, 0.75, 0.85, 1, 1.5, 2, 3, 4, 5];

% fill in missing input parameters with default values
if isempty(site.is_soil) % no soil type supplied
    S = 1;
elseif site.is_soil == 0
    S = 1; 
elseif site.is_soil == 1 || site.is_soil == 2
    S = 0;
end

if isempty(rup.HW) % no Hanging wall indicator supplied
    HW = 0;
else
    HW = rup.HW;
end

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
    
    [sa_low, sigma_low,~] = as_1997_vert(T_low,rup,site);
    [sa_hi, sigma_hi,~] = as_1997_vert(T_hi,rup,site);
    
    x = [T_low T_hi];
    Y_sa = [log(sa_low) log(sa_hi)];
    Y_sigma = [sigma_low sigma_hi];
    median = exp(interp1(x,Y_sa,period1(n)));
    sigma = interp1(x,Y_sigma,period1(n));
else
    index = find(period == period1(n));
    
    % get constants for the given index value
    V = get_abrahamson_silva_constants(index);
    R = sqrt(rup.Rrup^2 + V.c4^2);
    

    % compute the PGA term, if we need it
%     S = site.is_soil;
    if ((index ~= 0) || (S ~= 0))
        pga_constants = get_abrahamson_silva_constants(1);
        rock_S = 0;
        pga_rock = exp(calc_val(rup.M,R,pga_constants,FaultType,rock_S,HW,0));
    else
        pga_rock = 0;
    end
    
    median(n) = exp(calc_val(rup.M,R,V,FaultType,S,HW,pga_rock));
    sigma(n) = abrahamson_silva_sigma(rup.M,index);
end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [f1] = f_1(M, R, V)
% value of f1
if (M <= V.c1) 
        f1 = V.a1 + V.a2 * (M - V.c1) + V.a12 * (8.5 - M) ^ V.n + (V.a3 + V.a13 * (M - V.c1)) * log(R);
else 
        f1 = V.a1 + V.a4 * (M - V.c1) + V.a12 * (8.5 - M) ^ V.n + (V.a3 + V.a13 * (M - V.c1)) * log(R);
end

function [f3] = f_3(M, V)
% value of f_3
if M <= 5.8 
        f3 = V.a5;
elseif M < V.c1
%         f3 = V.a5 + (V.a6 - V.a5) / (V.c1 - 5.8);
        f3 = V.a5 + (V.a6 - V.a5) / (V.c1 - 5.8) * (M-5.8);
else
        f3 = V.a6;
end

function [f4] = f_4(M, R, V)
% value of f_4
    if (M <= 5.5) 
        f_HW_M = 0;
    elseif (M <= 6.5) 
        f_HW_M = M - 5.5;
    else
        f_HW_M = 1;
    end
    
    if (R <= 4) 
        f_HW_R = 0;
    elseif (R <= 8) 
        f_HW_R = V.a9 * (R-4)/4;
    elseif (R <= 18) 
        f_HW_R = V.a9;
    elseif (R <= 24) 
        f_HW_R = V.a9 * (1- (R-18)/7);
    else
        f_HW_R = 0;
    end
f4 = f_HW_M*f_HW_R;

function [f5] = f_5(pga_rock, V)
% value of f_5
f5 = V.a10 + V.a11 * log(pga_rock + V.c5);

function [X] = calc_val(M, R, constants, F, S, HW, pga_rock)
% calculate predicted value

%assume no hanging walls
X = f_1(M, R, constants) + F*f_3(M, constants) + HW*f_4(M, R, constants) + S*f_5(pga_rock, constants);

    
function [contants] = get_abrahamson_silva_constants(index)
% get relevant constants

% arrays with values by index
period = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.12, 0.15, 0.17, 0.2, 0.24, 0.3, 0.36, 0.4, 0.46, 0.5, 0.6, 0.75, 0.85, 1, 1.5, 2, 3, 4, 5];
c4 = [6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 6.00, 5.72, 5.35, 4.93, 4.42, 4.01, 3.77, 3.45, 3.26, 2.85, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50, 2.50];
a1 = [1.642, 1.642, 2.100, 2.420, 2.620, 2.710, 2.750, 2.730, 2.700, 2.480, 2.170, 1.960, 1.648, 1.312, 0.878, 0.617, 0.478, 0.271, 0.145, -0.087, -0.344, -0.469, -0.602, -0.966, -1.224, -1.581, -1.857, -2.053];
a2 = 0.909;
a3 = [-1.2520, -1.2520, -1.3168, -1.3700, -1.3700, -1.3700, -1.3700, -1.3700, -1.3700, -1.2986, -1.2113, -1.1623, -1.0987, -1.0274, -0.9400, -0.9004, -0.8776, -0.8472, -0.8291, -0.7896, -0.7488, -0.7451, -0.7404, -0.7285, -0.7200, -0.7200, -0.7200, -0.7200];
a4 = 0.275;
a5 = [0.390, 0.390, 0.432, 0.469, 0.496, 0.518, 0.545, 0.567, 0.580, 0.580, 0.580, 0.580, 0.580, 0.580, 0.580, 0.571, 0.539, 0.497, 0.471, 0.416, 0.348, 0.309, 0.260, 0.260, 0.260, 0.260, 0.260, 0.260];
a6 = [-0.050, -0.050, -0.050, -0.050, -0.050, -0.050, -0.050, -0.050, -0.050, -0.017, 0.024, 0.047, 0.076, 0.109, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.150, 0.058, -0.008, -0.100, -0.100, -0.100];
a9 = [0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.630, 0.604, 0.571, 0.533, 0.488, 0.450, 0.428, 0.400, 0.383, 0.345, 0.299, 0.273, 0.240, 0.240, 0.240, 0.240, 0.240, 0.240];
a10 = [-0.140, -0.140, -0.140, -0.140, -0.140, -0.140, -0.129, -0.119, -0.114, -0.104, -0.093, -0.087, -0.078, -0.069, -0.057, -0.048, -0.043, -0.035, -0.031, -0.022, -0.010, -0.004, 0.004, 0.025, 0.040, 0.040, 0.040, 0.040];
a11 = [-0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220, -0.220];
a12 = [0.0000, 0.0000, 0.0000, 0.0000, -0.0002, -0.0004, -0.0007, -0.0009, -0.0010, -0.0015, -0.0022, -0.0025, -0.0030, -0.0035, -0.0042, -0.0047, -0.0050, -0.0056, -0.0060, -0.0068, -0.0083, -0.0097, -0.0115, -0.0180, -0.0240, -0.0431, -0.0565, -0.0670];
a13 = 0.06;
mag1 = 6.4;
c5 = 0.3;
n = 3;
NPer = 27;

contants.period = period(index);
contants.c4 = c4(index);
contants.a1 = a1(index);
contants.a2 = a2;
contants.a3 = a3(index);
contants.a4 = a4;
contants.a5 = a5(index);
contants.a6 = a6(index);
contants.a9 = a9(index);
contants.a10 = a10(index);
contants.a11 = a11(index);
contants.a12 = a12(index);
contants.a13 = a13;
contants.c1 = mag1;
contants.c5 = c5;
contants.n = n;


function [sigma] = abrahamson_silva_sigma(M, index)
% calculate the sigma

b5 = [0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.76, 0.74, 0.72, 0.70, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.69, 0.72, 0.75, 0.78];
b6 = [0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.085, 0.075, 0.063, 0.056, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050];

if M <= 5 
    sigma = b5(index);
elseif M <= 7
    sigma = b5(index) - b6(index) * (M - 5);
else
    sigma = b5(index) - 2 * b6(index);
end
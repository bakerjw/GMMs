function [median, sigma, period1] = as_1997_active(T,rup,site)

% Created by Jack Baker, 2/1/2005, modified 1/12/2010
% Updated by Emily Mongold, 11/25/2020
%
% Purpose: Compute the Abrahamson and Silva ground motion prediction
%
% Source Model:
% Abrahamson, N. A., and Silva, W. J. (1997). "Empirical Response Spectral
% Attenuation Relations for Shallow Crustal Earthquakes." Seismological
% Research Letters, 68(1), 94-127.
%
% This script has also been modified to correct an error in the function
% f_3 that occured in the printed version of the model.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = period of vibration
%   rup           = rupture object input containing the following
%                   variables:
%       M               = earthquake magnitude
%       Rrup            = closest distance to fault rupture
%       lambda          = rake angle, used to set fault_type:
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
    fault_type = 1;
elseif rup.lambda >= 30 && rup.lambda <= 60 || rup.lambda >= 120 && rup.lambda <= 150
    fault_type = 0.5;
else
    fault_type = 0;
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

    [sa_low, sigma_low,~] = as_1997_active(T_low,rup,site);
    [sa_hi, sigma_hi,~] = as_1997_active(T_hi,rup,site);

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
        pga_rock = exp(calc_val(rup.M,R,pga_constants,fault_type,rock_S,HW,0));
    else
        pga_rock = 0;
    end

    median(n) = exp(calc_val(rup.M,R,V,fault_type,S,HW,pga_rock));
    sigma(n) = abrahamson_silva_sigma(rup.M,index);
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%    local functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
        f3 = V.a5 + (V.a6 - V.a5) / (V.c1 - 5.8) * (M-5.8); % includes correction of error in the paper
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

function [X] = calc_val(M,R,constants, F, S, HW, pga_rock)
% calculate predicted value

	X = f_1(M, R, constants) + F*f_3(M, constants) + HW*f_4(M, R, constants) + S*f_5(pga_rock, constants);

function [contants] = get_abrahamson_silva_constants(index)
% get relevant constants

% arrays with values by index
period = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.12, 0.15, 0.17, 0.2, 0.24, 0.3, 0.36, 0.4, 0.46, 0.5, 0.6, 0.75, 0.85, 1, 1.5, 2, 3, 4, 5];
c4 = [5.6, 5.6, 5.6, 5.6, 5.6, 5.6, 5.58, 5.54, 5.5, 5.39, 5.27, 5.2, 5.1, 4.97, 4.8, 4.62, 4.52, 4.38, 4.3, 4.12, 3.9, 3.81, 3.7, 3.55, 3.5, 3.5, 3.5, 3.5];
a1 = [1.64, 1.64, 1.69, 1.78, 1.87, 1.94, 2.037, 2.1, 2.16, 2.272, 2.407, 2.43, 2.406, 2.293, 2.114, 1.955, 1.86, 1.717, 1.615, 1.428, 1.16, 1.02, 0.828, 0.26, -0.15, -0.69, -1.13, -1.46];
a3 = [-1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.145, -1.135, -1.115, -1.079, -1.035, -1.005, -0.988, -0.965, -0.952, -0.922, -0.885, -0.865, -0.838, -0.772, -0.725, -0.725, -0.725, -0.725];
a5 = [0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.61, 0.592, 0.581, 0.557, 0.528, 0.512, 0.49, 0.438, 0.4, 0.4, 0.4, 0.4];
a6 = [0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.26, 0.232, 0.198, 0.17, 0.154, 0.132, 0.119, 0.091, 0.057, 0.038, 0.013, -0.049, -0.094, -0.156, -0.2, -0.2];
a9 = [0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.37, 0.331, 0.309, 0.281, 0.21, 0.16, 0.089, 0.039, 0];
a10 = [-0.417, -0.417, -0.47, -0.555, -0.62, -0.665, -0.628, -0.609, -0.598, -0.591, -0.577, -0.522, -0.445, -0.35, -0.219, -0.123, -0.065, 0.02, 0.085, 0.194, 0.32, 0.37, 0.423, 0.6, 0.61, 0.63, 0.64, 0.664];
a11 = [-0.23, -0.23, -0.23, -0.251, -0.267, -0.28, -0.28, -0.28, -0.28, -0.28, -0.28, -0.265, -0.245, -0.223, -0.195, -0.173, -0.16, -0.136, -0.121, -0.089, -0.05, -0.028, 0, 0.04, 0.04, 0.04, 0.04, 0.04];
a12 = [0, 0, 0.0143, 0.0245, 0.028, 0.03, 0.03, 0.03, 0.028, 0.018, 0.005, -0.004, -0.0138, -0.0238, -0.036, -0.046, -0.0518, -0.0594, -0.0635, -0.074, -0.0862, -0.0927, -0.102, -0.12, -0.14, -0.1726, -0.1956, -0.215];
mag1 = 6.4;
n = 2;
a2 = 0.512;
a4 = -0.144;
a13 = 0.17;
c5 = 0.03;

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

% use the published coefficients for the geometric mean
period = [0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.075, 0.09, 0.1, 0.12, 0.15, 0.17, 0.2, 0.24, 0.3, 0.36, 0.4, 0.46, 0.5, 0.6, 0.75, 0.85, 1, 1.5, 2, 3, 4, 5];
b5 = [0.7, 0.7, 0.7, 0.705, 0.713, 0.72, 0.728, 0.735, 0.739, 0.746, 0.754, 0.759, 0.765, 0.772, 0.78, 0.787, 0.791, 0.796, 0.799, 0.806, 0.814, 0.819, 0.825, 0.84, 0.851, 0.866, 0.877, 0.885];
b6 = [0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.135, 0.132, 0.13, 0.127, 0.123, 0.121, 0.118, 0.11, 0.105, 0.097, 0.092, 0.087];

if M <= 5
    sigma = b5(index);
elseif M <= 7
    sigma = b5(index) - b6(index) * (M - 5);
else
    sigma = b5(index) - 2 * b6(index);
end

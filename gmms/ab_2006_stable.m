function [median, sigma,period1] = ab_2006_stable(T,rup,site,stress)

% Created by Lynne Schleiffarth, 1/25/10, lynne.schleiffarth@stanford.edu
%
% Modified 3/3/2015 by Jack Baker to make sigma in units of natural log,
% fix units of PGV
%
% Modified 12/5/2016 by Jack Baker to fix an error in the implementation of
% equation 8b of the original paper
%
% Modified 12/4/2017 by Jack Baker to fix an error in a function name (line
% 153)
%
% Updated by Emily Mongold, 11/25/2020
%
% Source Model: 
% Atkinson, G. M., Boore, D. M. (2006). "Earthquake Ground-Motion 
% Prediction Equations for Eastern North America." Bulletin of the 
% Seismological Society of America, 96(6), 2181-2205.
%
% Revised February 19, 2010 to include stress drop factor erratum published in 
% Atkinson, G. M., and Boore, D. M. (2007). “Erratum: Earthquake Ground-
% Motion Prediction Equations for Eastern North America.” Bulletin of the 
% Seismological Society of America, 97(3), 1032–1032.

%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T            = period (s); use period = 0.0 for PGA computation
%                        use period = -1.0 for PGV computation
%   rup          = rupture object input containing the following
%                   variables:
%       M            = moment magnitude
%       Rrup         = closest distance to fault (km) (Rcd in paper)
%   site         = site object input containing the following
%                   variable:
%       Vs30         = shear wave velocity averaged over top 30m of soil (m/s)
%                    <= 760m/s for soil sites (typically; equations are not empirically
%                       constrained for sites with shear wave velocities over 760)
%                    >= 2000m/s for rock sites
%   stress        = fault stress parameter (bars)
%                 = 140 for typical faults
%                 = must be between 35 bars and 560 bars
% OUTPUT
%   median = median spectral acceleration prediction (g)
%   sigma  = NATURAL LOG standard deviation of spectral acceleration
%            prediction (g)
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%----------------------COEFFICIENTS FOR ROCK SITES (TABLE 6)--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
period = [5.0 4.0 3.13 2.5 2.0 1.59 1.25 1.0 0.794 0.629 0.5 0.397 0.315 0.251 0.199 0.158 0.125 0.1 0.079 0.063 0.05 0.04 0.031 0.025 0.0 -1.0];
cr = [-5.41E0 -5.79E0 -6.04E0 -6.17E0 -6.18E0 -6.04E0 -5.72E0 -5.27E0 -4.6E0 -3.92E0 -3.22E0 -2.44E0 -1.72E0 -1.12E0 -6.15E-1 -1.46E-1 2.14E-1 4.8E-1 6.91E-1 9.11E-1 1.11E0 1.26E0 1.44E0 1.52E0 9.07E-1 -1.44E0;...
    1.71E0 1.92E0 2.08E0 2.21E0 2.3E0 2.34E0 2.32E0 2.26E0 2.13E0 1.99E0 1.83E0 1.65E0 1.48E0 1.34E0 1.23E0 1.12E0 1.05E0 1.02E0 9.97E-1 9.8E-1 9.72E-1 9.68E-1 9.59E-1 9.6E-1 9.83E-1 9.91E-1;...
    -9.01E-2 -1.07E-1 -1.22E-1 -1.35E-1 -1.44E-1 -1.5E-1 -1.51E-1 -1.48E-1 -1.41E-1 -1.31E-1 -1.2E-1 -1.08E-1 -9.74E-2 -8.72E-2 -7.89E-2 -7.14E-2 -6.66E-2 -6.4E-2 -6.28E-2 -6.21E-2 -6.2E-2 -6.23E-2 -6.28E-2 -6.35E-2 -6.6E-2 -5.85E-2;...
    -2.54E0 -2.44E0 -2.37E0 -2.3E0 -2.22E0 -2.16E0 -2.1E0 -2.07E0 -2.06E0 -2.05E0 -2.02E0 -2.05E0 -2.08E0 -2.08E0 -2.09E0 -2.12E0 -2.15E0 -2.2E0 -2.26E0 -2.36E0 -2.47E0 -2.58E0 -2.71E0 -2.81E0 -2.7E0 -2.7E0;...
    2.27E-1 2.11E-1 2E-1 1.9E-1 1.77E-1 1.66E-1 1.57E-1 1.5E-1 1.47E-1 1.42E-1 1.34E-1 1.36E-1 1.38E-1 1.35E-1 1.31E-1 1.3E-1 1.3E-1 1.27E-1 1.25E-1 1.26E-1 1.28E-1 1.32E-1 1.4E-1 1.46E-1 1.59E-1 2.16E-1;...
    -1.27E0 -1.16E0 -1.07E0 -9.86E-1 -9.37E-1 -8.7E-1 -8.2E-1 -8.13E-1 -7.97E-1 -7.82E-1 -8.13E-1 -8.43E-1 -8.89E-1 -9.71E-1 -1.12E0 -1.3E0 -1.61E0 -2.01E0 -2.49E0 -2.97E0 -3.39E0 -3.64E0 -3.73E0 -3.65E0 -2.8E0 -2.44E0;...
    1.16E-1 1.02E-1 8.95E-2 7.86E-2 7.07E-2 6.05E-2 5.19E-2 4.67E-2 4.35E-2 4.3E-2 4.44E-2 4.48E-2 4.87E-2 5.63E-2 6.79E-2 8.31E-2 1.05E-1 1.33E-1 1.64E-1 1.91E-1 2.14E-1 2.28E-1 2.34E-1 2.36E-1 2.12E-1 2.66E-1;...
    9.79E-1 1.01E0 1E0 9.68E-1 9.52E-1 9.21E-1 8.56E-1 8.26E-1 7.75E-1 7.88E-1 8.84E-1 7.39E-1 6.1E-1 6.14E-1 6.06E-1 5.62E-1 4.27E-1 3.37E-1 2.14E-1 1.07E-1 -1.39E-1 -3.51E-1 -5.43E-1 -6.54E-1 -3.01E-1 8.48E-2;...
    -1.77E-1 -1.82E-1 -1.8E-1 -1.77E-1 -1.77E-1 -1.73E-1 -1.66E-1 -1.62E-1 -1.56E-1 -1.59E-1 -1.75E-1 -1.56E-1 -1.39E-1 -1.43E-1 -1.46E-1 -1.44E-1 -1.3E-1 -1.27E-1 -1.21E-1 -1.17E-1 -9.84E-2 -8.13E-2 -6.45E-2 -5.5E-2 -6.53E-2 -6.93E-2;...
    -1.76E-4 -2.01E-4 -2.31E-4 -2.82E-4 -3.22E-4 -3.75E-4 -4.33E-4 -4.86E-4 -5.79E-4 -6.95E-4 -7.7E-4 -8.51E-4 -9.54E-4 -1.06E-3 -1.13E-3 -1.18E-3 -1.15E-3 -1.05E-3 -8.47E-4 -5.79E-4 -3.17E-4 -1.23E-4 -3.23E-5 -4.85E-5 -4.48E-4 -3.73E-4];
%----------------------STRESS ADJUSTMENT FACTOR COEFFICIENTS (TABLE 7)----------------------------------------------------------------------------------------------------------------------------------------------------------------------------
delta = [0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.15 0.11];
M1 = [6.0 5.75 5.5 5.25 5.0 4.84 4.67 4.5 4.34 4.17 4.0 3.65 3.3 2.9 2.5 1.85 1.15 0.5 0.34 0.17 0.0 0.0 0.0 0.0 0.5 2.0];
Mh = [8.5 8.37 8.25 8.12 8.0 7.7 7.45 7.2 6.95 6.7 6.5 6.37 6.25 6.12 6.0 5.84 5.67 5.5 5.34 5.17 5.0 5.0 5.0 5.0 5.5 5.5];
%----------------------SOIL RESPONSE COEFFICIENTS (TABLE 8)---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
period_soil = [1/0.2 1/0.25 1/0.32 1/0.5 1/0.63 1/1 1/1.3 1/1.6 1/2 1/2.5 1/3.2 1/4 1/5 1/6.3 1/8 1/10 1/12.6 1/15.9 1/20 1/25 1/32 1/40 0 -1];
blin = [-0.752 -0.745 -0.74 -0.73 -0.726 -0.7 -0.69 -0.67 -0.6 -0.5 -0.445 -0.39 -0.306 -0.28 -0.26 -0.25 -0.232 -0.249 -0.286 -0.314 -0.322 -0.33 -0.361 -0.6];
b1 = [-0.3 -0.31 -0.33 -0.375 -0.395 -0.44 -0.465 -0.48 -0.495 -0.508 -0.513 -0.518 -0.521 -0.528 -0.56 -0.595 -0.637 -0.642 -0.643 -0.609 -0.618 -0.624 -0.641 -0.495];
b2 = [0.0 0.0 0.0 0.0 0.0 0.0 -0.002 -0.031 -0.06 -0.095 -0.13 -0.16 -0.185 -0.185 -0.14 -0.132 -0.117 -0.105 -0.105 -0.105 -0.108 -0.115 -0.144 -0.060];
%----------------------COEFFICIENTS FOR SOIL SITES (TABLE 9)--------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
cs = [-4.85E0 -5.26E0 -5.59E0 -5.8E0 -5.85E0 -5.75E0 -5.49E0 -5.06E0 -4.45E0 -3.75E0 -3.01E0 -2.28E0 -1.56E0 -8.76E-1 -3.06E-1 1.19E-1 5.36E-1 7.82E-1 9.67E-1 1.11E0 1.21E0 1.26E0 1.19E0 1.05E0 5.23E-1 -1.66E0;...
    1.58E0 1.79E0 1.97E0 2.13E0 2.23E0 2.29E0 2.29E0 2.23E0 2.12E0 1.97E0 1.8E0 1.63E0 1.46E0 1.29E0 1.16E0 1.06E0 9.65E-1 9.24E-1 9.03E-1 8.88E-1 8.83E-1 8.79E-1 8.88E-1 9.03E-1 9.69E-1 1.05E0;...
    -8.07E-2 -9.79E-2 -1.14E-1 -1.28E-1 -1.39E-1 -1.45E-1 -1.48E-1 -1.45E-1 -1.39E-1 -1.29E-1 -1.18E-1 -1.05E-1 -9.31E-2 -8.19E-2 -7.21E-2 -6.47E-2 -5.84E-2 -5.56E-2 -5.48E-2 -5.39E-2 -5.44E-2 -5.52E-2 -5.64E-2 -5.77E-2 -6.2E-2 -6.04E-2;...
    -2.53E0 -2.44E0 -2.33E0 -2.26E0 -2.2E0 -2.13E0 -2.08E0 -2.03E0 -2.01E0 -2.0E0 -1.98E0 -1.97E0 -1.98E0 -2.01E0 -2.04E0 -2.05E0 -2.11E0 -2.17E0 -2.25E0 -2.33E0 -2.44E0 -2.54E0 -2.58E0 -2.57E0 -2.44E0 -2.5E0;...
    2.22E-1 2.07E-1 1.91E-1 1.79E-1 1.69E-1 1.58E-1 1.5E-1 1.41E-1 1.36E-1 1.31E-1 1.27E-1 1.23E-1 1.21E-1 1.23E-1 1.22E-1 1.19E-1 1.21E-1 1.19E-1 1.22E-1 1.23E-1 1.3E-1 1.39E-1 1.45E-1 1.48E-1 1.47E-1 1.84E-1;...
    -1.43E0 -1.31E0 -1.2E0 -1.12E0 -1.04E0 -9.57E-1 -9.0E-1 -8.74E-1 -8.58E-1 -8.42E-1 -8.47E-1 -8.88E-1 -9.47E-1 -1.03E0 -1.15E0 -1.36E0 -1.67E0 -2.1E0 -2.53E0 -2.88E0 -3.04E0 -2.99E0 -2.84E0 -2.65E0 -2.34E0 -2.3E0;...
    1.36E-1 1.21E-1 1.1E-1 9.54E-2 8.0E-2 6.76E-2 5.79E-2 5.41E-2 4.98E-2 4.82E-2 4.7E-2 5.03E-2 5.58E-2 6.34E-2 7.38E-2 9.16E-2 1.16E-1 1.48E-1 1.78E-1 2.01E-1 2.13E-1 2.16E-1 2.12E-1 2.07E-1 1.91E-1 2.5E-1;... 
    6.34E-1 7.34E-1 8.45E-1 8.91E-1 8.67E-1 8.67E-1 8.21E-1 7.92E-1 7.08E-1 6.77E-1 6.67E-1 6.84E-1 6.5E-1 5.81E-1 5.08E-1 5.16E-1 3.43E-1 2.85E-1 1.0E-1 -3.19E-2 -2.1E-1 -3.91E-1 -4.37E-1 -4.08E-1 -8.7E-2 1.27E-1;...
    -1.41E-1 -1.56E-1 -1.72E-1 -1.8E-1 -1.79E-1 -1.79E-1 -1.72E-1 -1.7E-1 -1.59E-1 -1.56E-1 -1.55E-1 -1.58E-1 -1.56E-1 -1.49E-1 -1.43E-1 -1.5E-1 -1.32E-1 -1.32E-1 -1.15E-1 -1.07E-1 -9E-2 -6.75E-2 -5.87E-2 -5.77E-2 -8.29E-2 -8.7E-2;...
    -1.61E-4 -1.96E-4 -2.45E-4 -2.6E-4 -2.86E-4 -3.43E-4 -4.07E-4 -4.89E-4 -5.75E-4 -6.76E-4 -7.68E-4 -8.59E-4 -9.55E-4 -1.05E-3 -1.14E-3 -1.18E-3 -1.13E-3 -9.9E-4 -7.72E-4 -5.48E-4 -4.15E-4 -3.88E-4 -4.33E-4 -5.12E-4 -6.3E-4 -4.27E-4];
%-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------%
if length (T) == 1 && T == 1000
    median=zeros(1,length(period)-2);
    sigma=zeros(1,length(period)-2);
    period1=period(1:end-2);
else
    median=zeros(1,length(T));
    sigma=zeros(1,length(T));
    period1=T;
end

for n = 1:length(period1)

%**************************************************************************
% FIND C COEFFICIENTS BY INTERPOLATING BETWEEN PERIODS
%**************************************************************************
% find period directly above and below given period in the table
ilow = find(period<=period1(n), 1 );
T_low = period(ilow);
ihigh = find(period>=period1(n), 1, 'last' );
T_high = period(ihigh);
% if given period equals a period in the table, then no need to interpolate
if ihigh==ilow
    if site.Vs30>=2000
        c = cr(:,ihigh);
    else
        c = cs(:,ihigh);
    end
    delta_val = delta(ihigh);
    M1_val = M1(ihigh);
    Mh_val = Mh(ihigh);
% otherwise, interpolate between coeffients
else
    if site.Vs30>=2000
        c_high = cr(:,ihigh);
        c_low = cr(:,ilow);
    else
        c_high = cs(:,ihigh);
        c_low = cs(:,ilow);
    end
    for i=1:length(c_high)
        c(i) = interp1([T_low T_high], [c_low(i) c_high(i)], period1(n));
    end
    delta_val = interp1([T_low T_high], [delta(ilow) delta(ihigh)], period1(n));
    M1_val = interp1([T_low T_high], [M1(ilow) M1(ihigh)], period1(n));
    Mh_val = interp1([T_low T_high], [Mh(ilow) Mh(ihigh)], period1(n));
end

%**************************************************************************
% COMPUTE SOIL RESPONSE FACTOR USING EQN. (7)
%**************************************************************************
% define constants from Table (8)
Vref = 760;
v1 = 180;
v2 = 300;
% find period directly above and below given period in the soil response table
ilow = find(period_soil<=period1(n), 1 );
T_low = period_soil(ilow);
ihigh = find(period_soil>=period1(n), 1, 'last' );
T_high = period_soil(ihigh);
% if shear wave velocity indicates bedrock or is equal to the reference velocity, then no modfication is needed
if site.Vs30>=2000 || site.Vs30==Vref
    S = 0;
else   
    % if period given equals a period in the table, then no need to interpolate
    if ihigh==ilow
        blin = blin(ihigh);
        b1 = b1(ihigh);
        b2 = b2(ihigh);
    % otherwise, interpolate between coefficients
    else
        blin = interp1([T_low T_high], [blin(ilow) blin(ihigh)], period1(n));
        b1 = interp1([T_low T_high], [b1(ilow) b1(ihigh)], period1(n));
        b2 = interp1([T_low T_high], [b2(ilow) b2(ihigh)], period1(n));
    end    
    % compute bnl from eqn. (8)
    if site.Vs30<=v1
        bnl = b1;
    elseif v1<site.Vs30 && site.Vs30<=v2
        bnl = (b1-b2)*log(site.Vs30/v2)/log(v1/v2)+b2;
    elseif v2<site.Vs30 && site.Vs30<=Vref
        bnl = b2*log(site.Vs30/Vref)/log(v2/Vref);
    else
        bnl = 0;
    end
    % compute S from eqn. (7)
    siteref = site; 
    siteref.Vs30 = Vref;
    pgaBC = ab_2006_stable(0,rup,siteref,stress);
    if pgaBC<=60
        S = log10(exp(blin*log(site.Vs30/Vref) + bnl*log(60/100)));
    else
        S = log10(exp(blin*log(site.Vs30/Vref) + bnl*log(pgaBC/100)));
    end    
end

%**************************************************************************
% COMPUTE LOG MEAN SPECTRAL ACCELERAION USING EQN. (5)
%**************************************************************************
% define constants from eqn. (5)
R0 = 10;
R1 = 70;
R2 = 140;
f0 = max([log10(R0/rup.Rrup) 0]);
f1 = min([log10(rup.Rrup) log10(R1)]);
f2 = max([log10(rup.Rrup/R2) 0]);
% eqn. (5)
logSa = c(1) + c(2).*rup.M + c(3).*rup.M.^2 + (c(4)+c(5).*rup.M).*f1 + ...
    (c(6)+c(7).*rup.M).*f2 + (c(8)+c(9).*rup.M).*f0 + c(10).*rup.Rrup + S;

%**************************************************************************
% COMPUTE STRESS ADJUSTMENT FACTOR USING EQN. (6)
%**************************************************************************
logSF2 = min((delta_val+0.05), (0.05 + delta_val*max(rup.M-M1_val, 0)/(Mh_val-M1_val)));
% compute mean spectral acceleration including adjustment for stress 
% according to erratum published in Bulletin of Seismological Society of
% America, Vol. 96, No. 6, pp. 2181-2205, December 2006
if 35<=stress && stress<=560
    factor = log10(stress/140)/log10(2);
    logSa_adj = logSa + factor*logSF2;
else
    error('Stress parameter is out of suggested range');
end

%**************************************************************************
% CONVERT SPECTRAL ACCELERATION FROM LOG SCALE
%**************************************************************************
if T>=0
    median(n) = (10.^logSa_adj)/981; % convert to units of g
else
    median(n) = (10.^logSa_adj);
end

%**************************************************************************
% DEFINE SIGMA ACCORDING TO TABLE (6) & TABLE (9)
%**************************************************************************
sigma(n) = 0.3 * log(10);
end
end
function [median, sigma, period1] = ab_2003_subduction(T,rup,site,sub_ind,Zl)

% Created by James Bronder 06/09/2010, jbronder@stanford.edu
% Updated by Emily Mongold, 11/25/2020
%
% Purpose: Computes the mean and standard deviation of the PGA
% or psuedoacceleration, PSA, 5% damping. Additional modifications are also
% included for the regions of Cascadia and Japan.

% Source Model: 
% Atkinson, G. M., Boore, D. M. (2003). "Empirical Ground-Motion Relations 
% for Subduction-Zone Earthquakes and Their Application to Cascadia and 
% Other Regions." Bulletin of the Seismological Society
% of America, 93(4), 1703-1729.

% General Limitations: This equation is obtained from a global database of
% ground motions (~1200 horizontal records). It is highly recommended that
% a comparison be drawn using records that compare well with the region of
% interest.
%   The authors suggest applying the equation for these given bounds:

%   For Interface Events:
%   5.5 <= M < 6.5 and Df <= 80 km
%   6.5 <= M < 7.5 and Df <= 150 km
%   M >= 7.5 and Df <= 300 km

%   For Intraslab (In-slab) Events:
%   6.0 <= M < 6.5 and Df <= 100 km
%   M >= 6.5 and Df <= 200 km

%   The established bounds optimize the equations for seismic-hazard
%   analysis.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T        = Period (sec)
%   rup      = rupture object input containing the following variables:
%       M       = Moment Magnitude
%       h       = Focal (Hypocentral) Depth (km) (rup.Zhyp)
%       Rrup    = closest distance to the fault surface (km) (Df locally)
%   site     = site object input containing the following variable:
%       Vs30    = Shear Wave Velocity averaged over the top 30 meters of 
%                 soil of the soil profile (m/sec)
%   sub_ind  = subduction type indicator: 0 for interface, 1 for intraslab
%              'Zt' used locally within function  
%   Zl       = Cascadia or Japan indicator: Zl = 0 for General Cases
%                                           Zl = 1 for Cascadia
%                                           Zl = 2 for Japan
% OUTPUT
%   median   = Median spectral acceleration prediction (g)
%   sigma    = Logarithmic standard deviation of spectral acceleration
%            prediction
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zt = sub_ind;

%Period
period = [0 0.04 0.1 0.2 0.4 1 2 1/0.33];

if length (T) == 1 && T == 1000
    median=zeros(1,length(period));
    sigma=zeros(1,length(period));
    period1=period;
else
    median = zeros(1,length(T));
    sigma = zeros(1,length(T));
    period1 = T;
end

%Interslab Events Coefficients

c1_it = [2.991 2.8753 2.7789 2.6638 2.5249 2.1442 2.1907 2.301];
c1_it_jp = [3.14 3.05 2.95 2.84 2.58 2.18 2.14 2.27];
c1_it_cas = [2.79 2.60 2.50 2.54 2.50 2.18 2.33 2.36];
c2_it = [0.03525 0.07052 0.09841 0.12386 0.1477 0.1345 0.07148 0.02237];
c3_it = [0.00759 0.01004 0.00974 0.00884 0.00728 0.00521 0.00224 0.00012];
c4_it = [-0.00206 -0.00278 -0.00287 -0.0028 -0.00235 -0.0011 0 0];
c5_it = [0.19 0.15 0.15 0.15 0.13 0.1 0.1 0.1];
c6_it = [0.24 0.2 0.23 0.27 0.37 0.3 0.25 0.25];
c7_it = [0.29 0.2 0.2 0.25 0.38 0.55 0.4 0.36];
sigma_it = [0.23 0.26 0.27 0.28 0.29 0.34 0.34 0.36];
sigma1_it = [0.2 0.22 0.25 0.25 0.25 0.28 0.29 0.31];
sigma2_it = [0.11 0.14 0.1 0.13 0.15 0.19 0.18 0.18];

%In-slab Events Coefficients

c1_in = [-0.04713 0.50697 0.43928 0.51589 0.005445 -1.02133 -2.39234 -3.70012];
c1_in_jp = [0.10 0.68 0.61 0.70 0.07 -0.98 -2.44 -3.73];
c1_in_cas = [-0.25 0.23 0.16 0.40 -0.01 -0.98 -2.25 -3.64];
c2_in = [0.6909 0.63273 0.66675 0.69186 0.7727 0.8789 0.9964 1.1169];
c3_in = [0.0113 0.01275 0.0108 0.00572 0.00173 0.0013 0.00364 0.00615];
c4_in = [-0.00202 -0.00234 -0.00219 -0.00192 -0.00178 -0.00173 -0.00118 -0.00045];
c5_in = [0.19 0.15 0.15 0.15 0.13 0.1 0.1 0.1];
c6_in = [0.24 0.2 0.23 0.27 0.37 0.3 0.25 0.25];
c7_in = [0.29 0.2 0.2 0.25 0.38 0.55 0.4 0.36];
sigma_in = [0.27 0.25 0.28 0.28 0.28 0.29 0.3 0.3];
sigma1_in = [0.23 0.24 0.27 0.26 0.26 0.27 0.28 0.29];
sigma2_in = [0.14 0.07 0.07 0.1 0.1 0.11 0.11 0.08];

% Preliminary Inital Conditions and Variables Computation
if rup.Zhyp >= 100
    rup.Zhyp = 100;
else
    rup.Zhyp;
end

if Zt == 0 && rup.M >= 8.5
    rup.M = 8.5;
elseif Zt == 0 && rup.M < 8.5
    rup.M;
elseif Zt == 1 && rup.M >= 8.0
    rup.M = 8.0;
elseif Zt == 1 && rup.M < 8.0
    rup.M;
end

delta = 0.00724*(10^(0.507*rup.M));
R = sqrt(rup.Rrup^2 + delta^2);

if Zt == 0
    g = 10^(1.2 - 0.18*rup.M);
elseif Zt == 1
    g = 10^(0.301 - 0.01*rup.M);
end

if site.Vs30 > 760
    Sc = 0;
    Sd = 0;
    Se = 0;
elseif site.Vs30 > 360
    Sc = 1;
    Sd = 0;
    Se = 0;
elseif site.Vs30 >= 180
    Sc = 0;
    Sd = 1;
    Se = 0;
elseif site.Vs30 < 180
    Sc = 0;
    Sd = 0;
    Se = 1;
end

% Begin Computation of Ground Motion Parameter with the modifications
for n = 1:length(period1)
if Zt == 0
    
    if Zl == 0
        c1 = c1_it(1);
    elseif Zl == 1
        c1 = c1_it_cas(1);
    elseif Zl == 2
        c1 = c1_it_jp(1);
    end
    
    log_PGArx = c1 + c2_it(1)*rup.M + c3_it(1)*rup.Zhyp + c4_it(1)*R - ...
        g*log10(R);
    PGArx = 10^(log_PGArx);
    
    if PGArx <= 100 || (1/period1(n)) <= 1
        sl = 1;
    end
    
    if 1/period1(n) >= 2
        if PGArx < 500
            sl = 1 - (PGArx - 100)/400;
        elseif PGArx >= 500
            sl = 0;
        end
    elseif 1/period1(n) < 2
        if PGArx < 500
            sl = 1 - ((1/period1(n))-1)*(PGArx - 100)/400;
        elseif PGArx >= 500
            sl = 1 - ((1/period1(n))-1);
        end
    end
    
elseif Zt == 1
    
    if Zl == 0
        c1 = c1_in(1);
    elseif Zl == 1
        c1 = c1_in_cas(1);
    elseif Zl == 2
        c1 = c1_in_jp(1);
    end
    
    log_PGArx = c1 + c2_in(1)*rup.M + c3_in(1)*rup.Zhyp + c4_in(1)*R - ...
        g*log10(R);
    PGArx =10^(log_PGArx);
    
    if PGArx <= 100 || (1/period1(n)) <= 1
        sl = 1;
    end
    
    if 1/period1(n) >= 2
        if PGArx < 500
            sl = 1 - (PGArx - 100)/400;
        elseif PGArx >= 500
            sl = 0;
        end
    elseif 1/period1(n) < 2
        if PGArx < 500
            sl = 1 - ((1/period1(n))-1)*(PGArx - 100)/400;
        elseif PGArx >= 500
            sl = 1 - ((1/period1(n))-1);
        end
    end
end

if Zt == 0
    
    if isempty(find(period == period1(n)))
        
        i_lo = sum(period<period1(n));
        T_lo = period(i_lo);
        T_hi = period(i_lo + 1);
        
        [Sa_hi sigma_hi] = ab_2003_subduction(T_hi,rup,site,Zt,Zl);
        [Sa_lo sigma_lo] = ab_2003_subduction(T_lo,rup,site,Zt,Zl);
        
        x = [T_lo T_hi];
        Y_Sa = [Sa_lo Sa_hi];
        Y_sigma = [sigma_lo sigma_hi];
        median(n) = interp1(x,Y_Sa,period1(n));
        sigma(n) = interp1(x,Y_sigma,period1(n));
        
    else
        i=find(period == period1(n));
        
        if Zl == 0
            c1 = c1_it(i);
        elseif Zl == 1
            c1 = c1_it_cas(i);
        elseif Zl == 2
            c1 = c1_it_jp(i);
        end
        
        log_10_Y = c1 + c2_it(i)*rup.M + c3_it(i)*rup.Zhyp + c4_it(i)*R - ...
            g*log10(R) + c5_it(i)*sl*Sc + c6_it(i)*sl*Sd + ...
            c7_it(i)*sl*Se; % Log10 Sa in cm/s^2
        sigma_10 = sqrt((sigma1_it(i))^2 + (sigma2_it(i))^2);
        
        median(n) = 10.^(log_10_Y)/981; % Median Sa in g
        sigma(n) = log(10.^sigma_10);
    end
elseif Zt == 1
    if isempty(find(period == period1(n)))
        
        i_lo = sum(period<period1(n));
        T_lo = period(i_lo);
        T_hi = period(i_lo + 1);
        
        [Sa_hi sigma_hi] = ab_2003_subduction(T_hi,rup,site,Zt,Zl);
        [Sa_lo sigma_lo] = ab_2003_subduction(T_lo,rup,site,Zt,Zl);
        
        x = [T_lo T_hi];
        Y_Sa = [Sa_lo Sa_hi];
        Y_sigma = [sigma_lo sigma_hi];
        median(n) = interp1(x,Y_Sa,period1(n));
        sigma(n) = interp1(x,Y_sigma,period1(n));
        
    else
        i=find(period == period1(n));
        
        if Zl == 0
            c1 = c1_in(i);
        elseif Zl == 1
            c1 = c1_in_cas(i);
        elseif Zl == 2
            c1 = c1_in_jp(i);
        end
        
        log_10_Y = c1 + c2_in(i)*rup.M + c3_in(i)*rup.Zhyp + c4_in(i)*R - ...
            g*log10(R) + c5_in(i)*sl*Sc + c6_in(i)*sl*Sd + ...
            c7_in(i)*sl*Se; % Log10 Sa in cm/s^2
        sigma_10 = sqrt((sigma1_in(i))^2 + (sigma2_in(i))^2);
        
        median(n) = 10.^(log_10_Y)/981; % Median Sa in g
        sigma(n) = log(10.^sigma_10);
    end
end
end
end

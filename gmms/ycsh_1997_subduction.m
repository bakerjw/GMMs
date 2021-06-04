function [median, sigma, period1] = ycsh_1997_subduction(T,rup,site,sub_ind)

% Created by James Bronder 06/09/2010, jbronder@stanford.edu
% Updated by Emily Mongold, 11/27/20
%
% Purpose: Computes the median and logaritmic standard deviation of a
%          subduction zone earthquake with 5% damping.
%
% Source Model:
% Youngs, R.R., S.J. Chiou, W.J. Silva, J.R. Humphrey. (1997). "Strong 
% Ground Motion Attentuation Relationships for Subduction Zone 
% Earthquakes." Seismological Research Letters, 68(1), 58-73.

% General Limitations: According to the authors, "The attenuation
% relationships developed in this study are considered appropriate for
% earthquakes of magnitude M 5 and greater and for distances to the rupture
% surface of 10 to 500km."
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T           = Period of Interest (sec.). For PGA Computation, T = 0
%                 For Rock Sites, 0<=T<=3
%                 For Soil Sites, 0<=T<=4
%   rup         = rupture object input containing the following
%                 variables:
%       M           = Moment Magnitude
%       R           = Source to site distance to rupture surface (km)
%       Zhyp        = Hypocentral (Focal) depth from surface to focus (km)
%   site        = site object input containing the following
%                 variable:
%       is_soil     = 0 (soil),1(soft rock),2(hard rock)
%                   = soil type indicator, 'Zr' locally:
%                        For Rock Sites, Zr = 1
%                        Otherwise, Zr = 0
%   sub_ind     = Subduction Zone Type, 'Zt' locally: 
%                   For Intraslab, Zt = 1
%                   For Interface, Zt = 0
% OUTPUT
%   median      = Median spectral acceleration prediction (g)
%   sigma       = Logarithmic standard deviation of spectral acceleration
%                 prediction
%   period1     = periods for which the median and sigma values are
%                 provided. If T = 1000, then period1 = the full set of
%                 available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zt = sub_ind;
r = rup.R;

if site.is_soil == 0 % Mapping generic soil inputs to model-specific values
    Zr = 0; 
else
    Zr = 1;
end

%--------------------------Generic Rock Period----------------------------%
period_rock = [0 0.075 0.10 0.20 0.30 0.40 0.50 0.75 1.00 1.50 2.00 3.00];

%-------------------------Generic Rock Coefficients-----------------------%
grc1  = [0.2418 1.5168 1.4298 0.9638 0.4878 0.1268 -0.1582 -0.9072 -1.4942 -2.3922 -3.0862 -4.2692];
grc2  = [1.414 1.414 1.414 1.414 1.414 1.414 1.414 1.414 1.414 1.414 1.414 1.414];
grc3  = [0.0000 0.0000 -0.0011 -0.0027 -0.0036 -0.0043 -0.0048 -0.0057 -0.0064 -0.0073 -0.0080 -0.0089];
grc4  = [-2.552 -2.707 -2.655 -2.528 -2.454 -2.401 -2.360 -2.286 -2.234 -2.160 -2.107 -2.033];
grc5  = [0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607 0.00607];
grc6  = [0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846 0.3846];
grc7  = [1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818 1.7818];
grc8  = [0.554 0.554 0.554 0.554 0.554 0.554 0.554 0.554 0.554 0.554 0.554 0.554];
grc9  = [1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.50 1.55 1.65];
grc10 = [-0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1];
grc11 = [0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.700 0.750 0.850];

%--------------------------Generic Soil Period----------------------------%
period_soil = [0.0 0.075 0.10 0.20 0.30 0.40 0.50 0.75 1.00 1.50 2.00 3.00 4.00];

%-------------------------Generic Soil Coefficients-----------------------%
gsc1  = [-0.6687 1.7313 1.8473 0.8803 0.1243 -0.5247 -1.1067 -2.3727 -3.5387 -5.7697 -7.1017 -7.3407 -8.2867];
gsc2  = [1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438 1.438];
gsc3  = [0.0000 -0.0019 -0.0019 -0.0019 -0.0020 -0.0020 -0.0035 -0.0048 -0.0066 -0.0114 -0.0164 -0.0221 -0.0235];
gsc4  = [-2.329 -2.697 -2.697 -2.464 -2.327 -2.230 -2.140 -1.952 -1.785 -1.470 -1.290 -1.347 -1.272];
gsc5  = [0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648 0.00648];
gsc6  = [0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648 0.3648];
gsc7  = [1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970 1.0970];
gsc8  = [0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617 0.617];
gsc9  = [1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.45 1.50 1.55 1.65 1.65];
gsc10 = [-0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1 -0.1];
gsc11 = [0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.650 0.700 0.750 0.850 0.850];

% use the full period vector if T = 1000
if length (T) == 1 && T == 1000 && Zr ==0
    median=zeros(1,length(period_soil));
    sigma=zeros(1,length(period_soil));
    period1=period_soil;
elseif length (T) == 1 && T == 1000 && Zr ==1
    median=zeros(1,length(period_rock));
    sigma=zeros(1,length(period_rock));
    period1=period_rock;
else
    median=zeros(1,length(T));
    sigma=zeros(1,length(T));
    period1=T;
end

% iterate in the case that T is a vector
for n = 1:length(period1)
    
% Computation of Parameters
% For Sa Computation in Generic Soil
if Zr == 0
    
    if isempty(find(period_soil == period1(n), 1))
        
        i_low = sum(period_soil < period1(n));
        T_low = period_soil(i_low);
        T_high = period_soil(i_low + 1);
        
        [Sa_high, sigma_high,~] = ycsh_1997_subduction(T_high,rup,site,Zt);
        [Sa_low, sigma_low,~] = ycsh_1997_subduction(T_low,rup,site,Zt);
        
        x = [T_low T_high];
        Y_Sa = [Sa_low Sa_high];
        Y_sigma = [sigma_low sigma_high];
        median(n) = interp1(x,Y_Sa,period1(n));
        sigma(n) = interp1(x,Y_sigma,period1(n));
        
    else
        i = find(period_soil == period1(n));
        
        ln_Y = gsc1(i) + gsc2(i)*rup.M + gsc3(i)*(10-rup.M)^3 + gsc4(i)*(log((r)+...
            gsc7(i)*exp(gsc8(i)*rup.M))) + gsc5(i)*rup.Zhyp + gsc6(i)*Zt;
        
        sigma(n) = max((gsc9(i) + gsc10(i)*rup.M), gsc11(i));
        
        median(n) = exp(ln_Y);
        
    end
end

% For Sa Computation in Generic Rock
if Zr == 1
    
    if isempty(find(period_rock == period1(n), 1))
        
        i_low = sum(period_rock < period1(n));
        T_low = period_rock(i_low);
        T_high = period_rock(i_low + 1);
        
        [Sa_high, sigma_high,~] = ycsh_1997_subduction(T_high,rup,site,sub_ind);
        [Sa_low, sigma_low,~] = ycsh_1997_subduction(T_low, rup,site,sub_ind);
        
        x = [T_low T_high];
        Y_Sa = [Sa_low Sa_high];
        Y_sigma = [sigma_low sigma_high];
        median(n) = interp1(x,Y_Sa,period1(n));
        sigma(n) = interp1(x,Y_sigma,period1(n));
        
    else
        i = find(period_rock == period1(n));
        
        ln_Y = grc1(i) + grc2(i)*rup.M + grc3(i)*(10-rup.M)^3 + grc4(i)*(log((r)+...
            grc7(i)*exp(grc8(i)*rup.M))) + grc5(i)*rup.Zhyp + grc6(i)*Zt;
        
        sigma(n) = max((grc9(i) + grc10(i)*rup.M), grc11(i));
        
        median(n) = exp(ln_Y);
        
    end
end
end


function [median, sigma, period1] = aga_2016_subduction(T,rup,site,sub_ind,F_faba)

% Created by Reagan Chandramohan, circa 2017
% Modified by Jack Baker, 2/27/2019, to limit hypocentral depths to 120 km
% Updated by Emily Mongold, 11/27/20
% 
% To predict response spectra for subduction earthquakes
%
% Source Model:
% Abrahamson, N., Gregor, N., and Addo, K. (2016). "BC Hydro Ground Motion 
% Prediction Equations for Subduction Earthquakes." Earthquake Spectra, 
% 32(1), 23-44.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T            = vector of periods (s)
%   rup          = rupture object input containing the following
%                   variables:
%       M           = Moment magnitude
%       R           = Interface: Closest distance to rupture (km)
%                   = Intraslab: Distance to hypocenter (km)
%       Zhyp        = Hypocentral depth (km) (required only for intraslab events)
%   site         = site object input containing the following
%                   variable:
%       Vs30        = Average shear wave velocity over the top 30 m of the soil
%                   profile
%   sub_ind      = subduction indicator, reffered to locally as F_event
%                = 0 for interface
%                = 1 for intraslab events
%   F_faba       = 0 for forearc or unknown sites
%                = 1 for backarc sites
% OUTPUT
%   median          = vector of median response spectral ordinates (g)
%   sigma           = vector of logarithmic standard deviation
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F_event = sub_ind;

period = [0.000; 0.020; 0.050; 0.075; 0.100; 0.150; 0.200; 0.250; 0.300; 0.400; 0.500; 0.600; 0.750; 1.000; 1.500; 2.000; 2.500; 3.000; 4.000; 5.000; 6.000; 7.500; 10.000];
Vlin = [865.1; 865.1; 1053.5; 1085.7; 1032.5; 877.6; 748.2; 654.3; 587.1; 503.0; 456.6; 430.3; 410.5; 400.0; 400.0; 400.0; 400.0; 400.0; 400.0; 400.0; 400.0; 400.0; 400.0];
b = [-1.186; -1.186; -1.346; -1.471; -1.624; -1.931; -2.188; -2.381; -2.518; -2.657; -2.669; -2.599; -2.401; -1.955; -1.025; -0.299; 0.000; 0.000; 0.000; 0.000; 0.000; 0.000; 0.000];
n = 1.18;
c = 1.88;
C1 = 7.8;
C4 = 10;
theta1 = [4.2203; 4.2203; 4.5371; 5.0733; 5.2892; 5.4563; 5.2684; 5.0594; 4.7945; 4.4644; 4.0181; 3.6055; 3.2174; 2.7981; 2.0123; 1.4128; 0.9976; 0.6443; 0.0657; -0.4624; -0.9809; -1.6017; -2.2937];
theta2 = [-1.35; -1.35; -1.40; -1.45; -1.45; -1.45; -1.40; -1.35; -1.28; -1.18; -1.08; -0.99; -0.91; -0.85; -0.77; -0.71; -0.67; -0.64; -0.58; -0.54; -0.50; -0.46; -0.40];
theta3 = 0.1*ones(length(period), 1);
theta4 = 0.9*ones(length(period), 1);
theta5 = 0.0*ones(length(period), 1);
theta6 = [-0.0012; -0.0012; -0.0012; -0.0012; -0.0012; -0.0014; -0.0018; -0.0023; -0.0027; -0.0035; -0.0044; -0.0050; -0.0058; -0.0062; -0.0064; -0.0064; -0.0064; -0.0064; -0.0064; -0.0064; -0.0064; -0.0064; -0.0064];
theta7 = [1.0988; 1.0988; 1.2536; 1.4175; 1.3997; 1.3582; 1.1648; 0.9940; 0.8821; 0.7046; 0.5799; 0.5021; 0.3687; 0.1746; -0.0820; -0.2821; -0.4108; -0.4466; -0.4344; -0.4368; -0.4586; -0.4433; -0.4828];
theta8 = [-1.42; -1.42; -1.65; -1.80; -1.80; -1.69; -1.49; -1.30; -1.18; -0.98; -0.82; -0.70; -0.54; -0.34; -0.05; 0.12; 0.25; 0.30; 0.30; 0.30; 0.30; 0.30; 0.30];
theta9 = 0.4*ones(length(period), 1);
theta10 = [3.12; 3.12; 3.37; 3.37; 3.33; 3.25; 3.03; 2.80; 2.59; 2.20; 1.92; 1.70; 1.42; 1.10; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70; 0.70];
theta11 = [0.0130; 0.0130; 0.0130; 0.0130; 0.0130; 0.0130; 0.0129; 0.0129; 0.0128; 0.0127; 0.0125; 0.0124; 0.0120; 0.0114; 0.0100; 0.0085; 0.0069; 0.0054; 0.0027; 0.0005; -0.0013; -0.0033; -0.0060];
theta12 = [0.980; 0.980; 1.288; 1.483; 1.613; 1.882; 2.076; 2.248; 2.348; 2.427; 2.399; 2.273; 1.993; 1.470; 0.408; -0.401; -0.723; -0.673; -0.627; -0.596; -0.566; -0.528; -0.504];
theta13 = [-0.0135; -0.0135; -0.0138; -0.0142; -0.0145; -0.0153; -0.0162; -0.0172; -0.0183; -0.0206; -0.0231; -0.0256; -0.0296; -0.0363; -0.0493; -0.0610; -0.0711; -0.0798; -0.0935; -0.0980; -0.0980; -0.0980; -0.0980];
theta14 = [-0.40; -0.40; -0.40; -0.40; -0.40; -0.40; -0.35; -0.31; -0.28; -0.23; -0.19; -0.16; -0.12; -0.07; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00];
theta15 = [0.9996; 0.9996; 1.1030; 1.2732; 1.3042; 1.2600; 1.2230; 1.1600; 1.0500; 0.8000; 0.6620; 0.5800; 0.4800; 0.3300; 0.3100; 0.3000; 0.3000; 0.3000; 0.3000; 0.3000; 0.3000; 0.3000; 0.3000];
theta16 = [-1.00; -1.00; -1.18; -1.36; -1.36; -1.30; -1.25; -1.17; -1.06; -0.78; -0.62; -0.50; -0.34; -0.14; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00; 0.00];
phi = 0.6;
tau = 0.43;
sigma_ss = 0.6;

base = zeros(length(period), 1);
f_mag = zeros(length(period), 1);
f_dep = zeros(length(period), 1);
f_faba = zeros(length(period), 1);
f_site = zeros(length(period), 1);

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

if F_event
    deltaC1 = -0.3*ones(length(period), 1);
else
    deltaC1 = [0.2; interp1(log([1e-10; 0.3; 0.5; 1.0; 2.0; 3.0; 1e10]), [0.2; 0.2; 0.1; 0.0; -0.1; -0.2; -0.2], log(period(2:end)))];
end

if site.Vs30 > 1000
    Vs_star = 1000;
else
    Vs_star = site.Vs30;
end

for i = 1:length(period)
    base(i) = theta1(i) + theta4(i)*deltaC1(i) + (theta2(i) + theta14(i)*F_event + theta3(i)*(rup.M - C1))*log(rup.R + C4*exp((rup.M - 6)*theta9(i))) + theta6(i)*rup.R + theta10(i)*F_event;
    
    if rup.M <= C1 + deltaC1(i)
        f_mag(i) = theta4(i)*(rup.M - (C1 + deltaC1(i))) + theta13(i)*(10 - rup.M)^2;
    else
        f_mag(i) = theta5(i)*(rup.M - (C1 + deltaC1(i))) + theta13(i)*(10 - rup.M)^2;
    end
    
    if F_event
        f_dep(i) = theta11(i)*(min([120 rup.Zhyp]) - 60)*F_event; % modified by JWB, 2/27/2019, to include the "min" term, per equation 3 of the paper
    else
        f_dep(i) = 0;
    end
    
    if F_event
        f_faba(i) = (theta7(i) + theta8(i)*log(max(rup.R, 85)/40))*F_faba;
    else
        f_faba(i) = (theta15(i) + theta16(i)*log(max(rup.R, 100)/40))*F_faba;
    end
    
    if i == 1
        PGA1000 = exp(base(1) + f_mag(1) + f_dep(1) + f_faba(1) + theta12(1)*log(1000/Vlin(1)) + b(1)*n*log(1000/Vlin(1)));
    end
    if site.Vs30 < Vlin(i)
        f_site(i) = theta12(i)*log(Vs_star/Vlin(i)) - b(i)*log(PGA1000 + c) + b(i)*log(PGA1000 + c*(Vs_star/Vlin(i))^n);
    else
        f_site(i) = theta12(i)*log(Vs_star/Vlin(i)) + b(i)*n*log(Vs_star/Vlin(i));
    end
end

for n = 1:length(period1)
    if period1(n) == 0
        sa = exp(base + f_mag + f_dep + f_faba + f_site);
        median(n) = sa(1);
        sigma(n) = sqrt(phi.^2 + tau.^2);
    else
        sa = exp(base + f_mag + f_dep + f_faba + f_site);
        median(n) = interp1(log([1e-10; period(2:end)]), sa, log(period1(n)));
        sigma(n) = sqrt(phi.^2 + tau.^2);
    end
end
end

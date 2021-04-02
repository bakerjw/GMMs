function [median, sigma, period1] = ask_2014_active(T,rup,site)

% Created by Yue Hua, 5/19/10, yuehua@stanford.edu
% Edited by Jack Baker, 9/22/2017 to fix an error in FVs30 variable
% definition in the comments (no change to the function of the code).
% Updated by Emily Mongold, 11/27/20
%
% Source Model:
% Abrahamson, N. A., Silva, W. J., and Kamai, R. (2014). "Summary of the 
% ASK14 Ground Motion Relation for Active Crustal Regions." Earthquake 
% Spectra, 30(3), 1025-1055.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = Period (sec); Use Period = -1 for PGV computation
%                   Use 1000 for output the array of Sa with original period
%   rup           = rupture object input containing the following
%                   variables:
%       M             = Moment Magnitude
%       Rrup          = Closest distance (km) to the ruptured plane
%       Rjb           = Joyner-Boore distance (km); closest distance (km) to surface
%                       projection of rupture plane
%       Rx            =  Site coordinate (km) measured perpendicular to the fault strike
%                       from the fault line with down-dip direction to be positive
%       Ry0           = Horizontal distance off the end of the rupture measured parallel
%                       to strike
%       Ztor          = Depth(km) to the top of ruptured plane
%       AS            = Flag for aftershocks
%       HW            = Flag for hanging wall sites
%       W             = Down-dip rupture width (km) 
%       delta         = Fault dip angle (in degrees)
%       lambda        = Rake angle      (in degrees)
%   site          = site object input containing the following
%                   variables:
%       region        = 0 for global
%                     = 1 for California
%                     = 2 for Japan
%                     = 3 for China 
%                     = 4 for Italy 
%                     = 5 for Turkey
%                     = 6 for Taiwan
%       Z10           = Basin depth (km); depth from the groundsurface to the
%                       1km/s shear-wave horizon.
%       Vs30          = shear wave velocity averaged over top 30 m in m/s
%                     = ref: 1130
%       fvs30         = 1 for measured Vs30 
%                     = 0 for Vs30 inferred from geology
% OUTPUT
%   median     = Median spectral acceleration prediction
%   sigma      = logarithmic standard deviation of spectral acceleration
%                prediction
%   period1    = periods for which the median and sigma values are
%                provided. If T = 1000, then period1 = the full set of
%                available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0	-1];

frv = rup.lambda >= 30 & rup.lambda <= 150; % frv: 1 for lambda between 30 and 150, 0 otherwise
fnm = rup.lambda >= -150 & rup.lambda <= -30; % fnm: 1 for lambda between -150 and -30, 0 otherwise

if isempty(rup.Ztor)
    if frv == 1
        rup.Ztor = max(2.704 - 1.226 * max(rup.M-5.849,0),0)^2;
    else
        rup.Ztor = max(2.673 - 1.136 * max(rup.M-4.970,0),0)^2;
    end
end
W_empty_flag = 0;
if isempty(rup.W)
    W_empty_flag = 1; % To set rup.W = [] at the end
    rup.W = min(18/sin(deg2rad(rup.delta)),10^(-1.75+0.45*rup.M));
end

if isempty(site.Z10)
    if site.region == 2
        Z10 = exp(-5.23 / 2 * log((site.Vs30 ^ 2 + 412 ^ 2) / (1360 ^ 2 + 412 ^ 2))) / 1000;
    else % not Japan
        Z10 = exp((-7.67 / 4) * log((site.Vs30 ^ 4 + 610 ^ 4) / (1360 ^ 4 + 610 ^ 4))) / 1000;
    end
else
    Z10 = site.Z10;
end

if length (T) == 1 && T == 1000 % Compute median and sigma with pre-defined period
    median=zeros(1,length(period)-2);
    sigma=zeros(1,length(period)-2);
    period1=period(1:end-2);
    for ip=1:length(period)-2
        [median(ip),sigma(ip)]=ASK_2014_sub_1(rup.M,ip,rup.Rrup,rup.Rjb,rup.Rx,rup.Ry0,rup.Ztor,rup.delta,frv,fnm,rup.AS,rup.HW,rup.W,Z10,site.Vs30,site.fvs30,site.region);
    end
else                            % Compute median and sigma with user-defined period
    median=zeros(1, length(T));
    sigma=zeros(1, length(T));
    period1=T;
    for i=1:length(T)
        Ti = T(i);
        if (isempty(find(abs(period-Ti) < 0.0001, 1))) % The user defined period requires interpolation
            T_low = max(period(period < Ti));
            T_high = min(period(period > Ti));
            ip_low  = find(period==T_low);
            ip_high = find(period==T_high);
            
            [Sa_low, sigma_low] = ASK_2014_sub_1(rup.M,ip_low,rup.Rrup,rup.Rjb,rup.Rx,rup.Ry0,rup.Ztor,rup.delta,frv,fnm,rup.AS,rup.HW,rup.W,Z10,site.Vs30,site.fvs30,site.region);
            [Sa_high, sigma_high] = ASK_2014_sub_1(rup.M,ip_high,rup.Rrup,rup.Rjb,rup.Rx,rup.Ry0,rup.Ztor,rup.delta,frv,fnm,rup.AS,rup.HW,rup.W,Z10,site.Vs30,site.fvs30,site.region);
            x = [log(T_low) log(T_high)];
            Y_sa = [log(Sa_low) log(Sa_high)];
            Y_sigma = [sigma_low sigma_high];
            median(i) = exp(interp1(x, Y_sa, log(Ti)));
            sigma(i) = interp1(x, Y_sigma, log(Ti));
        else
            ip_T = find(abs((period- Ti)) < 0.0001);
            [median(i), sigma(i)] = ASK_2014_sub_1(rup.M,ip_T,rup.Rrup,rup.Rjb,rup.Rx,rup.Ry0,rup.Ztor,rup.delta,frv,fnm,rup.AS,rup.HW,rup.W,Z10,site.Vs30,site.fvs30,site.region);
        end
    end
    if W_empty_flag == 1
        rup.W = [];
    end
end



function [Sa, sigma]=ASK_2014_sub_1 (M, ip, R_RUP, R_JB, Rx, Ry0, Ztor, delta, F_RV, F_NM, F_AS, HW, W, Z10, Vs30, FVS30, region)

%% Coefficients
T = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	6	7.5	10	0	-1];
Vlin =[	660.0000	680.0000	770.0000	915.0000	960.0000	910.0000	740.0000	590.0000	495.0000	430.0000	360.0000	340.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	330.0000	660.0000	330.0000];
b	=[-1.4700	-1.4590	-1.3900	-1.2190	-1.1520	-1.2300	-1.5870	-2.0120	-2.4110	-2.7570	-3.2780	-3.5990	-3.8000	-3.5000	-2.4000	-1.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	-1.4700	-2.0200];
n	=[1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000	1.5000];
M1	=[6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.7500	6.8200	6.9200	7.0000	7.0600	7.1450	7.2500	6.7500	6.7500];
c	=[2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2.4000	2400.0000];
c4	=[4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000	4.5000];
a1	=[0.5870	0.5980	0.6020	0.7070	0.9730	1.1690	1.4420	1.6370	1.7010	1.7120	1.6620	1.5710	1.2990	1.0430	0.6650	0.3290	-0.0600	-0.2990	-0.5620	-0.8750	-1.3030	-1.9280	0.5870	5.9750];
a2	=[-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7900	-0.7650	-0.7110	-0.6340	-0.5290	-0.7900	-0.9190];
a3	=[0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750	0.2750];
a4	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000];
a5	=[-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100	-0.4100];
a6	=[2.1541	2.1461	2.1566	2.0845	2.0285	2.0408	2.1208	2.2241	2.3124	2.3383	2.4688	2.5586	2.6821	2.7630	2.8355	2.8973	2.9061	2.8888	2.8984	2.8955	2.8700	2.8431	2.1541	2.3657];
a7	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a8	=[-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0150	-0.0220	-0.0300	-0.0380	-0.0450	-0.0550	-0.0650	-0.0950	-0.1100	-0.1240	-0.1380	-0.1720	-0.1970	-0.2180	-0.2350	-0.2550	-0.2850	-0.0150	-0.0940];
a10	=[1.7350	1.7180	1.6150	1.3580	1.2580	1.3100	1.6600	2.2200	2.7700	3.2500	3.9900	4.4500	4.7500	4.3000	2.6000	0.5500	-0.9500	-0.9500	-0.9300	-0.9100	-0.8700	-0.8000	1.7350	2.3600];
a11	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a12	=[-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.1000	-0.2000	-0.2000	-0.2000	-0.1000	-0.1000];
a13	=[0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.6000	0.5800	0.5600	0.5300	0.5000	0.4200	0.3500	0.2000	0.0000	0.0000	0.0000	0.0000	0.0000	0.6000	0.2500];
a14	=[-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.3000	-0.2400	-0.1900	-0.1100	-0.0400	0.0700	0.1500	0.2700	0.3500	0.4600	0.5400	0.6100	0.6500	0.7200	0.8000	-0.3000	0.2200];
a15	=[1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.1000	1.0300	0.9200	0.8400	0.6800	0.5700	0.4200	0.3100	0.1600	0.0500	-0.0400	-0.1100	-0.1900	-0.3000	1.1000	0.3000];
a17	=[-0.0072	-0.0073	-0.0075	-0.0080	-0.0089	-0.0095	-0.0095	-0.0086	-0.0074	-0.0064	-0.0043	-0.0032	-0.0025	-0.0025	-0.0022	-0.0019	-0.0015	-0.0010	-0.0010	-0.0010	-0.0010	-0.0010	-0.0072	-0.0005];
a43	=[0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1000	0.1400	0.1700	0.2200	0.2600	0.3400	0.4100	0.5100	0.5500	0.4900	0.4200	0.1000	0.2800];
a44	=[0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0500	0.0700	0.1000	0.1400	0.1700	0.2100	0.2500	0.3000	0.3200	0.3200	0.3200	0.2750	0.2200	0.0500	0.1500];
a45 =[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0300	0.0600	0.1000	0.1400	0.1700	0.2000	0.2200	0.2300	0.2300	0.2200	0.2000	0.1700	0.1400	0.0000	0.0900];
a46	=[-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0500	-0.0300	0.0000	0.0300	0.0600	0.0900	0.1300	0.1400	0.1600	0.1600	0.1600	0.1400	0.1300	0.1000	0.0900	0.0800	-0.0500	0.0700];
a25	=[-0.0015	-0.0015	-0.0016	-0.0020	-0.0027	-0.0033	-0.0035	-0.0033	-0.0029	-0.0027	-0.0023	-0.0020	-0.0010	-0.0005	-0.0004	-0.0002	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	-0.0015	-0.0001];
a28	=[0.0025	0.0024	0.0023	0.0027	0.0032	0.0036	0.0033	0.0027	0.0024	0.0020	0.0010	0.0008	0.0007	0.0007	0.0006	0.0003	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0025	0.0005];
a29	=[-0.0034	-0.0033	-0.0034	-0.0033	-0.0029	-0.0025	-0.0025	-0.0031	-0.0036	-0.0039	-0.0048	-0.0050	-0.0041	-0.0032	-0.0020	-0.0017	-0.0020	-0.0020	-0.0020	-0.0020	-0.0020	-0.0020	-0.0034	-0.0037];
a31	=[-0.1503	-0.1479	-0.1447	-0.1326	-0.1353	-0.1128	0.0383	0.0775	0.0741	0.2548	0.2136	0.1542	0.0787	0.0476	-0.0163	-0.1203	-0.2719	-0.2958	-0.2718	-0.2517	-0.1400	-0.0216	-0.1503	-0.1462];
a36	=[0.2650	0.2550	0.2490	0.2020	0.1260	0.0220	-0.1360	-0.0780	0.0370	-0.0910	0.1290	0.3100	0.5050	0.3580	0.1310	0.1230	0.1090	0.1350	0.1890	0.2150	0.1500	0.0920	0.2650	0.3770];
a37	=[0.3370	0.3280	0.3200	0.2890	0.2750	0.2560	0.1620	0.2240	0.2480	0.2030	0.2320	0.2520	0.2080	0.2080	0.1080	0.0680	-0.0230	0.0280	0.0310	0.0240	-0.0700	-0.1590	0.3370	0.2120];
a38	=[0.1880	0.1840	0.1800	0.1670	0.1730	0.1890	0.1080	0.1150	0.1220	0.0960	0.1230	0.1340	0.1290	0.1520	0.1180	0.1190	0.0930	0.0840	0.0580	0.0650	0.0000	-0.0500	0.1880	0.1570];
a39	=[0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000	0.0000];
a40	=[0.0880	0.0880	0.0930	0.1330	0.1860	0.1600	0.0680	0.0480	0.0550	0.0730	0.1430	0.1600	0.1580	0.1450	0.1310	0.0830	0.0700	0.1010	0.0950	0.1330	0.1510	0.1240	0.0880	0.0950];
a41	=[-0.1960	-0.1940	-0.1750	-0.0900	0.0900	0.0060	-0.1560	-0.2740	-0.2480	-0.2030	-0.1540	-0.1590	-0.1410	-0.1440	-0.1260	-0.0750	-0.0210	0.0720	0.2050	0.2850	0.3290	0.3010	-0.1960	-0.0380];
a42	=[0.0440	0.0610	0.1620	0.4510	0.5060	0.3350	-0.0840	-0.1780	-0.1870	-0.1590	-0.0230	-0.0290	0.0610	0.0620	0.0370	-0.1430	-0.0280	-0.0970	0.0150	0.1040	0.2990	0.2430	0.0440	0.0650];
s1	=[0.7540	0.7600	0.7810	0.8100	0.8100	0.8100	0.8010	0.7890	0.7700	0.7400	0.6990	0.6760	0.6310	0.6090	0.5780	0.5550	0.5480	0.5270	0.5050	0.4770	0.4570	0.4290	0.7540	0.6620];
s2	=[0.5200	0.5200	0.5200	0.5300	0.5400	0.5500	0.5600	0.5650	0.5700	0.5800	0.5900	0.6000	0.6150	0.6300	0.6400	0.6500	0.6400	0.6300	0.6300	0.6300	0.6300	0.6300	0.5200	0.5100];
s3	=[0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.4700	0.3800];
s4	=[0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3600	0.3800];
s1_m =[0.7410	0.7470	0.7690	0.7980	0.7980	0.7950	0.7730	0.7530	0.7290	0.6930	0.6440	0.6160	0.5660	0.5410	0.5060	0.4800	0.4720	0.4470	0.4250	0.3950	0.3780	0.3590	0.7410	0.6600];
s2_m =[0.5010	0.5010	0.5010	0.5120	0.5220	0.5270	0.5190	0.5140	0.5130	0.5190	0.5240	0.5320	0.5480	0.5650	0.5760	0.5870	0.5760	0.5650	0.5680	0.5710	0.5750	0.5850	0.5010	0.5100];
s5_JP=[0.5400	0.5400	0.5500	0.5600	0.5700	0.5700	0.5800	0.5900	0.6100	0.6300	0.6600	0.6900	0.7300	0.7700	0.8000	0.8000	0.8000	0.7600	0.7200	0.7000	0.6700	0.6400	0.5400	0.5800];
s6_JP=[0.6300	0.6300	0.6300	0.6500	0.6900	0.7000	0.7000	0.7000	0.7000	0.7000	0.7000	0.7000	0.6900	0.6800	0.6600	0.6200	0.5500	0.5200	0.5000	0.5000	0.5000	0.5000	0.6300	0.5300];

M2 = 5;
CRjb = 999.9; %   CRjb   = Centroid CRjb, assumed to be 999.9 here -> assume no aftershock

%% Term f1 - Basic form
if M > 5
    c4m = c4(ip);
elseif M > 4 && M <= 5
    c4m = c4(ip)-(c4(ip)-1)*(5-M);
else
    c4m = 1;
end

R = sqrt(R_RUP^2+ c4m^2);

if M > M1(ip)
    f1 = a1(ip) + a5(ip)*(M - M1(ip)) + a8(ip)*(8.5 - M)^2 + (a2(ip) + a3(ip)*(M - M1(ip)))*log(R) + a17(ip)*R_RUP;
elseif M >= M2 && M <= M1(ip)
    f1 = a1(ip) + a4(ip)*(M - M1(ip)) + a8(ip)*(8.5 - M)^2 + (a2(ip) + a3(ip)*(M - M1(ip)))*log(R) + a17(ip)*R_RUP;
else
    f1 = a1(ip) + a4(ip)*(M2 - M1(ip)) + a8(ip)*(8.5 - M2)^2 + a6(ip)*(M - M2) + a7(ip)*(M - M2)^2 + ...
        (a2(ip) + a3(ip)*(M2 - M1(ip)))*log(R) + a17(ip)*R_RUP;
end

%% term f4 - Hanging wall model
R1 = W * cos(deg2rad(delta));
R2 = 3 * R1;
Ry1 = Rx * tan(deg2rad(20));
h1 = 0.25;
h2 = 1.5;
h3 = -0.75;

if delta > 30
    T1 = (90- delta)/45;
else
    T1 = 60/45;
end

a2hw = 0.2;

if M > 6.5
    T2 = 1 + a2hw * (M - 6.5);
elseif M > 5.5
    T2 = 1 + a2hw * (M - 6.5) - (1 - a2hw) * (M - 6.5)^2;
else
    T2 = 0;
end

if Rx <= R1
    T3 = h1 + h2*(Rx/R1) + h3*(Rx/R1)^2;
elseif Rx < R2
    T3 = 1 - (Rx - R1)/(R2 - R1);
else 
    T3 = 0;
end

if Ztor < 10
    T4 = 1 - Ztor^2/100;
else
    T4 = 0;
end

if isempty(Ry0) || Ry0 == 0
    if R_JB == 0
        T5 = 1;
    elseif R_JB <30
        T5 = 1 - R_JB/30;
    else
        T5 = 0;
    end
else
    if Ry0 - Ry1 <= 0
        T5 = 1;
    elseif Ry0 - Ry1 < 5
        T5 = 1- (Ry0-Ry1)/5;
    else
        T5 = 0;
    end
end

if HW == 1
    f4 = a13(ip) * T1 * T2 * T3 * T4 * T5;
else
    f4 = 0;
end

%% Term f6 - Depth to top rupture model

if Ztor < 20
    f6 = a15(ip)*Ztor/20;
else
    f6 = a15(ip);
end

%% Term: f7 and f8 - Style of Faulting

if M > 5
    f7 = a11(ip);
    f8 = a12(ip);
elseif M >= 4 && M <= 5
    f7 = a11(ip)*(M - 4);
    f8 = a12(ip)*(M - 4);
else 
    f7 = 0;
    f8 = 0;
end

if T(ip) <= 0.5
    V1 = 1500;
elseif T(ip) < 3
    V1 = exp(-0.35*log(T(ip)/0.5)+log(1500));
else 
    V1 = 800;
end  

if Vs30 < V1
    Vs30s = Vs30;
else
    Vs30s = V1;
end

if 1180 >= V1
    Vs30star1180 = V1;
else
    Vs30star1180 = 1180;
end

%% term  Regional:
Ftw = (region == 6);
Fcn = (region == 3);
Fjp = (region == 2);

% Japan
if Vs30 < 150 
    y1 = a36(ip);
    y2 = a36(ip);
    x1 = 50;
    x2 = 150;
elseif Vs30 < 250 
    y1 = a36(ip);
    y2 = a37(ip);
    x1 = 150;
    x2 = 250;
elseif Vs30 < 350 
    y1 = a37(ip);
    y2 = a38(ip);
    x1 = 250;
    x2 = 350;
elseif Vs30 < 450
    y1 = a38(ip);
    y2 = a39(ip);
    x1 = 350;
    x2 = 450;
elseif Vs30 < 600 
    y1 = a39(ip);
    y2 = a40(ip);
    x1 = 450;
    x2 = 600;
elseif Vs30 < 850
    y1 = a40(ip);
    y2 = a41(ip);
    x1 = 600;
    x2 = 850;
elseif Vs30 < 1150
    y1 = a41(ip);
    y2 = a42(ip);
    x1 = 850;
    x2 = 1150;
else
    y1 = a42(ip);
    y2 = a42(ip);
    x1 = 1150;
    x2 = 3000;
end
    f13Vs30 = y1 + (y2 - y1) / (x2 - x1) * (Vs30 - x1);

% Taiwan
    f12Vs30 = a31(ip) * log(Vs30s/Vlin(ip));
    f12Vs30_1180 = a31(ip) * log(Vs30star1180/Vlin(ip));

Regional = Ftw * (f12Vs30 + a25(ip) * R_RUP) + Fcn * (a28(ip) * R_RUP) + Fjp * (f13Vs30 + a29(ip) * R_RUP);
Regional_1180 = Ftw * (f12Vs30_1180 + a25(ip) * R_RUP) + Fcn * (a28(ip) * R_RUP) + Fjp * (f13Vs30 + a29(ip) * R_RUP);

%% Term f5 - site response model

%Sa 1180
f5_1180 = (a10(ip) + b(ip) * n(ip)) * log(Vs30star1180 / Vlin(ip));

Sa1180 = exp(f1 + f6 + F_RV*f7 + F_NM * f8 +  HW * f4 + f5_1180 + Regional_1180);

if Vs30 >= Vlin(ip)
    f5 = (a10(ip)+ b(ip)*n(ip))*log(Vs30s/Vlin(ip));
else 
    f5 = a10(ip)*log(Vs30s/Vlin(ip)) - b(ip)*log(Sa1180 + c(ip)) + b(ip)*log(Sa1180 + c(ip)*(Vs30s/Vlin(ip))^n(ip));
end

%% Term f10 - soil depth model

if region ~= 2 % california
        Z1ref = 1/1000* exp(-7.67/4 * log((Vs30^4 + 610^4)/(1360^4 + 610^4)));
else  % Japan
        Z1ref =  1/1000* exp(-5.23/2 * log((Vs30^2 + 412^2)/(1360^2 + 412^2)));
end

if Vs30 <= 150
    y1z = a43(ip);
    y2z = a43(ip);
    x1z = 50;
    x2z = 150;
    elseif Vs30 <= 250
    y1z = a43(ip);
    y2z = a44(ip);
    x1z = 150;
    x2z = 250;
    elseif Vs30 <= 400 
    y1z = a44(ip);
    y2z = a45(ip);
    x1z = 250;
    x2z = 400;
    elseif Vs30 <= 700 
    y1z = a45(ip);
    y2z = a46(ip);
    x1z = 400;
    x2z = 700;
    else
    y1z = a46(ip);
    y2z = a46(ip);
    x1z = 700;
    x2z = 1000;
end

%'f10 term goes to zero at 1180 m/s (reference)
f10 = (y1z + (Vs30 - x1z) * (y2z - y1z) / (x2z - x1z)) * log((Z10 + 0.01) / (Z1ref + 0.01));
    
%% term f11 - Aftershock scaling

if CRjb <= 5
    f11 = a14(ip);
elseif CRjb < 15
    f11 = a14(ip)*(1-(CRjb-5)/10);
else
    f11 = 0;
end

if F_AS == 0
    f11 = 0;
end

%% Sa

lnSa = f1+ f6 + F_RV*f7 + F_NM*f8 +  HW*f4 + F_AS*f11 + f5 + f10 + Regional;
 
Sa= exp(lnSa);
 
 
%% Standard deviation 

if FVS30 == 1 % measured
    s1 = s1_m;
    s2 = s2_m;
end

if M < 4
    phi_AL = s1(ip);
elseif M <= 6
    phi_AL = s1(ip)+(s2(ip)-s1(ip))/2*(M-4);
else
    phi_AL = s2(ip);
end

if M < 5
    tau_AL = s3(ip);
elseif M <= 7
    tau_AL = s3(ip) + (s4(ip)-s3(ip))/2*(M-5);
else 
    tau_AL = s4(ip);
end
tau_B = tau_AL;

if Fjp == 1
    if R_RUP < 30 
        phi_AL = s5_JP(ip);
    elseif R_RUP <= 80
        phi_AL = s5_JP(ip) + (s6_JP(ip) - s5_JP(ip))/50 * (R_RUP - 30);
    else 
        phi_AL = s6_JP(ip);
    end
end

phi_amp = 0.4;

phi_B = sqrt(phi_AL^2 - phi_amp^2);

if Vs30 >= Vlin(ip)
    dln = 0;
else
    dln = -b(ip)*Sa1180/(Sa1180 + c(ip)) + b(ip)*Sa1180/(Sa1180 + c(ip)*(Vs30/Vlin(ip))^n(ip));
end

phi = sqrt(phi_B^2 * (1 + dln)^2 + phi_amp^2);

tau = tau_B*(1+ dln);

sigma = sqrt(phi^2+ tau^2);



function [median, sigma, period1] = ba_2008_active(T, rup, site)

% Created by Yoshifumi Yamamoto, 11/10/08, yama4423@stanford.edu
%   updated 2010/05/20
%   updated 2009/05/05
% Updated by Emily Mongold, 11/27/20
%
% Source Model: 
% Boore D. M. and Atkinson, G. M.(2008). "Ground-Motion Prediction
% Equations for the Average Horizontal Component of PGA, PGV, and 5%-Damped
% PSA at Spectral Periods between 0.01s and 10.0s." Earthquake Spectra,
% 24(1), 99-138. 
% 
% This script has been modified to correct an error based on the website
% <http://www.daveboore.com/pubs_online.php>.
% Table 6 should read "Distance-scaling coefficients (Mref=4.5 and Rref=1.0
% km for all periods)"
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T            = Period (sec); Use Period = -1 for PGV computation
%                  Use 1000 for output the array of Sa with period
%   rup          = rupture object input containing the following
%                   variables:
%       M           = Moment Magnitude
%       Rjb         = Joyner-Boore distance (km)
%       lambda      = rake angle, used to set FaultType:
%                       = 1 for unspecified fault 
%                       = 2 for strike-slip fault
%                       = 3 for normal fault
%                       = 4 for reverse fault
%   site         = site object input containing the following
%                   variable:
%       Vs30        = shear wave velocity averaged over top 30 m in m/s
% OUTPUT
%   median      = Median spectral acceleration prediction
%   sigma       = logarithmic standard deviation of spectral acceleration
%                 prediction
%   period1     = periods for which the median and sigma values are
%                 provided. If T = 1000, then period1 = the full set of
%                 available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Setting fault type based on rake angle
if rup.lambda <= 30 && rup.lambda >= -30 || rup.lambda <= 180 && rup.lambda >= 150 || rup.lambda <= -150 && rup.lambda >= -180  
    FaultType = 2; %Strike-slip
elseif rup.lambda <= -60 && rup.lambda >= -120
    FaultType = 3; %Normal
elseif rup.lambda <= 120 && rup.lambda >= 60
    FaultType = 4; %Reverse
else
    FaultType = 1; %Other/unspecified
end   

% Coefficients
period = [-1	0	0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10];

U = (FaultType == 1);
S = (FaultType == 2);
N = (FaultType == 3);
R = (FaultType == 4);

nT=length(T);
iflg=0;
if(nT==1)
    if(T==1000)
        iflg=1;
%     Computation with original periods
        nperi=length(period);
        median=zeros(1,nperi-2);
        sigma=zeros(1,nperi-2);
        period1=period(3:nperi);
        for i=3:1:nperi
            [median(i-2), sigma(i-2)] = BA_2008_nga_sub (rup.M, i, rup.Rjb, U,S,N,R, site.Vs30);
        end
    end
end
%     Computation with input periods
if(iflg==0)
    median=zeros(1,nT);
    sigma=zeros(1,nT);
    period1=T;
    for it=1:1:nT
        Teach=T(it);
        % interpolate between periods if neccesary    
        if (isempty(find(abs((period - Teach)) < 0.0001)))
            T_low = max(period(period<Teach));
            T_hi  = min(period(period>Teach));

            [sa_low, sigma_low] = ba_2008_active(T_low,rup,site);
            [sa_hi,  sigma_hi]  = ba_2008_active(T_hi,rup,site);

            x = [log(T_low) log(T_hi)];
            Y_sa = [log(sa_low) log(sa_hi)];
            Y_sigma = [sigma_low sigma_hi];
            median(it) = exp(interp1(x,Y_sa,log(Teach)));
            sigma(it) = interp1(x,Y_sigma,log(Teach));
        else
            i = find(abs((period - Teach)) < 0.0001); % Identify the period
            [median(it), sigma(it)] = BA_2008_nga_sub (rup.M, i,rup.Rjb, U,S,N,R,site.Vs30);
        end
    end
end 
    
function [median, sigma] = BA_2008_nga_sub (M, ip, Rjb, U,S,N,R, Vs30)
        
% Coefficients
e01 = [5.00121	-0.53804	-0.52883	-0.52192	-0.45285	-0.28476	0.00767	0.20109	0.46128	0.5718	0.51884	0.43825	0.3922	0.18957	-0.21338	-0.46896	-0.86271	-1.22652	-1.82979	-2.24656	-1.28408	-1.43145	-2.15446];
e02 = [5.04727	-0.5035	-0.49429	-0.48508	-0.41831	-0.25022	0.04912	0.23102	0.48661	0.59253	0.53496	0.44516	0.40602	0.19878	-0.19496	-0.43443	-0.79593	-1.15514	-1.7469	-2.15906	-1.2127	-1.31632	-2.16137];
e03 = [4.63188	-0.75472	-0.74551	-0.73906	-0.66722	-0.48462	-0.20578	0.03058	0.30185	0.4086	0.3388	0.25356	0.21398	0.00967	-0.49176	-0.78465	-1.20902	-1.57697	-2.22584	-2.58228	-1.50904	-1.81022	-2.53323];
e04 = [5.0821	-0.5097	-0.49966	-0.48895	-0.42229 -0.26092	0.02706	0.22193	0.49328	0.61472	0.57747	0.5199	0.4608	0.26337	-0.10813	-0.3933	-0.88085	-1.27669	-1.91814	-2.38168	-1.41093	-1.59217	-2.14635];
e05 = [0.18322	0.28805	0.28897	0.25144	0.17976	0.06369	0.0117	0.04697	0.1799	0.52729	0.6088	0.64472	0.7861	0.76837	0.75179	0.6788	0.70689	0.77989	0.77966	1.24961	0.14271	0.52407	0.40387];
e06 = [-0.12736	-0.10164	-0.10019	-0.11006	-0.12858	-0.15752	-0.17051	-0.15948	-0.14539	-0.12964	-0.13843	-0.15694	-0.07843	-0.09054	-0.14053	-0.18257	-0.2595	-0.29657	-0.45384	-0.35874	-0.39006	-0.37578	-0.48492];
e07 = [0	0	0	0	0	0	0	0	0	0.00102	0.08607	0.10601	0.02262	0	0.10302	0.05393	0.19082	0.29888	0.67466	0.79508	0	0	0];
mh = [8.5	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	6.75	8.5	8.5	8.5];
c01 = [-0.8737	-0.6605	-0.6622	-0.666	-0.6901	-0.717	-0.7205	-0.7081	-0.6961	-0.583	-0.5726	-0.5543	-0.6443	-0.6914	-0.7408	-0.8183	-0.8303	-0.8285	-0.7844	-0.6854	-0.5096	-0.3724	-0.09824];
c02 = [0.1006	0.1197	0.12	0.1228	0.1283	0.1317	0.1237	0.1117	0.09884	0.04273	0.02977	0.01955	0.04394	0.0608	0.07518	0.1027	0.09793	0.09432	0.07282	0.03758	-0.02391	-0.06568	-0.138];
c03 = [-0.00334	-0.01151	-0.01151	-0.01151	-0.01151	-0.01151	-0.01151	-0.01151	-0.01113	-0.00952	-0.00837	-0.0075	-0.00626	-0.0054	-0.00409	-0.00334	-0.00255	-0.00217	-0.00191	-0.00191	-0.00191	-0.00191	-0.00191];
h = [2.54	1.35	1.35	1.35	1.35	1.35	1.55	1.68	1.86	1.98	2.07	2.14	2.24	2.32	2.46	2.54	2.66	2.73	2.83	2.89	2.93	3	3.04];

blin = [-0.6	-0.36	-0.36	-0.34	-0.33	-0.29	-0.23	-0.25	-0.28	-0.31	-0.39	-0.44	-0.5	-0.6	-0.69	-0.7	-0.72	-0.73	-0.74	-0.75	-0.75	-0.692	-0.65];
b1 = [-0.5	-0.64	-0.64	-0.63	-0.62	-0.64	-0.64	-0.6	-0.53	-0.52	-0.52	-0.52	-0.51	-0.5	-0.47	-0.44	-0.4	-0.38	-0.34	-0.31	-0.291	-0.247	-0.215];
b2 = [-0.06	-0.14	-0.14	-0.12	-0.11	-0.11	-0.11	-0.13	-0.18	-0.19	-0.16	-0.14	-0.1	-0.06	0	0	0	0	0	0	0	0	0];
sig1 = [0.5	0.502	0.502	0.502	0.507	0.516	0.513	0.52	0.518	0.523	0.527	0.546	0.541	0.555	0.571	0.573	0.566	0.58	0.566	0.583	0.601	0.626	0.645];
sig2u = [0.286	0.265	0.267	0.267	0.276	0.286	0.322	0.313	0.288	0.283	0.267	0.272	0.267	0.265	0.311	0.318	0.382	0.398	0.41	0.394	0.414	0.465	0.355];
sigtu = [0.576	0.566	0.569	0.569	0.578	0.589	0.606	0.608	0.592	0.596	0.592	0.608	0.603	0.615	0.649	0.654	0.684	0.702	0.7	0.702	0.73	0.781	0.735];
sig2m = [0.256	0.26	0.262	0.262	0.274	0.286	0.32	0.318	0.29	0.288	0.267	0.269	0.267	0.265	0.299	0.302	0.373	0.389	0.401	0.385	0.437	0.477	0.477];
sigtm = [0.56	0.564	0.566	0.566	0.576	0.589	0.606	0.608	0.594	0.596	0.592	0.608	0.603	0.615	0.645	0.647	0.679	0.7	0.695	0.698	0.744	0.787	0.801];
a1 = 0.03;
pga_low = 0.06;
a2 = 0.09;
v1 = 180;
v2 = 300;
vref = 760;

mref=4.5;
rref=1;

% Magnitude Scaling

if M <= mh (ip)
    Fm = e01(ip) * U + e02(ip) * S + e03(ip) * N + e04(ip) * R + e05(ip) * (M - mh(ip)) + e06(ip) * (M - mh(ip))^2;
else
    Fm = e01(ip) * U + e02(ip) * S + e03(ip) * N + e04(ip) * R + e07(ip) * (M - mh(ip));
end

% Distance Scaling

r = sqrt (Rjb^2 + h(ip)^2);
Fd = (c01(ip) + c02(ip) * (M - mref)) * log (r / rref) + c03(ip) * (r - rref);

if Vs30~=vref || ip~=2
    pga4nl = BA_2008_nga_sub (M, 2, Rjb, U,S,N,R, vref);
else
    % Compute median and sigma
    lny = Fm + Fd;
    median = exp (lny); 
    sigma = (U==1) * sigtu(2) + (U~=1) * sigtm(2);
    return;
end

% Site Amplification

% Linear term
Flin = blin(ip) * log (Vs30 / vref);

% Nonlinear term

% Computation of nonlinear factor
if Vs30 <= v1
    bnl = b1(ip);
elseif Vs30 <= v2
    bnl = b2(ip) + (b1(ip) - b2(ip)) * log (Vs30 / v2) / log (v1 / v2);
elseif Vs30 <= vref
    bnl = b2(ip) * log (Vs30 / vref) / log (v2 / vref);
else
    bnl = 0.0;
end
deltax = log(a2/a1);
deltay = bnl * log(a2/pga_low);
c = (3 *deltay - bnl * deltax) / (deltax^2);
d = - (2 * deltay - bnl * deltax) / (deltax^3);

if pga4nl <= a1
    Fnl = bnl * log (pga_low / 0.1);
elseif (pga4nl <= a2)
    Fnl = bnl * log (pga_low / 0.1) + c * (log (pga4nl / a1))^2 + d * (log (pga4nl / a1))^3;
else
    Fnl = bnl * log (pga4nl / 0.1);
end
Fs = Flin + Fnl;

% Compute median and sigma
lny = Fm + Fd + Fs;
median = exp (lny); 
sigma = (U==1) * sigtu(ip) + (U~=1) * sigtm(ip);


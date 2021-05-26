function [median, sigma, period1] = cb_2008_active(T,rup,site,A1100)

% Created by Yoshifumi Yamamoto,  5/8/10, yama4423@stanford.edu
% Updated by Emily Mongold, 11/27/20
%
% Source Model:
% Campbell, K.W., Bozorgnia, Y. (2008). "NGA Ground Motino Model for the 
% Geometric Mean Horizontal Component of PGA, PGV, PGD, and 5% Damped 
% Linear Elastic Response Spectra for Periods Ranging from 0.01 to 10s."
% Earthquake Spectra, 24(1), 139-171.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = Period (sec); Use -1 for PGV computation and -10 for PGD
%                   computation
%                   Use 1000 for output the array of Sa with period
%   rup           = rupture object input containing the following
%                   variables:
%       M             = Magnitude
%       Rrup          = Closest distance coseismic rupture (km)
%       Rjb           = Joyner-Boore distance (km)
%       Ztor          = Depth to the top of coseismic rupture (km)
%       delta         = average dip of the rupture place (degree)
%       lambda        = rake angle (degree)
%   site          = site object input containing the following
%                   variables:
%       Vs30          = shear wave velocity averaged over top 30 m (m/s)
%       Z25           = Depth to the 2.5 km/s shear-wave velocity horizon (km)
%   A1100         = median estimate of PGA on reference rock with Vs30 =
%                   1100 m/s
% OUTPUT
%   median        = Median spectral acceleration prediction
%   sigma         = logarithmic standard deviation of spectral acceleration
%                   prediction
%   period1       = periods for which the median and sigma values are
%                   provided. If T = 1000, then period1 = the full set of
%                   available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if (nargin < 4) % no Hanging wall indicator supplied
    A1100 = 0;
end

% Coefficients
period = [0.01	0.02	0.03	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75	1	1.5	2	3	4	5	7.5	10 0 -1 -10];
        % Style of faulting
if rup.Ztor < 1
    ffltz = rup.Ztor;
else    
    ffltz = 1;
end
if rup.Rjb == 0
    fhngr = 1;
else    
    if rup.Ztor < 1
        fhngr = (max(rup.Rrup,sqrt(rup.Rjb^2+1))-rup.Rjb)/max(rup.Rrup,sqrt(rup.Rjb^2+1));
    else    
        fhngr = (rup.Rrup - rup.Rjb)/rup.Rrup;
    end
end

Frv = (rup.lambda > 30 & rup.lambda < 150);
Fnm = (rup.lambda > -150 & rup.lambda < -30);
if rup.M <= 6
    fhngm = 0;
elseif rup.M < 6.5
    fhngm = 2 * (rup.M - 6);
else    
    fhngm = 1;
end

nT=length(T);
iflg=0;
if(nT==1)
    if(T==1000)
        iflg=1;
        %     Computation with original periods
        nperi=length(period)-3;
        median=zeros(1,nperi);
        sigma=zeros(1,nperi);
        period1=period(1:nperi);
        for i=1:1:nperi
            [median(i), sigma(i)] = CB_2008_nga_sub(rup,i,site,ffltz,fhngr,Frv,Fnm,fhngm,A1100);
        end
    end
end
if(iflg==0)
    median=zeros(1,nT);
    sigma=zeros(1,nT);
    period1=T;
    for it=1:1:nT
        Teach=T(it);
        % interpolate between periods if neccesary    
        if (isempty(find(period == Teach, 1)))
            T_low = max(period(period<Teach));
            T_hi = min(period(period>Teach));

            [sa_low, sigma_low] = cb_2008_active (T_low,rup,site,A1100);
            [sa_hi, sigma_hi] = cb_2008_active (T_hi,rup,site,A1100);

            x = [log(T_low) log(T_hi)];
            Y_sa = [log(sa_low) log(sa_hi)];
            Y_sigma = [sigma_low sigma_hi];
            median(it) = exp(interp1(x,Y_sa,log(Teach)));
            sigma(it) = interp1(x,Y_sigma,log(Teach));
        else
            i = find(period == Teach); 
            [median(it), sigma(it)] = CB_2008_nga_sub (rup,i,site,ffltz,fhngr,Frv,Fnm,fhngm,A1100);
        end
    end
end

function [median, sigma] = CB_2008_nga_sub (rup,ip,sitevar,ffltz,fhngr,Frv,Fnm,fhngm,A1100)

c0 = [-1.715	-1.68	-1.552	-1.209	-0.657	-0.314	-0.133	-0.486	-0.89	-1.171	-1.466	-2.569	-4.844	-6.406	-8.692	-9.701	-10.556	-11.212	-11.684	-12.505	-13.087	-1.715	0.954	-5.27];
c1 = [0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.5	0.656	0.972	1.196	1.513	1.6	1.6	1.6	1.6	1.6	1.6	0.5	0.696	1.6];
c2 = [-0.53	-0.53	-0.53	-0.53	-0.53	-0.53	-0.53	-0.446	-0.362	-0.294	-0.186	-0.304	-0.578	-0.772	-1.046	-0.978	-0.638	-0.316	-0.07	-0.07	-0.07	-0.53	-0.309	-0.07];
c3 = [-0.262	-0.262	-0.262	-0.267	-0.302	-0.324	-0.339	-0.398	-0.458	-0.511	-0.592	-0.536	-0.406	-0.314	-0.185	-0.236	-0.491	-0.77	-0.986	-0.656	-0.422	-0.262	-0.019	0];
c4 = [-2.118	-2.123	-2.145	-2.199	-2.277	-2.318	-2.309	-2.22	-2.146	-2.095	-2.066	-2.041	-2	-2	-2	-2	-2	-2	-2	-2	-2	-2.118	-2.016	-2];
c5 = [0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17	0.17];
c6 = [5.6	5.6	5.6	5.74	7.09	8.05	8.79	7.6	6.58	6.04	5.3	4.73	4	4	4	4	4	4	4	4	4	5.6	4	4];
c7 = [0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.28	0.255	0.161	0.094	0	0	0	0	0	0.28	0.245	0];
c8 = [-0.12	-0.12	-0.12	-0.12	-0.12	-0.099	-0.048	-0.012	0	0	0	0	0	0	0	0	0	0	0	0	0	-0.12	0	0];
c9 = [0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.49	0.371	0.154	0	0	0	0	0.49	0.358	0];
c10 = [1.058	1.102	1.174	1.272	1.438	1.604	1.928	2.194	2.351	2.46	2.587	2.544	2.133	1.571	0.406	-0.456	-0.82	-0.82	-0.82	-0.82	-0.82	1.058	1.694	-0.82];
c11 = [0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.04	0.077	0.15	0.253	0.3	0.3	0.3	0.3	0.3	0.3	0.04	0.092	0.3];
c12 = [0.61	0.61	0.61	0.61	0.61	0.61	0.61	0.61	0.70	0.75	0.85	0.883	1	1	1	1	1	1	1	1	1	0.61	1	1];
k1 = [865	865	908	1054	1086	1032	878	748	654	587	503	457	410	400	400	400	400	400	400	400	400	865	400	400];
k2 = [-1.186	-1.219	-1.273	-1.346	-1.471	-1.624	-1.931	-2.188	-2.381	-2.518	-2.657	-2.669	-2.401	-1.955	-1.025	-0.299	0	0	0	0	0	-1.186	-1.955	0];
k3 = [1.839	1.84	1.841	1.843	1.845	1.847	1.852	1.856	1.861	1.865	1.874	1.883	1.906	1.929	1.974	2.019	2.11	2.2	2.291	2.517	2.744	1.839	1.929	2.744];
slny = [0.478   0.480   0.489   0.510   0.520   0.531   0.532   0.534   0.534   0.544   0.541   0.550   0.568   0.568   0.564   0.571   0.558   0.576   0.601   0.628   0.667   0.478   0.484   0.667];
tlny = [0.219   0.219   0.235   0.258   0.292   0.286   0.280   0.249   0.240   0.215   0.217   0.214   0.227   0.255   0.296   0.296   0.326   0.297   0.359   0.428   0.485   0.219   0.203   0.485];
sigmac = [0.166	0.166	0.165	0.162	0.158	0.17	0.18	0.186	0.191	0.198	0.206	0.208	0.221	0.225	0.222	0.226	0.229	0.237	0.237	0.271	0.29	0.166	0.19	0.29];
sigmat   = [0.526 0.528 0.543 0.572 0.596 0.603 0.601 0.589 0.585 0.585 0.583 0.590 0.612 0.623 0.637 0.643 0.646 0.648 0.700 0.760 0.825 0.526 0.525 0.825];
sigmaarb = [0.551 0.553 0.567 0.594 0.617 0.627 0.628 0.618 0.616 0.618 0.618 0.626 0.650 0.662 0.675 0.682 0.686 0.690 0.739 0.807 0.874 0.551 0.558 0.874];
roh  = [1.000   0.999   0.989   0.963   0.922   0.898   0.890   0.871   0.852   0.831   0.785   0.735   0.628   0.534   0.411   0.331   0.289   0.261   0.200   0.174   0.174   1.000   0.691   0.174];

c = 1.88;
n = 1.18;

% Magnitude dependence
if rup.M <= 5.5
    fmag = c0(ip) + c1(ip) * rup.M;
elseif rup.M<=6.5
    fmag = c0(ip) + c1(ip) * rup.M + c2(ip) * (rup.M - 5.5); 
else    
    fmag = c0(ip) + c1(ip) * rup.M + c2(ip) * (rup.M - 5.5) + c3(ip) * (rup.M - 6.5); 
end

% Distance dependence
fdis = (c4(ip) + c5(ip) * rup.M) * log(sqrt(rup.Rrup^2 + c6(ip)^2));

% Style of faulting
fflt = c7(ip) * Frv * ffltz + c8(ip) * Fnm;

% Hanging-wall effects
fhngz = ((20 - rup.Ztor)/20) * (rup.Ztor >= 0 && rup.Ztor < 20);

fhngdelta = (rup.delta <= 70) + ((90 - rup.delta)/20) * (rup.delta > 70);

fhng = c9(ip) * fhngr * fhngm * fhngz * fhngdelta;

% Site conditions
site1100 = site(sitevar.is_soil,1100,0,sitevar.Z25,sitevar.Z10,sitevar.Zbot,sitevar.region);
if sitevar.Vs30 < k1(ip)
    if A1100==0
        A1100 = cb_2008_active (0,rup,site1100); 
    end
    fsite = c10(ip) * log(sitevar.Vs30/k1(ip)) + k2(ip) * (log(A1100 + c * (sitevar.Vs30/k1(ip))^n) - log(A1100 + c));
elseif sitevar.Vs30 < 1100
    fsite = (c10(ip) + k2(ip) * n) * log(sitevar.Vs30/k1(ip));
else    
    fsite = (c10(ip) + k2(ip) * n) * log(1100/k1(ip));
end

% Sediment effects
if sitevar.Z25 < 1
    fsed = c11(ip) * (sitevar.Z25 - 1);
elseif sitevar.Z25 <= 3
    fsed = 0;
else    
    fsed = c12(ip) * k3(ip) * exp(-0.75) * (1 - exp(-0.25 * (sitevar.Z25 - 3)));
end

% Median value
median = exp(fmag + fdis + fflt + fhng + fsite + fsed);

% Standard deviation computations
if (sitevar.Vs30 < k1(ip))
    alpha = k2(ip) * A1100 * ((A1100 + c*(sitevar.Vs30/k1(ip))^n)^(-1) - (A1100+c)^(-1));
else
    alpha = 0;
end

slnaf=0.3;
slnpga=slny(length(slny)-2);
slnab=sqrt(slnpga^2  -slnaf^2);
slnyb=sqrt(slny(ip)^2  -slnaf^2);
% from equation in the paper
sigmat = sqrt(slny(ip)^2 + alpha^2*slnab^2 + 2*alpha*roh(ip)*slnyb*slnab);
% from graph in the paper
% sigmat = sqrt(slny(ip)^2 + alpha^2*slnab^2 + 2*alpha*roh(ip)*slny(ip)*slnab);
tau = tlny(ip);

sigma=sqrt(sigmat^2 + tau^2);


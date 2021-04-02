function [median, sigma, period1] = cy_2008_active(T,rup,site)

% Created by Yoshifumi Yamamoto, 11/10/08, yama4423@stanford.edu
%   updated 2009/05/05
%   updated 2010/05/20  fixing bug about nonlinear effect
%   updated 2012/07/17  by Jack Baker to complete documentation of inputs and improve interpolation
% Updated by Emily Mongold, 11/27/20
%
% Source Model:
% Chiou, B. S.-J. and Youngs R. R.(2008). "An NGA Model for the Average 
% Horizontal Component of Peak Ground Motion and Response Spectra."  
% Earthquake Spectra, 24(1), 173-215.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T           = Period (0.01 - 10s)
%   rup         = rupture object input containing the following
%                 variables:
%       M           = Magnitude 
%       Rrup        = Closest distance to the rupture plane (km)
%       Rjb         = Joyner-Boore distance to the rupture plane (km)
%       Rx          = Distance from surface projection of updip edge of fault
%                   rupture, perpendicular to strike (km)
%       Ztor        = Depth to top of rupture (km)
%       AS          = 1 if Aftershock
%       delta       = Dip
%       lambda      = Rake angle
%   site        = site object input containing the following
%                 variables:
%       Vs30        = Average shear velocity for top 30 m (m/s)
%       fvs30       = 1 for Vs30 is inferred from geology
%                   = 0 for measured  Vs30
%       Z10         = Depth to 1 km/s shear wave velocity
% OUTPUT
%   median      = median spectral acceleration prediction
%   sigma       = Intra-event standard error + Inter-event standard error
%   period1     = periods for which the median and sigma values are
%                 provided. If T = 1000, then period1 = the full set of
%                 available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Coefficients
period = [0 -1 0.01	0.02	0.03	0.04	0.05	0.075   0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.75 1.0 1.5 2.0 3.0 4.0 5.0 7.5 10.0];

deltar=rup.delta*pi()/180.0;
frv = rup.lambda >= 30 & rup.lambda <= 150; % frv: 1 for lambda between 30 and 150, 0 otherwise
fnm = rup.lambda >= -120 & rup.lambda <= -60; % fnm: 1 for lambda between -120 and -60, 0 otherwise
HW = rup.Rx>=0;
% rup.AS = rup.AS==1;
nT=length(T);
iflg=0;
if(nT==1)
    if(T==1000)
        iflg=1;
	    nperi=length(period);
	    median=zeros(1,nperi-2);
	    sigma=zeros(1,nperi-2);
	    tau=zeros(1,nperi-2);
	    period1=period(3:nperi);
	    for i=3:1:nperi
	        [median(i-2),sigma(i-2),tau(i-2)] = CY_2008_nga_sub(rup.M,i,rup.Rrup,rup.Rjb,rup.Ztor,site.Vs30,deltar, frv, fnm,rup.AS,rup.Rx,rup.HW,site.Z10,site.fvs30);
        end
    end
end
if(iflg==0)
    median=zeros(1,nT);
    sigma=zeros(1,nT);
    tau=zeros(1,nT);
    period1=T;
    for it=1:1:nT
        Teach=T(it);
        % interpolate between periods if neccesary    
        if (~any(period == Teach))
            T_low = max(period(find(period<Teach)));
            T_hi =  min(period(find(period>Teach)));

            [sa_low, sigma_low, tau_low] = cy_2008_active(T_low,rup,site);
            [sa_hi,  sigma_hi,  tau_hi]  = cy_2008_active(T_hi,rup,site);
            
            x = [log(T_low) log(T_hi)];
            Y_sa = [log(sa_low) log(sa_hi)];
            Y_sigma = [sigma_low sigma_hi];
            %Y_tau = [tau_low tau_hi];
            median(it) = exp(interp1(x,Y_sa,log(Teach)));
            sigma(it) = interp1(x,Y_sigma,log(Teach));
            %tau(it) = interp1(x,Y_tau,log(Teach));
        else
            i = find(period == Teach); % Identify the period
            [median(it), sigma(it), tau(it)] = CY_2008_nga_sub(rup.M,i,rup.Rrup,rup.Rjb,rup.Ztor,site.Vs30,deltar,frv,fnm,rup.AS,rup.Rx,rup.HW,site.Z10,site.fvs30);
        end
    end
end

function [median, sigma, tau] = CY_2008_nga_sub(M, ip, Rrup, Rjb, Ztor, Vs30, deltar, frv, fnm, AS, Rx, HW, Z10, FVS30)
c1     = [-1.2687   2.2884  -1.2687 -1.2515	-1.1744	-1.0671	-0.9464	-0.7051	-0.5747	-0.5309	-0.6352	-0.7766	-0.9278	-1.2176	-1.4695	-1.9278	-2.2453	-2.7303	-3.1413	-3.7413	-4.1814	-4.5178	-5.1224	-5.5872];
c1a    = [0.1       0.1094  0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.1     0.0999  0.0997  0.0991  0.0936  0.0766  0.0022  -0.0591 -0.0931 -0.0982 -0.0994 -0.0999 -0.1000];
c1b    = [-0.2550   -0.0626 -0.2550 -0.2550 -0.2550 -0.2550 -0.2550 -0.2540 -0.2530 -0.2500 -0.2449 -0.2382 -0.2313 -0.2146 -0.1972 -0.1620 -0.1400 -0.1184 -0.1100 -0.1040 -0.1020 -0.1010 -0.1010 -0.1000];
c2     = 1.06;
c3     = 3.45;
cn     = [2.996     1.648   2.996   3.292   3.514   3.563   3.547   3.448   3.312   3.044   2.831   2.658   2.505   2.261   2.087   1.812   1.648   1.511   1.470   1.456   1.465   1.478   1.498   1.502];
cm     = [4.1840    4.2979  4.1840  4.1879  4.1556  4.1226  4.1011  4.0860  4.1030  4.1717  4.2476  4.3184  4.3844  4.4979  4.5881  4.7571  4.8820  5.0697  5.2173  5.4385  5.5977  5.7276  5.9891  6.1930];
c4     = -2.1;
c4a    = -0.5;
crb    = 50.0;
c5     = [6.1600    5.1760  6.1600  6.1580  6.1550  6.1508  6.1441  6.1200  6.0850  5.9871  5.8699  5.7547  5.6527  5.4997  5.4029  5.2900  5.2480  5.2194  5.2099  5.2040  5.2020  5.2010  5.2000  5.2000];
c6     = [0.4893    0.4407  0.4893  0.4892  0.4890  0.4888  0.4884  0.4872  0.4854  0.4808  0.4755  0.4706  0.4665  0.4607  0.4571  0.4531  0.4517  0.4507  0.4504  0.4501  0.4501  0.4500  0.4500  0.4500];
chm    = 3.0;
c7     = [0.0512    0.0207  0.0512  0.0512  0.0511  0.0508  0.0504  0.0495  0.0489  0.0479  0.0471  0.0464  0.0458  0.0445  0.0429  0.0387  0.0350  0.0280  0.0213  0.0106  0.0041  0.0010  0.0000  0.0000];
c7a    = [0.0860    0.0437  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0860  0.0850  0.0830  0.0690  0.0450  0.0134  0.0040  0.0010  0.0000  0.0000  0.0000  0.0000];
c9     = [0.7900    0.3079  0.7900  0.8129  0.8439  0.8740  0.8996  0.9442  0.9677  0.9660  0.9334  0.8946  0.8590  0.8019  0.7578  0.6788  0.6196  0.5101  0.3917  0.1244  0.0086  0.0000  0.0000  0.0000];
c9a    = [1.5005    2.6690  1.5005  1.5028  1.5071  1.5138  1.5230  1.5597  1.6104  1.7549  1.9157  2.0709  2.2005  2.3886  2.5000  2.6224  2.6690  2.6985  2.7085  2.7145  2.7164  2.7172  2.7177  2.7180];
c10    = [-0.3218   -0.1166 -0.3218 -0.3323 -0.3394 -0.3453 -0.3502 -0.3579 -0.3604 -0.3565 -0.3470 -0.3379 -0.3314 -0.3256 -0.3189 -0.2702 -0.2059 -0.0852 0.0160  0.1876  0.3378  0.4579  0.7514  1.1856];
cy1    = [-0.00804  -0.00275    -0.00804    -0.00811    -0.00839    -0.00875    -0.00912    -0.00973    -0.00975    -0.00883    -0.00778    -0.00688    -0.00612    -0.00498    -0.00420    -0.00308    -0.00246    -0.00180    -0.00147    -0.00117    -0.00107    -0.00102    -0.00096    -0.00094];
cy2    = [-0.00785  -0.00625    -0.00785    -0.00792    -0.00819    -0.00855    -0.00891    -0.00950    -0.00952    -0.00862    -0.00759    -0.00671    -0.00598    -0.00486    -0.00410    -0.00301    -0.00241    -0.00176    -0.00143    -0.00115    -0.00104    -0.00099    -0.00094    -0.00091];
cy3    = 4.0;
phi1   = [-0.4417   -0.7861 -0.4417 -0.4340 -0.4177 -0.4000 -0.3903 -0.4040 -0.4423 -0.5162 -0.5697 -0.6109 -0.6444 -0.6931 -0.7246 -0.7708 -0.7990 -0.8382 -0.8663 -0.9032 -0.9231 -0.9222 -0.8346 -0.7332];
phi2   = [-0.1417   -0.0699 -0.1417 -0.1364 -0.1403 -0.1591 -0.1862 -0.2538 -0.2943 -0.3113 -0.2927 -0.2662 -0.2405 -0.1975 -0.1633 -0.1028 -0.0699 -0.0425 -0.0302 -0.0129 -0.0016 0.0000	0.000   0.000];
phi3   = [-0.007010 -0.008444   -0.007010   -0.007279   -0.007354   -0.006977   -0.006467   -0.005734   -0.005604   -0.005845   -0.006141   -0.006439   -0.006704   -0.007125   -0.007435   -0.008120   -0.008444   -0.007707   -0.004792   -0.001828   -0.001523   -0.01440    -0.001369   -0.001361];
phi4   = [0.102151	5.41000     0.102151    0.108360    0.119888    0.133641    0.148927    0.190596    0.230662    0.266468    0.255253    0.231541    0.207277    0.165464    0.133828    0.085153    0.058595    0.031787    0.019716    0.009643    0.005379    0.003223    0.001134    0.000515];
phi5   = [0.2289    0.2899  0.2289  0.2289  0.2289  0.2289	0.2290  0.2292  0.2297  0.2326  0.2386  0.2497  0.2674  0.3120  0.3610  0.4353  0.4629  0.4756  0.4785  0.4796  0.4799  0.4799  0.4800  0.4800];
phi6   = [0.014996  0.006718    0.014996    0.014996    0.014996    0.014996    0.014996    0.014996    0.014996    0.014988    0.014964    0.014881    0.014639    0.013493    0.011133	0.006739    0.005749    0.005544 	0.005521    0.005517    0.005517    0.005517    0.005517    0.005517];
phi7   = [580.0     459.0   580.0   580.0   580.0   579.9   579.9   579.6   579.2   577.2   573.9   568.5   560.5   540.0   512.9   441.9   391.8   348.1   332.5   324.1   321.7   320.9   320.3   320.1];
phi8   = [0.0700    0.1138  0.0700  0.0699  0.0701  0.0702  0.0701  0.0686  0.0646  0.0494  -0.0019 -0.0479 -0.0756 -0.0960 -0.0998 -0.0765 -0.0412 0.0140  0.0544  0.1232  0.1859  0.2295  0.2660  0.2682];
tau1   = [0.3437    0.2539  0.3437  0.3471  0.3603  0.3718  0.3848  0.3878  0.3835  0.3719  0.3601  0.3522  0.3438  0.3351  0.3353  0.3429  0.3577  0.3769  0.4023  0.4406  0.4784  0.5074  0.5328  0.5542];
tau2   = [0.2637    0.2381  0.2637  0.2671  0.2803  0.2918  0.3048  0.3129  0.3152  0.3128  0.3076  0.3047  0.3005  0.2984  0.3036  0.3205  0.3419  0.3703  0.4023  0.4406  0.4784  0.5074  0.5328  0.5542];
sigma1 = [0.4458    0.4496  0.4458  0.4458  0.4535  0.4589  0.4630  0.4702  0.4747  0.4798  0.4816  0.4815  0.4801  0.4758  0.4710  0.4621  0.4581  0.4493  0.4459  0.4433  0.4424  0.4420  0.4416  0.4414];
sigma2 = [0.3459    0.3554  0.3459  0.3459  0.3537  0.3592  0.3635  0.3713  0.3769  0.3847  0.3902  0.3946  0.3981  0.4036  0.4079  0.4157  0.4213  0.4213  0.4213  0.4213  0.4213  0.4213  0.4213  0.4213];
sigma3 = [0.8000    0.7504  0.8000  0.8000  0.8000  0.8000  0.8000  0.8000  0.8000  0.8000  0.8000  0.7999  0.7997  0.7988  0.7966  0.7792  0.7504  0.7136  0.7035  0.7006  0.7001  0.7000  0.7000  0.7000];
sigma4 = [0.0663    0.0133  0.0663  0.0663  0.0663  0.0663  0.0663  0.0663  0.0663  0.0612  0.0530  0.0457  0.0398  0.0312  0.0255  0.0175  0.0133  0.0090  0.0068  0.0045  0.0034  0.0027  0.0018  0.0014];
term1 = c1(ip);   
term2 = (c1a(ip)*frv + c1b(ip)*fnm + c7(ip)*(Ztor-4))*(1-AS) + (c10(ip) + c7a(ip)*(Ztor-4))*AS;
term5 = c2*(M - 6);    
term6 = ((c2 - c3)/cn(ip)) * log (1 + exp(cn(ip) * (cm(ip) - M)));    
term7 = c4 * log(Rrup + c5(ip) * cosh(c6(ip) * max(M - chm,0)));    
term8 = (c4a - c4) * log (sqrt(Rrup^2 + crb^2));    
term9 = (cy1(ip) + cy2(ip)/cosh(max(M-cy3,0)))*Rrup;    
term10 = c9(ip)*HW * tanh(Rx*cos(deltar)^2/c9a(ip))*(1-sqrt(Rjb^2+Ztor^2)/(Rrup + 0.001));

Sa1130 = exp(term1 + term2 + term5 + term6 + term7 + term8 + term9 + term10);

term11 = phi1(ip) * min(log(Vs30/1130),0);    
term12 = phi2(ip) * (exp(phi3(ip) * (min(Vs30,1130) - 360)) - exp(phi3(ip) * (1130 - 360))) * log((Sa1130 + phi4(ip))/phi4(ip));
term13 = phi5(ip) * (1-1/cosh(phi6(ip)*max(0,Z10-phi7(ip)))) + phi8(ip)/cosh(0.15*max(0,Z10-15));

% Compute median
median = exp(log(Sa1130) + term11 + term12 + term13);

% Compute standard deviation
Finferred=(FVS30==1); % 1: Vs30 is inferred from geology.
Fmeasured=(FVS30==0); % 1: Vs30 is measured.
b=phi2(ip)*( exp(phi3(ip)*(min(Vs30,1130)-360))-exp(phi3(ip)*(1130-360)) );
c=phi4(ip);
NL0=b*Sa1130/(Sa1130+c);
sigma = (sigma1(ip)+(sigma2(ip) - sigma1(ip))/2*(min(max(M,5),7)-5) + sigma4(ip)*AS)*sqrt((sigma3(ip)*Finferred + 0.7* Fmeasured) + (1+NL0)^2);
tau = tau1(ip) + (tau2(ip)-tau1(ip))/2 * (min(max(M,5),7)-5);
sigma=sqrt((1+NL0)^2*tau^2+sigma^2);


function [median, sigma, period1] = sp_2016_stable(T,rup)

% Created by Abhineet Gupta, 2016-07-03, Stanford University
% lightly edited by Jack Baker to add documentation, 12/4/2017
% Updated by Emily Mongold, 11/27/20
%
% Source Model:
% Shahjouei, A., and Pezeshk, S. (2016). "Alternative Hybrid Empirical 
% Ground-Motion Model for Central and Eastern North America Using Hybrid 
% Simulations and NGA-West2 Models." Bulletin of the Seismological Society 
% of America, 106(2), 734-754.
%
% Applicable for 5 <= M <= 8, 2 <= Rjb <= 1000 and Vs30 = 3000 m/s
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = period of vibration, use 0 for PGA, -1 for PGV
%   rup           = rupture object input containing the following
%                   variables:
%       M               = moment magnitude
%       Rjb             = Joyner-Boore distance, closest horizontal distance 
%                         to the vertical projection of the rupture plane (km)
% OUTPUT   
%   median        = median spectral acceleration, units of g, except for PGV in cm/s
%   sigma         = natural logarithm standard deviation of spectral acceleration
%   period1       = periods for which the median and sigma values are
%                   provided. If T = 1000, then period1 = the full set of
%                   available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% for the given period T, get the index for the constants
period = [-1           0           0.010       0.020       0.030       0.040       0.050       0.075       0.100       0.150       0.200       0.250       0.300       0.400       0.500       0.750       1.000       1.500       2.000       3.000       4.000       5.000       7.500       10.000     ];
c1     = [-2.3891       -0.3002     -0.3472     0.8320      1.1850      1.2460      1.1793      0.8045      0.3500      -0.5264     -1.2884     -1.9422     -2.5071     -3.4360     -4.1699     -5.4797     -6.3464     -7.4087     -8.0057     -8.5793     -8.8246     -8.9855     -9.3927     -9.7350    ];
c2     = [1.259E00      5.066E-01   4.838E-01   1.934E-01   1.064E-01   8.986E-02   1.037E-01   1.866E-01   2.871E-01   4.782E-01   6.413E-01   7.789E-01   8.961E-01   1.085E00    1.231E00    1.482E00    1.641E00    1.823E00    1.916E00    1.985E00    1.990E00    1.975E00    1.925E00    1.879E00   ];
c3     = [-7.901E-02    -4.526E-02  -4.093E-02  -2.060E-02  -1.423E-02  -1.268E-02  -1.321E-02  -1.788E-02  -2.381E-02  -3.519E-02  -4.486E-02  -5.295E-02  -5.976E-02  -7.059E-02  -7.878E-02  -9.245E-02  -1.006E-01  -1.093E-01  -1.130E-01  -1.146E-01  -1.131E-01  -1.105E-01  -1.032E-01  -9.666E-02 ];
c4     = [-2.9386       -3.2240     -3.0832     -3.1134     -3.1029     -3.0785     -3.0488     -2.9697     -2.8940     -2.7610     -2.6504     -2.5573     -2.4780     -2.3495     -2.2510     -2.0865     -1.9931     -1.9162     -1.9173     -2.0184     -2.1475     -2.2496     -2.3572     -2.4139    ];
c5     = [3.034E-01     2.998E-01   2.712E-01   2.786E-01   2.792E-01   2.773E-01   2.744E-01   2.660E-01   2.576E-01   2.426E-01   2.301E-01   2.196E-01   2.107E-01   1.961E-01   1.849E-01   1.659E-01   1.546E-01   1.438E-01   1.418E-01   1.499E-01   1.635E-01   1.764E-01   1.973E-01   2.117E-01  ];
c6     = [-9.290E-03    -1.283E00   -9.676E-01  -1.133E00   -1.078E00   -9.743E-01  -8.635E-01  -6.122E-01  -4.123E-01  -1.319E-01  4.637E-02   1.631E-01   2.407E-01   3.244E-01   3.544E-01   3.284E-01   2.530E-01   9.019E-02   -3.828E-02  -1.744E-01  -1.844E-01  -1.043E-01  3.465E-01   1.010E00   ];
c7     = [-4.605E-02    1.045E-01   4.983E-02   5.994E-02   5.239E-02   4.160E-02   3.077E-02   7.491E-03   -1.012E-02  -3.338E-02  -4.690E-02  -5.478E-02  -5.919E-02  -6.197E-02  -6.046E-02  -4.979E-02  -3.709E-02  -1.551E-02  -1.252E-03  9.393E-03   3.919E-03   -1.187E-02  -7.832E-02  -1.678E-01 ];
c8     = [-2.7548       -3.0856     -2.9695     -3.5023     -3.5722     -3.5083     -3.3986     -3.0852     -2.7947     -2.3312     -1.9927     -1.7399     -1.5470     -1.2793     -1.1111     -0.9131     -0.8641     -0.9200     -1.0327     -1.2453     -1.3849     -1.4511     -1.3728     -1.0631    ];
c9     = [3.467E-01     2.778E-01   2.693E-01   2.901E-01   2.865E-01   2.769E-01   2.659E-01   2.391E-01   2.163E-01   1.818E-01   1.576E-01   1.398E-01   1.265E-01   1.085E-01   9.757E-02   8.570E-02   8.405E-02   9.103E-02   1.016E-01   1.214E-01   1.357E-01   1.446E-01   1.490E-01   1.370E-01  ];
c10    = [-7.623E-04    -7.711E-04  -6.695E-04  -5.857E-04  -6.220E-04  -6.818E-04  -7.439E-04  -8.801E-04  -9.848E-04  -1.125E-03  -1.209E-03  -1.258E-03  -1.286E-03  -1.304E-03  -1.294E-03  -1.219E-03  -1.123E-03  -9.407E-04  -7.926E-04  -5.919E-04  -4.855E-04  -4.439E-04  -5.176E-04  -7.420E-04 ];
c11    = [-4.598E00     3.810E00    -4.434E00   -4.412E00   -4.353E00   -4.303E00   -4.266E00   -4.214E00   4.201E00    4.239E00    4.325E00    4.438E00    4.571E00    -4.872E00   -5.211E00   -6.154E00   -7.174E00   -9.253E00   -1.122E01   1.438E01    1.619E01    1.671E01    1.458E01    1.123E01   ];
c12    = [-5.54E-02     -4.10E-02   -5.60E-02   -5.59E-02   -5.77E-02   -5.77E-02   -5.78E-02   -5.61E-02   -5.65E-02   -5.59E-02   -5.60E-02   -5.37E-02   -5.11E-02   -4.70E-02   -4.42E-02   -3.84E-02   -3.14E-02   -2.27E-02   -1.84E-02   -1.89E-02   -1.60E-02   -1.53E-02   -1.43E-02   -1.70E-02  ];
c13    = [9.78E-01      8.76E-01    9.82E-01    9.83E-01    1.00E00     1.01E00     1.03E00     1.03E00     1.05E00     1.04E00     1.03E00     1.02E00     1.01E00     9.87E-01    9.81E-01    9.67E-01    9.33E-01    8.83E-01    8.57E-01    8.59E-01    8.30E-01    8.26E-01    8.15E-01    8.22E-01   ];
c14    = [6.63E-01      6.11E-01    6.64E-01    6.65E-01    6.76E-01    6.88E-01    7.01E-01    7.21E-01    7.32E-01    7.24E-01    7.15E-01    7.12E-01    7.18E-01    7.25E-01    7.36E-01    7.60E-01    7.70E-01    7.76E-01    7.78E-01    7.77E-01    7.66E-01    7.66E-01    7.62E-01    7.52E-01   ];
sigReg = [1.00E-01      1.94E-01    1.32E-01    9.28E-02    8.33E-02    7.98E-02    7.76E-02    7.38E-02    7.17E-02    7.16E-02    7.43E-02    7.79E-02    8.15E-02    8.76E-02    9.23E-02    9.91E-02    1.02E-01    1.05E-01    1.06E-01    1.07E-01    1.07E-01    1.07E-01    1.13E-01    1.40E-01   ];
sigPar = [2.88E-01      3.73E-01    2.81E-01    2.81E-01    2.77E-01    2.79E-01    2.72E-01    2.52E-01    2.65E-01    2.76E-01    2.58E-01    2.68E-01    2.84E-01    3.40E-01    3.57E-01    3.74E-01    3.92E-01    4.26E-01    4.40E-01    5.80E-01    5.89E-01    6.31E-01    7.21E-01    7.39E-01   ];

if length (T) == 1 && T == 1000
    median=zeros(1,length(period)-2);
    sigma=zeros(1,length(period)-2);
    period1=period(3:end);
else
    median=zeros(1,length(T));
    sigma=zeros(1,length(T));
    period1=T;
end

for i = 1:length(T)
    Ti = T(i);
    if Ti < 0 && Ti ~= -1
        error('T must be non-negative or -1');
    end
    indx = find(abs(period - Ti) < 1e-4);
    if isempty(indx)
        % The user defined period requires interpolation
        indx_high = find(period > Ti, 1, 'first');
        if isempty(indx_high)
            error(['T exceeds the maximum period supplied for this model (' num2str(period(end)) ' s)']);
        end
        T_high = period(indx_high);
        indx_low = indx_high - 1;
        T_low = period(indx_low);
        [Sa_high, sigma_high] = sp_2016_stable(T_high,rup);
        [Sa_low, sigma_low] = sp_2016_stable(T_low,rup);

        % interpolate
        x = [log(T_low) log(T_high)];
        Y_sa = [log(Sa_low) log(Sa_high)];
        Y_sigma = [sigma_low sigma_high];
        median(i) = exp(interp1(x, Y_sa, log(Ti)));
        sigma(i) = interp1(x, Y_sigma, log(Ti));
    else
        R = sqrt(rup.Rjb^2 + c11(indx)^2);
        log10Sa = c1(indx) + c2(indx)*rup.M + c3(indx)*(rup.M^2) + ...
            (c4(indx) + c5(indx)*rup.M)*min([log10(R), log10(60)]) + ...
            (c6(indx) + c7(indx)*rup.M)*max([min([log10(R/60), log10(120/60)]), 0]) + ...
            (c8(indx) + c9(indx)*rup.M)*max([log10(R/120), 0]) + c10(indx)*R;
        median(i) = 10^(log10Sa);

        % Aleatory uncertainty
        if Ti == -1
            psi = -3.054E-05;
        else
            psi = -6.898E-03;
        end
        if rup.M <= 6.5
            sigma_1 = c12(indx)*rup.M + c13(indx);
        else
            sigma_1 = psi*rup.M + c14(indx);
        end
        sigma_Al = sqrt(sigma_1^2 + sigReg(indx)^2);

        % Epistemic uncertainty
        if Ti < 1
            if rup.M < 7
                sigma_2 = 0.072;
            else
                sigma_2 = 0.0665*(rup.M - 7) + 0.072;
            end
        else
            if rup.M < 7
                sigma_2 = 0.072 + 0.0217*log(Ti);
            else
                sigma_2 = 0.0665*(rup.M - 7) + 0.072 + 0.0217*log(Ti);
            end
        end
        sigma_Ep = sqrt(sigma_2^2 + sigPar(indx)^2);
        % disregard sigma_Ep since authors state - The epistemic uncertainty for an individual GMM is infrequently employed...
        sigma_Ep = zeros(size(sigma_Ep));
        
        % Total uncertainty
        sigma(i) = sqrt(sigma_Al^2 + sigma_Ep^2);
    end
end
function [median,sigma, period1] = a_2015_active(T,rup)

% Created by Jack Baker, 2/26/2015, bakerjw@stanford.edu
% Updated by Emily Mongold, 11/25/2020
%
% updated 3/3/2015 to correct an error in the original electronically
% posted version of the model
%
% updated 4/14/2015 to change default h_eff from Yenier and Atkinson (2014)
% to alternate h_eff proposed in Atkinson (2015) based on personal 
% communication with Atkinson.
%
% Atkinson 2015 GMPE, as defined in the following document. Note that 
% output standard deviations are converted to natural log values from
% the base-10 values reported in the paper.
%
% Source Model:
% Atkinson, G. M. (2015). “Ground-Motion Prediction Equation for 
% Small-to-Moderate Events at Short Hypocentral Distances, with Application 
% to Induced-Seismicity Hazards.” Bulletin of the Seismological Society of 
% America.
% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% INPUT
%   T             = Period (sec); Use Period = -1 for PGV 
%   rup           = rupture object input containing the following
%                   variables:
%       M             = Moment Magnitude
%       Rhyp          = hypocentral distance (km)
%       h_eff         = effective depth (optional). Will be estimated
%                       if not provided by the user
% OUTPUT
%   median                = Median amplitude prediction (units of g)
%   sigma                 = NATURAL LOG standard deviation 
%   period1         = periods for which the median and sigma values are
%                     provided. If T = 1000, then period1 = the full set of
%                     available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(rup.Rhyp)
    disp("Missing Rhyp value")
end

% % Coefficients
% period = [          -1        0       0.03	0.05	0.1     0.2     0.3     0.5     1       2       3       5       ];
% c0	 = [            -4.198	-2.427	-2.313	-2.337	-2.839	-3.918	-2.076	-4.128	-2.009	-4.503	-3.869	-4.374	];
% c1	 = [            1.818    1.877	 1.84	1.902	1.905	2.112	1.889	1.792	1.89	1.532	1.11	1.134	];
% c2	 = [            -0.1009	-0.1214	-0.1119	-0.1252	-0.1134	-0.1266	-0.1257	-0.0791	-0.1248	-0.043	0.0039	0.0038	];
% c3	 = [            -1.721	-1.806	-1.708	-1.838	-1.658	-1.591	-1.886	-1.526	-1.828	-1.404	-1.447	-1.426	];
% c4   = [           -0.0006  -0.002  -0.002  -0.002  -0.002  -0.0014	-0.00105 -0.0006 0       0       0       0      ];
% sigmaIntra	 = [	0.28	0.29	0.29	0.29	0.3     0.31	0.31	0.3     0.27	0.25	0.25	0.26	];
% sigmaInter	 = [	0.18	0.24	0.26	0.29	0.25    0.2     0.18	0.19	0.21	0.22	0.21	0.17	];
% sigmaTotal	 = [	0.33	0.37	0.39	0.41	0.39    0.37	0.36	0.35	0.34	0.33	0.33	0.31	];

period	= [ 	-1	0	0.03	0.05	0.1	0.2	0.3	0.5	1	2	3	5	];
c0	= [ 	-4.151	-2.376	-2.283	-2.018	-1.954	-2.266	-2.794	-3.873	-4.081	-4.462	-3.827	-4.321	];
c1	= [ 	1.762	1.818	1.842	1.826	1.83	1.785	1.852	2.06	1.742	1.485	1.06	1.08	];
c2	= [ 	-0.09509	-0.1153	-0.1189	-0.1192	-0.1185	-0.1061	-0.1078	-0.1212	-0.07381	-0.03815	0.009086	0.009376	];
c3	= [ 	-1.669	-1.752	-1.785	-1.831	-1.774	-1.657	-1.608	-1.544	-1.481	-1.361	-1.398	-1.378	];
c4	= [ 	-0.0006	-0.002	-0.002	-0.002	-0.002	-0.0014	-0.001	-0.0006	0	0	0	0	];
sigmaIntra	= [ 	0.27	0.28	0.28	0.28	0.29	0.3	0.3	0.29	0.26	0.24	0.24	0.25	];
sigmaInter	= [ 	0.19	0.24	0.27	0.3	0.25	0.21	0.19	0.2	0.22	0.23	0.22	0.18	];
sigmaTotal	= [ 	0.33	0.37	0.39	0.41	0.39	0.37	0.36	0.35	0.34	0.33	0.32	0.31	];

if isempty(rup.h_eff) % use estimated effective depth
    h_eff = max(1, 10.^(-0.28+0.19.*rup.M)); % effective depth based on alternative method proposed in Atkinson (2015)
    % rup.h_eff = max(1, 10.^(-1.72+0.43.*M)); % alternate effective depth from Yenier and Atkinson (2014) -- not preferred by Atkinson
end

if length(T) == 1 && T == 1000
    median=zeros(1,length(period)-2);
    sigma=zeros(1,length(period)-2);
    period1=period(3:end);
else
    period1 = T;
    sigma = zeros(1,length(period1));
    median = zeros(1,length(period1));
end

for n = 1:length(period1)
    
% interpolate between periods if neccesary    
if (isempty(find(period == period1(n), 1))) % if the period of interest doesn't match one of the provided periods
    index_low = sum(period<period1(n));
    T_low = period(index_low);
    T_high = period(index_low+1);
    
    % get predictions at closest periods to the period of interest
    [sa_low, sigma_low] = a_2015_active(T_low,rup);
    [sa_high, sigma_high] = a_2015_active(T_high,rup);
    
    % interpolate results
    x = [log(T_low) log(T_high)];
    Y_sa = [log(sa_low) log(sa_high)];
    Y_sigma = [sigma_low sigma_high];
    median(n) = exp(interp1(x,Y_sa,log(period1(n))));
    sigma(n) = interp1(x,Y_sigma,log(period1(n)));
    
else % compute median and sigma for the given index value
    i = find(period == period1(n));
    
    % effective distance
    R = sqrt(rup.Rhyp.^2 + h_eff.^2);
    
    % Compute median and sigma
    logY = c0(i) + c1(i)*rup.M + c2(i)*rup.M.^2 + c3(i).*log10(R) + c4(i).*R;
    if period1(n)>=0
        median(n) = 10.^logY/981; % convert to units of g
    else
        median(n) = 10.^logY;
    end
    sigma(n) = sigmaTotal(i) * log(10) * ones(size(median(n))); % convert to natural log, and make a vector if needed
end
end
end

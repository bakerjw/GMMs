function [rho] = lb_2013_spatial_corr(T1, T2, h)
% Compute the spatial correlation of epsilons for the NGA ground motion models
%
% Created by Christophe Loth, 12/18/2012
% Updated by Jack Baker, 6/28/2019 to resolve an error in the original
% calculation (see Erratum below)
% Updated by Jack Baker, 4/20/2022 to refine the interpolation of B
% matrices and avoid creating saddle points along the diagonal
% Updated by Jack Baker, 5/9/2022 to add significant digits to the B3
% matrix, to ensure that it is positive definite
%
% The function is strictly empirical, fitted over the range 0.01s <= T1, T2 <= 10s
%
% Documentation is provided in the following documents:
% Loth, C., and Baker, J. W. (2013). "A spatial cross-correlation model of 
% ground motion spectral accelerations at multiple periods." 
% Earthquake Engineering & Structural Dynamics, 42(3), 397-417. 
%
% Loth, C., and Baker, J. W. (2020). "Erratum: A spatial cross-correlation 
% model for ground motion spectral accelerations at multiple periods." 
% Earthquake Engineering & Structural Dynamics, 49(3), 315-316.
%
%
% INPUT
%
%   T1, T2      = The two periods of interest. The periods may be equal,
%                 and there is no restriction on which one is larger.
%
%   h           = The separation distance between two sites (units of km)
%
% OUTPUT
%
%   rho         = The predicted correlation coefficient


% Check the validity of input arguments
if min(T1,T2)<0.01
    error('The period must be greater than or equal to 0 s')
end
if max(T1,T2)>10
    error('The periods must be less than or equal to 10 s')
end



Tlist=[0.01 0.1 0.2 0.5 1 2 5 7.5 10.0001];

% matrices of Tlist values corresponding to the B matrix entries
[TT1, TT2] = meshgrid(Tlist, Tlist); 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model coefficients--note that these were updated in 2019 due to an error
% in the original publication. The original publication's coefficients are
% at the bottom of this function.
%
% Table II. Short range coregionalization matrix, B1
B{1}=[0.29	0.25	0.23	0.23	0.18	0.1	0.06	0.06	0.06;
0.25	0.30	0.2	0.16	0.1	0.04	0.03	0.04	0.05;
0.23	0.20	0.27	0.18	0.1	0.03	0	0.01	0.02;
0.23	0.16	0.18	0.31	0.22	0.14	0.08	0.07	0.07;
0.18	0.10	0.1	0.22	0.33	0.24	0.16	0.13	0.12;
0.10	0.04	0.03	0.14	0.24	0.33	0.26	0.21	0.19;
0.06	0.03	0	0.08	0.16	0.26	0.37	0.3	0.26;
0.06	0.04	0.01	0.07	0.13	0.21	0.3	0.28	0.24;
0.06	0.05	0.02	0.07	0.12	0.19	0.26	0.24	0.23];

% Table III. Long range coregionalization matrix, B2
B{2}=[0.47	0.4	0.43	0.35	0.27	0.15	0.13	0.09	0.12;
0.4	0.42	0.37	0.25	0.15	0.03	0.04	0	0.03;
0.43	0.37	0.45	0.36	0.26	0.15	0.09	0.05	0.08;
0.35	0.25	0.36	0.42	0.37	0.29	0.2	0.16	0.16;
0.27	0.15	0.26	0.37	0.48	0.41	0.26	0.21	0.21;
0.15	0.03	0.15	0.29	0.41	0.55	0.37	0.33	0.32;
0.13	0.04	0.09	0.2	0.26	0.37	0.51	0.49	0.49;
0.09	0.00	0.05	0.16	0.21	0.33	0.49	0.62	0.6;
0.12	0.03	0.08	0.16	0.21	0.32	0.49	0.6	0.68];

% Table IV. Nugget effect coregionalization matrix, B3
B{3}=[0.240000000000000	0.219983028675722	0.209991239369580	0.0899940658151642	-0.0199982490874490	0.0100004273375877	0.0299729607606612	0.0200291990885140	0.00995702711846606
0.219983028675722	0.280000000000000	0.199999710563431	0.0400020556476041	-0.0500003168664929	-5.02841885169300e-07	0.0100141900421009	0.00994747690890486	-0.00996663511790072
0.209991239369580	0.199999710563431	0.280000000000000	0.0500007637487926	-0.0600002196848805	-1.80663938055364e-07	0.0399992445541918	0.0299487635912810	0.0100035235224244
0.0899940658151642	0.0400020556476041	0.0500007637487926	0.270000000000000	0.139999321454879	0.0499996979574019	0.0499981807238188	0.0499227563681531	0.0399858842409999
-0.0199982490874490	-0.0500003168664929	-0.0600002196848805	0.139999321454879	0.190000000000000	0.0700000354290215	0.0499897414826383	0.0499443288879162	0.0499652709189241
0.0100004273375877	-5.02841885169300e-07	-1.80663938055364e-07	0.0499996979574019	0.0700000354290215	0.120000000000000	0.0799859494118349	0.0699172702775962	0.0599608152312721
0.0299729607606612	0.0100141900421009	0.0399992445541918	0.0499981807238188	0.0499897414826383	0.0799859494118349	0.120000000000000	0.0997643834727755	0.0800031285024676
0.0200291990885140	0.00994747690890486	0.0299487635912810	0.0499227563681531	0.0499443288879162	0.0699172702775962	0.0997643834727755	0.100000000000000	0.0896690207890228
0.00995702711846606	-0.00996663511790072	0.0100035235224244	0.0399858842409999	0.0499652709189241	0.0599608152312721	0.0800031285024676	0.0896690207890228	0.0900000000000000];



% Find the interval containing each input period 
index1 = find((Tlist<=T1),1,'last');
index2 = find((Tlist<=T2),1,'last');


% Interpolate each coregionalization matrix coefficient
if index1 == index2 
    % period pair is close to the diagonal. Uses a careful interpolation to 
    % keep the ridgeo on the diagonal and avoid creating saddle points 
    for i=1:3
        Bcoeff{i} = interpolate_B(Tlist, B{i}, T1, T2);
    end

else 
    % linearly interpolate using the nearest four points
    for i=1:3
        Bcoeff{i} = griddata(TT1(index1:index1+1,index2:index2+1), ...
                             TT2(index1:index1+1,index2:index2+1), ...
                             B{i}(index1:index1+1,index2:index2+1), T2, T1);
        % the indexing above is to pass in only the adjacent entries of the
        % matrix, to speed the interpolation calculation
    end
end


% Compute the correlation coefficient (Equation 42)
rho=Bcoeff{1}*exp(-3*h/20) + Bcoeff{2}*exp(-3*h/70) + Bcoeff{3}*(h==0);


end



function Bcoeff = interpolate_B(Tlist, B, T1, T2)

% interpolate the matrix B, recognizing that the diagonal of the matrix is
% the peak

% Find the interval containing each input period 
index1 = find((Tlist<=T1),1,'last');
index2 = find((Tlist<=T2),1,'last');

% take just the adjacent cells of the B matrix
T1vals = Tlist(index1:index1+1);
T2vals = Tlist(index2:index2+1);
Bvals = B(index1:index1+1, index2:index2+1);

% interpolate along the diagonal
Tavgtarg = mean([T1; T2]);

Tavgvals = mean([T1vals; T2vals]);
Bavgvals = diag(Bvals);
Bdiagval = interp1(Tavgvals,Bavgvals,Tavgtarg);

% then interpolate between the diagonal and the corner value of B
Tdifftarg = abs(T1-T2);

Tdiffvals = [0 abs(min(T1vals) - max(T1vals))];
Bdiffvals = [Bdiagval Bvals(1,2)];
Bcoeff = interp1(Tdiffvals, Bdiffvals, Tdifftarg);


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Below are the original tables from the 2013 manuscript, but they are 
% superceded due to an error in the original calculations. 
% Uncomment these and move them up in the function to reproduce the 
% original paper's results.

% % Table II. Short range coregionalization matrix, B{1}
% B{1}=[0.30 0.24 0.23 0.22 0.16 0.07 0.03 0 0;
% 0.24 0.27 0.19 0.13 0.08 0 0 0 0;
% 0.23 0.19 0.26 0.19 0.12 0.04 0 0 0;
% 0.22 0.13 0.19 0.32 0.23 0.14 0.09 0.06 0.04;
% 0.16 0.08 0.12 0.23 0.32 0.22 0.13 0.09 0.07;
% 0.07 0 0.04 0.14 0.22 0.33 0.23 0.19 0.16;
% 0.03 0 0 0.09 0.13 0.23 0.34 0.29 0.24;
% 0 0 0 0.06 0.09 0.19 0.29 0.30 0.25;
% 0 0 0 0.04 0.07 0.16 0.24 0.25 0.24];
% 
% % Table III. Long range coregionalization matrix, B{2}
% B{2}=[0.31 0.26 0.27 0.24 0.17 0.11 0.08 0.06 0.05;
% 0.26 0.29 0.22 0.15 0.07 0 0 0 -0.03;
% 0.27 0.22 0.29 0.24 0.15 0.09 0.03 0.02 0;
% 0.24 0.15 0.24 0.33 0.27 0.23 0.17 0.14 0.14;
% 0.17 0.07 0.15 0.27 0.38 0.34 0.23 0.19 0.21;
% 0.11 0 0.09 0.23 0.34 0.44 0.33 0.29 0.32;
% 0.08 0 0.03 0.17 0.23 0.33 0.45 0.42 0.42;
% 0.06 0 0.02 0.14 0.19 0.29 0.42 0.47 0.47;
% 0.05 -0.03 0 0.14 0.21 0.32 0.42 0.47 0.54];
% 
% % Table IV. Nugget effect coregionalization matrix, B{3}
% B{3}=[0.38 0.36 0.35 0.17 0.04 0.04 0 0.03 0.08;
% 0.36 0.43 0.35 0.13 0 0.02 0 0.02 0.08;
% 0.35 0.35 0.45 0.11 -0.04 -0.02 -0.04 -0.02 0.03;
% 0.17 0.13 0.11 0.35 0.2 0.06 0.02 0.04 0.02;
% 0.04 0 -0.04 0.20 0.30 0.14 0.09 0.12 0.04;
% 0.04 0.02 -0.02 0.06 0.14 0.22 0.12 0.13 0.09;
% 0 0 -0.04 0.02 0.09 0.12 0.21 0.17 0.13;
% 0.03 0.02 -0.02 0.04 0.12 0.13 0.17 0.23 0.10;
% 0.08 0.08 0.03 0.02 0.04 0.09 0.13 0.10 0.22];



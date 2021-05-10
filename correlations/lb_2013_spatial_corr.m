function [rho] = lb_2013_spatial_corr(T1, T2, h)

%
% Created by Christophe Loth, 12/18/2012
% Updated by Jack Baker, 6/28/2019 to resolve an error in the original
% calculation (see Erratum below)
%
% Compute the spatial correlation of epsilons for the NGA ground motion models
%
% The function is strictly empirical, fitted over the range 0.01s <= T1, T2 <= 10s
%
% Documentation is provided in the following document:
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


% Verify the validity of input arguments
if min(T1,T2)<0.01
    error('The period must be greater than or equal to 0 s')
end
if max(T1,T2)>10
    error('The periods must be less than or equal to 10 s')
end



Tlist=[0.01 0.1 0.2 0.5 1 2 5 7.5 10.0001];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Model coefficients--note that these were updated due to an error
% in the original publication. The original publication's coefficients are
% at the bottom of this function.
%
% Table II. Short range coregionalization matrix, B1
B1=[0.29	0.25	0.23	0.23	0.18	0.1	0.06	0.06	0.06;
0.25	0.30	0.2	0.16	0.1	0.04	0.03	0.04	0.05;
0.23	0.20	0.27	0.18	0.1	0.03	0	0.01	0.02;
0.23	0.16	0.18	0.31	0.22	0.14	0.08	0.07	0.07;
0.18	0.10	0.1	0.22	0.33	0.24	0.16	0.13	0.12;
0.10	0.04	0.03	0.14	0.24	0.33	0.26	0.21	0.19;
0.06	0.03	0	0.08	0.16	0.26	0.37	0.3	0.26;
0.06	0.04	0.01	0.07	0.13	0.21	0.3	0.28	0.24;
0.06	0.05	0.02	0.07	0.12	0.19	0.26	0.24	0.23];

% Table III. Long range coregionalization matrix, B2
B2=[0.47	0.4	0.43	0.35	0.27	0.15	0.13	0.09	0.12;
0.4	0.42	0.37	0.25	0.15	0.03	0.04	0	0.03;
0.43	0.37	0.45	0.36	0.26	0.15	0.09	0.05	0.08;
0.35	0.25	0.36	0.42	0.37	0.29	0.2	0.16	0.16;
0.27	0.15	0.26	0.37	0.48	0.41	0.26	0.21	0.21;
0.15	0.03	0.15	0.29	0.41	0.55	0.37	0.33	0.32;
0.13	0.04	0.09	0.2	0.26	0.37	0.51	0.49	0.49;
0.09	0.00	0.05	0.16	0.21	0.33	0.49	0.62	0.6;
0.12	0.03	0.08	0.16	0.21	0.32	0.49	0.6	0.68];

% Table IV. Nugget effect coregionalization matrix, B3
B3=[0.24	0.22	0.21	0.09	-0.02	0.01	0.03	0.02	0.01;
0.22	0.28	0.2	0.04	-0.05	0	0.01	0.01	-0.01;
0.21	0.20	0.28	0.05	-0.06	0	0.04	0.03	0.01;
0.09	0.04	0.05	0.26	0.14	0.05	0.05	0.05	0.04;
-0.02	-0.05	-0.06	0.14	0.20	0.07	0.05	0.05	0.05;
0.01	0.00	0.00	0.05	0.07	0.12	0.08	0.07	0.06;
0.03	0.01	0.04	0.05	0.05	0.08	0.12	0.1	0.08;
0.02	0.01	0.03	0.05	0.05	0.07	0.1	0.1	0.09;
0.01	-0.01	0.01	0.04	0.05	0.06	0.08	0.09	0.09];




% Find in which interval each input period is located
for i=1:length(Tlist)-1
    if T1-Tlist(i)>=0 & T1-Tlist(i+1)<0
        index1=i;
    end
    if T2-Tlist(i)>=0 & T2-Tlist(i+1)<0
        index2=i;
    end
end


% Linearly interpolate the corresponding value of each coregionalization
% matrix coefficient

B1coeff1=B1(index1,index2)+(B1(index1+1,index2)-B1(index1,index2))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B1coeff2=B1(index1,index2+1)+(B1(index1+1,index2+1)-B1(index1,index2+1))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B1coeff0=B1coeff1+(B1coeff2-B1coeff1)./...
    (Tlist(index2+1)-Tlist(index2)).*(T2-Tlist(index2));




B2coeff1=B2(index1,index2)+(B2(index1+1,index2)-B2(index1,index2))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B2coeff2=B2(index1,index2+1)+(B2(index1+1,index2+1)-B2(index1,index2+1))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B2coeff0=B2coeff1+(B2coeff2-B2coeff1)./...
    (Tlist(index2+1)-Tlist(index2)).*(T2-Tlist(index2));




B3coeff1=B3(index1,index2)+(B3(index1+1,index2)-B3(index1,index2))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B3coeff2=B3(index1,index2+1)+(B3(index1+1,index2+1)-B3(index1,index2+1))./...
    (Tlist(index1+1)-Tlist(index1)).*(T1-Tlist(index1));

B3coeff0=B3coeff1+(B3coeff2-B3coeff1)./...
    (Tlist(index2+1)-Tlist(index2)).*(T2-Tlist(index2));


% Compute the correlation coefficient (Equation 42)

rho=B1coeff0*exp(-3*h/20)+B2coeff0*exp(-3*h/70);

if h==0
    rho=B1coeff0*exp(-3*h/20)+B2coeff0*exp(-3*h/70)+B3coeff0;
end


end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Below are the original tables from the 2013 manuscript, but they are 
% superceded due to an error in the original calculations. 
% Uncomment these and move them up in the function to reproduce the 
% original paper's results.

% % Table II. Short range coregionalization matrix, B1
% B1=[0.30 0.24 0.23 0.22 0.16 0.07 0.03 0 0;
% 0.24 0.27 0.19 0.13 0.08 0 0 0 0;
% 0.23 0.19 0.26 0.19 0.12 0.04 0 0 0;
% 0.22 0.13 0.19 0.32 0.23 0.14 0.09 0.06 0.04;
% 0.16 0.08 0.12 0.23 0.32 0.22 0.13 0.09 0.07;
% 0.07 0 0.04 0.14 0.22 0.33 0.23 0.19 0.16;
% 0.03 0 0 0.09 0.13 0.23 0.34 0.29 0.24;
% 0 0 0 0.06 0.09 0.19 0.29 0.30 0.25;
% 0 0 0 0.04 0.07 0.16 0.24 0.25 0.24];
% 
% % Table III. Long range coregionalization matrix, B2
% B2=[0.31 0.26 0.27 0.24 0.17 0.11 0.08 0.06 0.05;
% 0.26 0.29 0.22 0.15 0.07 0 0 0 -0.03;
% 0.27 0.22 0.29 0.24 0.15 0.09 0.03 0.02 0;
% 0.24 0.15 0.24 0.33 0.27 0.23 0.17 0.14 0.14;
% 0.17 0.07 0.15 0.27 0.38 0.34 0.23 0.19 0.21;
% 0.11 0 0.09 0.23 0.34 0.44 0.33 0.29 0.32;
% 0.08 0 0.03 0.17 0.23 0.33 0.45 0.42 0.42;
% 0.06 0 0.02 0.14 0.19 0.29 0.42 0.47 0.47;
% 0.05 -0.03 0 0.14 0.21 0.32 0.42 0.47 0.54];
% 
% % Table IV. Nugget effect coregionalization matrix, B3
% B3=[0.38 0.36 0.35 0.17 0.04 0.04 0 0.03 0.08;
% 0.36 0.43 0.35 0.13 0 0.02 0 0.02 0.08;
% 0.35 0.35 0.45 0.11 -0.04 -0.02 -0.04 -0.02 0.03;
% 0.17 0.13 0.11 0.35 0.2 0.06 0.02 0.04 0.02;
% 0.04 0 -0.04 0.20 0.30 0.14 0.09 0.12 0.04;
% 0.04 0.02 -0.02 0.06 0.14 0.22 0.12 0.13 0.09;
% 0 0 -0.04 0.02 0.09 0.12 0.21 0.17 0.13;
% 0.03 0.02 -0.02 0.04 0.12 0.13 0.17 0.23 0.10;
% 0.08 0.08 0.03 0.02 0.04 0.09 0.13 0.10 0.22];



function [ rho ] = a_2011_corr( T1, T2 )
% Created by Jack Baker, circa 2015
% Updated 5/7/2021 to improve documentation and allow vector inputs
%
%
% Horizontal spectral correlation coefficients from subduction motions.
% Documentation at:
%
% Al Atik, L. (2011). "Correlation of spectral acceleration values for 
% subduction and crustal models, COSMOS Technical Session. Emeryville, 
% California, 4 November 2011.
% 
% Carlton, B., and Abrahamson, N. (2014). “Issues and Approaches for 
% Implementing Conditional Mean Spectra in Practice.” Bulletin of the 
% Seismological Society of America, 104(1), 503–512.
%
% The function is fitted over the range the range 0.01s <= T1, T2 <= 5s
%
% INPUT
%
%   T1, T2      = The two periods of interest. The periods may be equal,
%                 and there is no restriction on which one is larger.
%
% OUTPUT
%
%   rho         = The predicted correlation matrix. If length(T1)=n and
%                   length(T2)=m, then size(rho) = [n m]


%% tabulated model parameters

% Periods reported in Table 6 (PGA mapped to 0.001s)
Periods = [0.01	0.05	0.075	0.1	0.15	0.2	0.25	0.3	0.4	0.5	0.6	0.75	1	1.5	2	2.5	3	4	5];

% digitized correlation coefficients from Table 6
RHO = [1	0.9594	0.9376	0.9303	0.93	0.9102	0.8866	0.8579	0.8169	0.7289	0.6482	0.5801	0.4641	0.3913	0.309	0.2744	0.2795	0.2231	0.2307; ...
0.9594	1	0.971	0.9415	0.8996	0.8548	0.8198	0.783	0.7166	0.6342	0.5478	0.4768	0.4085	0.2995	0.218	0.2028	0.1865	0.1487	0.1679; ...
0.9376	0.971	1	0.9677	0.9089	0.8487	0.7984	0.7543	0.6795	0.5842	0.4975	0.427	0.3532	0.2471	0.1751	0.1618	0.1558	0.1208	0.1347; ...
0.9303	0.9415	0.9677	1	0.9346	0.8737	0.8058	0.7534	0.6988	0.5761	0.4869	0.4159	0.2836	0.2295	0.1669	0.1462	0.1622	0.1123	0.1313; ...
0.93	0.8996	0.9089	0.9346	1	0.9338	0.8702	0.8124	0.7246	0.624	0.5358	0.4643	0.3842	0.2693	0.1903	0.1746	0.1709	0.1404	0.1468; ...
0.9102	0.8548	0.8487	0.8737	0.9338	1	0.9411	0.8847	0.8101	0.6998	0.6138	0.5354	0.4188	0.3256	0.278	0.2282	0.2452	0.181	0.184; ...
0.8866	0.8198	0.7984	0.8058	0.8702	0.9411	1	0.9516	0.8601	0.7659	0.6785	0.5977	0.5028	0.383	0.301	0.2763	0.256	0.2248	0.22; ...
0.8579	0.783	0.7543	0.7534	0.8124	0.8847	0.9516	1	0.9174	0.8256	0.7427	0.6643	0.5636	0.4369	0.3542	0.3278	0.3007	0.2694	0.2549; ...
0.8169	0.7166	0.6795	0.6988	0.7246	0.8101	0.8601	0.9174	1	0.9213	0.8425	0.7619	0.6272	0.5306	0.4328	0.4118	0.3713	0.3315	0.3061; ...
0.7289	0.6342	0.5842	0.5761	0.624	0.6998	0.7659	0.8256	0.9213	1	0.9392	0.8637	0.7561	0.6156	0.5217	0.4832	0.4361	0.3919	0.3637; ...
0.6482	0.5478	0.4975	0.4869	0.5358	0.6138	0.6785	0.7427	0.8425	0.9392	1	0.9307	0.8176	0.6817	0.5811	0.5316	0.4812	0.4357	0.4; ...
0.5801	0.4768	0.427	0.4159	0.4643	0.5354	0.5977	0.6643	0.7619	0.8637	0.9307	1	0.9043	0.772	0.6765	0.6241	0.5707	0.5264	0.4983; ...
0.4641	0.4085	0.3532	0.2836	0.3842	0.4188	0.5028	0.5636	0.6272	0.7561	0.8176	0.9043	1	0.8793	0.7771	0.7287	0.6714	0.6288	0.5863; ...
0.3913	0.2995	0.2471	0.2295	0.2693	0.3256	0.383	0.4369	0.5306	0.6156	0.6817	0.772	0.8793	1	0.9243	0.8663	0.8156	0.754	0.7093; ...
0.309	0.218	0.1751	0.1669	0.1903	0.278	0.301	0.3542	0.4328	0.5217	0.5811	0.6765	0.7771	0.9243	1	0.9474	0.901	0.8455	0.7984; ...
0.2744	0.2028	0.1618	0.1462	0.1746	0.2282	0.2763	0.3278	0.4118	0.4832	0.5316	0.6241	0.7287	0.8663	0.9474	1	0.9653	0.901	0.8609; ...
0.2795	0.1865	0.1558	0.1622	0.1709	0.2452	0.256	0.3007	0.3713	0.4361	0.4812	0.5707	0.6714	0.8156	0.901	0.9653	1	0.9463	0.9062; ...
0.2231	0.1487	0.1208	0.1123	0.1404	0.181	0.2248	0.2694	0.3315	0.3919	0.4357	0.5264	0.6288	0.754	0.8455	0.901	0.9463	1	0.9611; ...
0.2307	0.1679	0.1347	0.1313	0.1468	0.184	0.22	0.2549	0.3061	0.3637	0.4	0.4983	0.5863	0.7093	0.7984	0.8609	0.9062	0.9611	1];

%% compute correlations
% for i = 1:length(T1)
%     for j = 1:length(T2)
%           
%         if (min([T1(i) T2(j)])<0.01 || max([T1(i) T2(j)])>5)
%             % periods are out of allowable range
%             rho(i,j) = nan;
%             
%             
%         else % periods are in range--compute correlations
%             
%             rho(i,j) = griddata(log(Periods), log(Periods), RHO, log(T1(i)), log(T2(j)));
%             
%         end
%     end
% end
        
rho = nan*ones(length(T1), length(T2)); % initialize nans

% find indices of input periods that are within range
t1Idx = find(T1>=0.01 & T1<=5);
t2Idx = find(T2>=0.01 & T2<=5);

[X, Y] = meshgrid(log(T2(t2Idx)), log(T1(t1Idx)));

% compute correlations at periods that are within range
rho(t1Idx, t2Idx) = interp2(log(Periods), log(Periods), RHO, X, Y);


end


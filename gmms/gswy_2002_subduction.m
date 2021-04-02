function [median, sigma, period1] = gswy_2002_subduction(T,rup,site)

% Created by Jack Baker, 3/29/2016
% Updated by Emily Mongold, 11/27/20
%
% Predict ground motions for Cascadia Megathrust earthquakes
%
% Source Model:
% Gregor, N. J., Silva, W. J., Wong, I. G., and Youngs, R. R. (2002).
% "Ground-Motion Attenuation Relationships for Cascadia Subduction Zone
% Megathrust Earthquakes Based on a Stochastic Finite-Fault Model."
% Bulletin of the Seismological Society of America, 92(5), 1923-1932.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%INPUT
%   T           = Period of Interest (sec). For PGA Computation, T = 0.001,
%                 for the whole period spectra, T = 1000.
%   rup         = rupture object input containing the following
%                 variables:
%       M           = Moment Magnitude
%       R           = Source to site distance to rupture surface (km)
%   site        = site object input containing the following
%                 variable:
%       is_soil     = Soil type indicator (Zr used locally: 
%                        For Rock Sites, Zr = 1
%                        Otherwise, Zr = 0)
%OUTPUT
%   median      = Median spectral acceleration prediction (g)
%   sigma       = Logarithmic standard deviation of spectral acceleration
%                 prediction
%   period1     = periods for which the median and sigma values are
%                 provided. If T = 1000, then period1 = the full set of
%                 available periods. Else period1 = T
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if site.is_soil == 0 % Mapping generic soil inputs to model-specific values
    Zr = 0; 
else
    Zr = 1;
end

% get coefficients
[period, C1, C2, C3, C4, C5, C6, TotalSigma] = fn_get_coefficients(Zr);

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

% iterate in the case that T is a vector
for n = 1:length(period1)

i = find(abs(period - period1(n)) < 1e-4); % find appropriate period index

if isempty(i) % period not included, so interpolate
    
    i_low = sum(period < period1(n));
    T_low = period(i_low);
    T_high = period(i_low + 1);
    
    [Sa_high, sigma_high,~] = gswy_2002_subduction(T_high,rup,site);
    [Sa_low, sigma_low,~]   = gswy_2002_subduction(T_low,rup,site);
    
    % interpolate
    median(n) = exp(interp1(log([T_low T_high]),log([Sa_low Sa_high]),log(period1(n))));
    sigma(n) = exp(interp1(log([T_low T_high]),log([sigma_low sigma_high]),log(period1(n))));
    
else % evaluate using the appropriate index
    
    ln_Y = C1(i) + C2(i)*rup.M + ( C3(i) +  C4(i)*rup.M)*log(rup.R+exp(C5(i))) ...
        + C6(i)*(rup.M-10).^3;
    
    sigma(n) = TotalSigma(i);
    median(n) = exp(ln_Y);
    
end
end
end


% % LOCAL FUNCTION % %

function [period, C1, C2, C3, C4, C5, C6, TotalSigma] = fn_get_coefficients(Zr)
% return the appropriate set of coefficients for a rock or soil site

if Zr == 1
    
    %------------------------- Rock Coefficients (from Table 2)-----------------------%
    period = [  	0.001	0.01	0.02	0.025	0.032	0.04	0.05	0.056	0.063	0.071	0.083	0.1	0.125	0.143	0.167	0.2	0.25	0.333	0.4	0.5	0.769	1	1.667	2	2.5	5];
    C1 = [ 	21.0686	20.9932	21.072	21.152	21.366	17.525	19.347	20.774	21.331	24.221	24.95	30.005	39.719	43.414	39.579	39.345	37.69	34.787	33.393	29.159	15.279	6.528	7.467	8.657	6.637	8.013];
    C2 = [ 	-1.7712	-1.7658	-1.772	-1.779	-1.797	-1.339	-1.519	-1.625	-1.672	-1.924	-1.979	-2.349	-3.09	-3.385	-2.957	-3.087	-2.96	-2.899	-2.776	-2.424	-1.22	-0.406	-0.676	-0.851	-0.651	-0.943];
    C3 = [ 	-5.0631	-5.0404	-5.0529	-5.0663	-5.1036	-4.8602	-4.9731	-5.1875	-5.2561	-5.625	-5.6696	-6.3862	-7.8541	-8.3122	-7.9723	-7.6002	-7.379	-6.7855	-6.9595	-6.2114	-4.324	-3.1991	-2.6465	-2.7398	-2.3124	-2.4087];
    C4 = [ 	0.4153	0.4132	0.4142	0.4154	0.4187	0.3868	0.396	0.4118	0.4173	0.4478	0.4493	0.5009	0.6161	0.6513	0.6139	0.5972	0.5842	0.5616	0.5863	0.5216	0.3618	0.2589	0.2193	0.2339	0.1879	0.2154];
    C5 = [ 	4.2	4.2	4.2	4.2	4.2	4.2	4.2	4.3	4.3	4.4	4.4	4.7	5.1	5.2	5.2	5.1	5.1	4.9	4.9	4.7	3.9	3.2	2.8	2.8	2.8	2.3];
    C6 = [ 	0.0017	0.0226	0.0025	0.0023	0.0017	-0.0318	-0.0155	-0.0155	-0.0146	-0.0071	-0.0018	-0.0019	-0.0064	-0.0001	-0.0264	0.006	-0.0023	0.0256	-0.0039	0.0161	-0.0011	-0.0225	0.0416	0.037	0.0364	0.0647];
    ParamSigma = [ 	0.6083	0.6031	0.6036	0.6042	0.6062	0.5836	0.5908	0.5974	0.6028	0.6116	0.6337	0.6448	0.6654	0.6769	0.681	0.7034	0.7121	0.7372	0.711	0.6745	0.6111	0.5898	0.4931	0.4666	0.4163	0.3931];
    ModelSigma = [	0.3926	0.3926	0.3926	0.3983	0.3926	0.3818	0.3925	0.4052	0.4132	0.4042	0.4584	0.4668	0.5461	0.5225	0.505	0.5089	0.4539	0.4764	0.5187	0.4382	0.5611	0.4751	0.4889	0.4247	0.5198	0.6656];
    TotalSigma = [	0.724	0.7195	0.7195	0.7235	0.7221	0.6969	0.7086	0.7215	0.7302	0.7326	0.7815	0.7954	0.8605	0.8544	0.8478	0.8679	0.8444	0.8776	0.8801	0.8039	0.8295	0.7567	0.6943	0.6305	0.6657	0.773];
    
elseif Zr == 0
    
    %------------------------- Soil Coefficients (from Table 3)-----------------------%
    period = [0.001	0.01	0.02	0.025	0.032	0.04	0.05	0.056	0.0625	0.07	0.083	0.1	0.125	0.143	0.167	0.2	0.25	0.33	0.4	0.5	0.77	1	1.67	2	2.5	5];
    C1  = [	23.8613	25.4516	25.4339	25.42	25.3849	22.7042	23.2948	23.2165	24.7067	24.9425	26.5395	29.9693	35.666	50.7368	55.6402	75.8218	100.3357	71.7967	67.372	56.0088	26.3013	17.233	11.9971	17.9124	16.1666	7.4856];
    C2  = [	-2.2742	-2.4206	-2.4185	-2.4168	-2.4127	-2.1004	-2.1619	-2.1528	-2.2814	-2.3045	-2.4402	-2.7254	-3.1853	-4.5292	-4.9662	-6.8396	-9.0324	-6.499	-6.1755	-5.1176	-2.4482	-1.5506	-1.118	-1.7505	-1.5091	-0.836];
    C3  = [	-4.8803	-5.1071	-5.1044	-5.1026	-5.0977	-4.9006	-4.8855	-4.8744	-5.0947	-5.0672	-5.3025	-5.8054	-6.6251	-8.7213	-9.5555	-12.0687	-15.3511	-11.6056	-11.1567	-9.5083	-5.3818	-4.3287	-2.9451	-3.815	-3.7101	-2.0627];
    C4  = [	0.4399	0.4605	0.4602	0.46	0.4594	0.4353	0.4332	0.4319	0.4509	0.4476	0.4677	0.5098	0.5769	0.7649	0.8435	1.0753	1.3731	1.0415	1.0167	0.8632	0.4957	0.393	0.2639	0.3574	0.3344	0.1779];
    C5  = [	4.7	4.8	4.8	4.8	4.8	4.8	4.8	4.8	4.9	4.9	5	5.2	5.5	5.9	6	6.3	6.6	6.2	6.1	5.9	4.8	4.2	3.7	4.1	4.1	-0.2];
    C6  = [	0.0366	0.0372	0.037	0.0369	0.0366	0.0164	0.0263	0.0255	0.0245	0.0295	0.0276	0.0226	0.0123	0.0108	-0.007	0.0096	-0.0043	0.0102	0.0035	0.0164	0.0259	0.0133	0.0538	0.0583	0.0473	0.0821];
    ParamSigma  = [	0.376	0.3742	0.3742	0.3743	0.3743	0.359	0.3592	0.3598	0.3607	0.3609	0.3617	0.3654	0.3821	0.3923	0.3927	0.4231	0.4472	0.4324	0.4243	0.4305	0.4601	0.4599	0.4781	0.4628	0.4193	0.4802];
    ModelSigma  = [	0.3926	0.3926	0.3926	0.3983	0.3926	0.3818	0.3925	0.4052	0.4132	0.4042	0.4584	0.4668	0.5461	0.5225	0.505	0.5089	0.4539	0.4764	0.5187	0.4382	0.5611	0.4751	0.4889	0.4247	0.5198	0.6656];
    TotalSigma	 = [0.5436	0.5422	0.5422	0.5464	0.5422	0.5241	0.5319	0.5413	0.548	0.5413	0.5835	0.5926	0.6665	0.6532	0.6393	0.6618	0.6371	0.6431	0.6699	0.6139	0.7256	0.6606	0.6837	0.6276	0.6676	0.8207];
    
end

end

function [rho] = hm_2019_spatial_corr(T, h)

% Created by Jack Baker, 5/10/2021
%
% Compute the spatial correlation of epsilons. Documentation is provided 
% in the following document:
%
% Heresi, P., and Miranda, E. (2019). "Uncertainty in intraevent spatial 
% correlation of elastic pseudo-acceleration spectral ordinates." Bulletin 
% of Earthquake Engineering, 17(3), 1099-1115.
%
% The model was fitted over the range 0 <= T <= 10s.
% The paper proposes a model that randomly samples correlation functions
% for a given rupture. Only the median parameters are implemented in this
% function, for ease of comparisons with other models.
%
% INPUT
%
%   T               = The period of interest. 
%
%   h               = The separation distance between two sites (units of km)
%
% OUTPUT
%
%   rho             = The predicted correlation coefficient



% Check input parameters
if T<0
    error('The period must be greater than or equal to 0 s')
end
if T>10
    error('The periods must be less than or equal to 10 s')
end

% compute beta parameter (equation 9)
if T<1.37
    beta = 4.231*T^2 - 5.180*T + 13.392;
else
    beta = 0.140*T^2 - 2.249*T + 17.050;
end
    
% compute correlation (equation 8)
rho = exp( -(h/beta).^0.55 );

end
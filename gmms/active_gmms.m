function [median,sigma, period1] = active_gmms(T,rup,site,GMM_name)
%ACTIVE_CALLER Function to call active crustal GMM functions
%   Created by Emily Mongold 4/23/2021
%   
%   Inputs:
%       T = Period
%       rup = rup object 
%       site = site object
%       GMM_name = one of the following options:
%                       'a_2015'
%                       'as_1997'
%                       'as_2008'
%                       'ask_2014'
%                       'ba_2008'
%                       'bjf_1997'
%                       'bssa_2014'
%                       'cb_2008'
%                       'cb_2014'
%                       'cy_2008'
%                       'cy_2014'
%                       'i_2008'
%                       'i_2014'
%                       'scemy_1997'
%                       'z_2006'
% 

% addpath('./gmms/')
% addpath('./testing/')

switch GMM_name
    
    case 'a_2015'
        [median,sigma, period1] = a_2015_active(T,rup);
    
    case 'as_1997'
        if T == 0
            [median, sigma, period1] = as_1997_active(0.01,rup,site);
        else
            [median, sigma, period1] = as_1997_active(T,rup,site);
        end
    
    case 'as_2008'
        [median, sigma, period1] = as_2008_active(T,rup,site);
    
    case 'ask_2014'
        [median, sigma, period1] = ask_2014_active(T,rup,site);
    
    case 'ba_2008'
        [median, sigma, period1] = ba_2008_active(T, rup, site);
        
    case 'bjf_1997'
        if T == 0
            [median, sigma, period1] = bjf_1997_active(0.001,rup,site);
        else
            [median, sigma, period1] = bjf_1997_active(T,rup,site);
        end
        
    case 'bssa_2014'
        [median, sigma, period1] = bssa_2014_active(T, rup, site);
        
    case 'cb_2008'
        [median, sigma, period1] = cb_2008_active(T,rup,site);
        
    case 'cb_2014'
        [median, sigma, period1] = cb_2014_active(T,rup,site);
        
    case 'cy_2008'
        [median, sigma, period1] = cy_2008_active(T,rup,site);

    case 'cy_2014'
        [median, sigma, period1] = cy_2014_active(T,rup,site);
        
    case 'i_2008'
        if T == 0
            [median, sigma, period1] = i_2008_active(0.01,rup,site);
        else
            [median, sigma, period1] = i_2008_active(T,rup,site);
        end
        
    case 'i_2014'
        [median, sigma, period1] = i_2014_active(T,rup,site);
        
    case 'scemy_1997'
        if T == 0
            [median, sigma, period1] = scemy_1997_active(0.001,rup,site);
        else
            [median, sigma, period1] = scemy_1997_active(T,rup,site);
        end
        
    case 'z_2006'
        [median, sigma, period1] = z_2006_active(T,rup,site,[],0);
        
end


end


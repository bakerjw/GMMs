%% Define site class
classdef site < handle
    % Defines site class and properties
    properties
        is_soil % Flag for soil conditions
%               = 0 for soil
%               = 1 for soft rock
%               = 2 for hard rock
        Vs30 % Average shear wave velocity over top 30 m
        fvs30 % Flag for shear wave velocity
%               = 0 for Vs30 inferred from geology
%               = 1 for measured Vs30
        Z25 % Depth to Vs = 2.5km/s
        Z10 % Depth to Vs = 1.0km/s
        Zbot % Depth to bottom of the seismogenic crust
        region % Country
%               = 0 for global
%               = 1 for California
%               = 2 for Japan
%               = 3 for China
%               = 4 for Italy
%               = 5 for Turkey
%               = 6 for Taiwan
    end
        methods
        function self = site(is_soil,Vs30,fvs30,Z25,Z10,Zbot,region)
            self.is_soil = is_soil; % Flag for soil conditions
            self.Vs30 = Vs30; % Average shear wave velocity over top 30 m
            self.fvs30 = fvs30; % Flag for measured or inferred Vs30
            self.Z25 = Z25; % Depth to Vs = 2.5km/s
            self.Z10 = Z10; % Depth to Vs = 1.0km/s
            self.Zbot = Zbot; % Depth to bottom of seismogenic crust
            self.region = region; % Country
        end
    
    end
    
end
   
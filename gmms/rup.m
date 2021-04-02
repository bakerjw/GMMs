%% Define rupture class
classdef rup < handle
    % Defines rupture class and properties
    properties
        M % Magnitude
        R % Generic distance
        Rrup % Closest distance to fault rupture
        Rjb % Joyner-Boore distance
        Rhyp % Hypocentral distance
        Rx % Perpendicular distance to projection of rupture edge
        Ry0 % Parallel distance off end of rupture
        HW % Hanging wall indicator = 1 for hanging wall, 0 otherwise
        AS % Aftershock indicator   = 1 for aftershock, 0 otherwise
        Ztor % Depth to top of rupture
        Zhyp % Hypocentral depth
        h_eff % Effective height
        W % Down-dip width
        delta % Average dip angle
        lambda % Rake angle (used to determine fault types)
    end
    
    methods
        function self = rup(M,R,Rrup,Rjb,Rhyp,Rx,Ry0,HW,AS,Ztor,Zhyp,h_eff,W,delta,lambda)
            self.M = M; % Magnitude
            self.R = R; % Generic distance
            self.Rrup = Rrup; % Closest distance to fault rupture
            self.Rjb = Rjb; % Joyner-Boore distance
            self.Rhyp = Rhyp; % Hypocentral distance
            self.Rx = Rx; % Perpendicular distance from surface projection of undip edge of rupture
            self.Ry0 = Ry0; % Horizontal distance off end of rupture measured parallel to strike
            self.HW = HW; % Hanging wall indicator
            self.AS = AS; % Aftershock indicator
            self.Ztor = Ztor; % Depth to top of rupture
            self.Zhyp = Zhyp; % Hypocentral depth
            self.h_eff = h_eff; % Effective height
            self.W = W; % Down-dip rupture width
            self.delta = delta; % Average dip angle
            self.lambda = lambda; % Rake angle
        end
    end
    
end
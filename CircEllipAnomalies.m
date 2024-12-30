function [M, E]= CircEllipAnomalies(e, TA)

%   Converts true anomaly to mean and eccentric anomaly.
%
%   INPUT:  e -- eccentricity
%          TA -- deg, True Anomaly of Elliptical or Circular Orbit
%
%   OUTPUT: M -- deg, Mean Anomaly
%           E -- deg, Eccentric Anomaly
%



E = 2*atan2d(sqrt(1-e)*tand(TA/2),sqrt(1+e));

M = deg2rad(E) - e*sind(E);


end
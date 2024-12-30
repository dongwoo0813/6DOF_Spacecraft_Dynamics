function [TA] = MeanAnomalytoTrueAnomaly(M, e)

%
% Converts Mean Anomaly to True Anomaly using Newton's Method
%
%   INPUT: M -- deg, Mean Anomaly of Elliptical or Circular Orbit
%          e -- nd, eccentricity
%
%   OUTPUT: TA -- deg, True Anomaly of Elliptical or Circular Orbit
%
%
%


eps = 1e-10;

E_new = 0;

% Converts Mean Anomaly from degrees to radians.
M = deg2rad(M);

if M < pi
    E = M + e/2;
else
    E = M - e/2;
end


f_E = E - e*sin(E) - M;

f_prime_E = 1 - e*cos(E);


while abs(E_new - E) > eps
    E = E_new;

    f_E = E - e*sin(E) - M;
    f_prime_E = 1 - e*cos(E);

    E_new = E - f_E/f_prime_E;
end

% Converts Eccentric Anomaly from radians to degrees.
E = rad2deg(E);

TA = 2*atan2d(sqrt(1+e)*tand(E/2),sqrt(1-e));


if TA < 0
    TA = 360 + TA;
elseif TA > 360
    TA = TA - 360;
end


end
function [r_ijk, v_ijk]= COEstoRV(COEs, mu)

% Classical Orbital Elements to position and velocity vectors in 
% Cartesian coordinates.
%
%
% Input: COEs = [a, e, i, OMEGA, omega, theta]
%
%           a = km, semi-major axis
%           e = nd, eccentricity
%           i = deg, inclination
%           RAAN (OMEGA) = deg, RAAN
%           AoP (omega) = deg, Argument of Periapsis
%           MA (theta) = deg, Mean Anomaly
%
%
%   SINGULARITY NOT CONSIDERED IN THIS FUNCTION YET.
%
%
% 

format long g

a = COEs(1);
e = COEs(2);
i = COEs(3);
RAAN = COEs(4);
omega = COEs(5);
theta = COEs(6);


deg2rad = pi/180;

% Energy
E = -mu/(2*a);

p = a*(1 - e^2);


r_pqw = (p/(1 + e*cosd(theta)))*[cosd(theta); sind(theta); 0];

v_pqw = sqrt(mu/p)*[-sind(theta); e + cosd(theta); 0];

e1 = [1;0;0];
e3 = [0;0;1];

R3_RAAN = rotationMatrix(e3, -RAAN*deg2rad);

R1_i = rotationMatrix(e1, -i*deg2rad);

R3_AoP = rotationMatrix(e3, -omega*deg2rad);

Q_pqw2ijk = R3_RAAN*R1_i*R3_AoP;

r_ijk = Q_pqw2ijk*r_pqw;

v_ijk = Q_pqw2ijk*v_pqw;



end
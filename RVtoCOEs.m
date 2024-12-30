function [COEs]= RVtoCOEs(r_ijk, v_ijk)

% Cartesian coordinates 
%
% Input: r_ijk -- Position vector of the spacecraft in Cartesian coordinate
%        v_ijk -- Velocity vector of the spacecraft in Cartesian coordinate
%           mu -- Gravitational Parameter (km^3/s^2)  
%
% Output: COEs = [a, e, i, OMEGA, omega, theta]
%
%           a = km, semi-major axis
%           e = nd, eccentricity
%           i = deg, inclination
%           RAAN (OMEGA) = deg, RAAN
%           AoP (omega) = deg, Argument of Periapsis
%           TA (theta) = deg, True Anomaly
%
%    Function from Curtis, H., Orbital Mechanics for Engineering Students
%    Appendix D8

global mu

format long g

eps = 1e-10;

norm_r = norm(r_ijk, 2);

norm_v = norm(v_ijk, 2);

E = norm_v^2/2 - mu/norm_r;

a = -mu/(2*E);

% Eccentricity Vector
e_vec = 1/mu*( (norm_v^2 - mu/norm_r)*r_ijk - dot(r_ijk, v_ijk)*v_ijk);


% Eccentricity
if norm(e_vec, 2) < eps
    e = 0;
else
    e = norm(e_vec, 2);
end

% Angular Momentum Vector
h_vec = cross(r_ijk, v_ijk);

% Angular Momentum
h = norm(h_vec, 2);

if abs(acosd(h_vec(3)/h)) < eps
    i = 0;
else
    i = acosd(h_vec(3)/h);
end

i_ijk = [1; 0; 0];

k_ijk = [0; 0; 1];

n_vec = cross(k_ijk, h_vec);

n = norm(n_vec,2);


%...Equation 4.13a (incorporating the case e = 0):
%   True Anomaly
if e > eps
    TA = acosd(dot(r_ijk, e_vec)/(norm_r*e))
    if dot(r_ijk, v_ijk) < 0
        TA = 360 - TA;
    end
else

if (dot(n_vec, r_ijk)/(norm_r*n) - 1 ) < eps
    cosTA = 1;
else
    cosTA = cosTA;
end

    cp = cross(n_vec, r_ijk)

    if cp(3) >= 0
        TA = acosd(cosTA);
    else
        TA = 360 - acos(cosTA);
    end
end


% Right Ascension
if n ~= 0
    RAAN = acosd(dot(n_vec, i_ijk)/n);
    if n_vec(2) < 0
        RAAN = 360 - RAAN;
    end
else
    RAAN = 0;
end


% Equation 4.12 in Orbital Mechanics for Engineering Students        
% Argument of Periapsis
if n ~= 0
    if e > eps
        w = acosd(dot(n_vec, e_vec)/n/e);
        if e_vec(3) < 0
            w = 360 - w;
        end
    else
        w = 0;
    end
else
    w = 0;
end





COEs = [a, e, i, RAAN, w, TA];

end
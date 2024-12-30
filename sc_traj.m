function [XDOT]= sc_traj(t, X, mu)

%
% Spacecraft Orbit Dynamics Equation
%
% Input:
%
%
% Output:
%
%
%

global mu

format long g


r = X(1:3);
v = X(4:6);


norm_r = norm(r,2);


rdot = v;
vdot = -mu/norm_r^3*r;


XDOT = [rdot; vdot];
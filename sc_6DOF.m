function [XDOT]= sc_6DOF(t, X, mu)

%
% Spacecraft Orbit Dynamics Equation
%
% Input: t -- sec or TU, time
%        X -- State [r, v, omega]
%
% Output: XDOT -- Derivative of state
%
%
%
%
%
% Inertial Frame
%
% omega -- angular velocity of the body frame with respect to inertial
% frame, expressed in inertial frame

global mu J dt

format long g


r = X(1:3);
v = X(4:6);
omega = X(7:9);

norm_r = norm(r,2);


% beta = omega*dt;



% Control Input
u = zeros(3,1);

rdot = v;
vdot = -mu/norm_r^3*r;
omegadot = (J\eye(3))*(-cross(omega, J*omega) + u); % Assume no control input


XDOT = [rdot; vdot; omegadot];
clc
clear all

global mu J dt

mu = 3.986004418e14/(1000^3);

J = diag([3.05; 2.66; 2.18]);

% I1 = 3;
% I2 = 4;
% I3 = 2;
% J = [I1 0 0; 0 I2 0; 0 0 I3];
% Iinv  = inv(J);
% 
% r0 = [6678; 0; 0];
% 
% v0 = sqrt(mu/norm(r0))*[0;1;0]; %km/s

% perigee altitude
pa = 6378 + 813.2; %km
e = 0.7;
a = pa/(1 - e);

i = 63.4;
RAAN = 329.6;
AoP = 270;
TA = 180;

COEs = [a; e; i; RAAN; AoP; TA]

[r0, v0] = COEstoRV(COEs, mu)


% COEs = RVtoCOEs(r0, v0)

% e_vec = 1/mu*( (norm(v0)^2 - mu/norm(r0))*r0 - dot(r0, v0)*v0)



% omega0 = [0; 0; 0]; %rad/s

omega0 = [1; 2; 3]*1e-3; %rad/s

X0 = [r0;v0; omega0];

dt = 0.005;

time = 0:dt:5000;



options = odeset('RelTol', 1e-12, 'abstol', 1e-20);

[t, X] = ode45(@sc_6DOF, time, X0, options);

r_Mat = X(:,1:3);

x_Mat = X(:,1);
y_Mat = X(:,2);
z_Mat = X(:,3);

v_Mat = X(:,4:6);

omega_Mat= X(:,7:9);




beta = omega_Mat*dt;

nums = length(time);

STM = zeros(4,4,nums);

% q0 = [1; 0; 0; 0];

R = -pi/4;
P = pi/4;
Y = pi/3;

q0 = [cos(R/2)*cos(P/2)*cos(Y/2) - sin(R/2)*sin(P/2)*sin(Y/2);
         sin(R/2)*cos(P/2)*cos(Y/2) + cos(R/2)*sin(P/2)*sin(Y/2);
         cos(R/2)*sin(P/2)*cos(Y/2) - sin(R/2)*cos(P/2)*sin(Y/2);
         cos(R/2)*cos(P/2)*sin(Y/2) + sin(R/2)*sin(P/2)*cos(Y/2)];



q_Mat = zeros(4,nums + 1);

q_Mat(:,1) = q0;

for ii = 1:nums
    u = beta(ii, 1);
    v = beta(ii, 2);
    w = beta(ii, 3);

    rho = (1/2)*[0, -u, -v, -w;
                 u, 0, w, -v;
                 v, -w, 0, u;
                 w, v, -u, 0];


    if norm(beta(ii, :)) < 1e-12
        Var = 0;
    else
        Var = sin(norm(beta(ii, :),2))/norm(beta(ii, :), 2);
    end

    STM(:,:,ii) = cos(norm(beta(ii, :),2))*eye(4) + 2*Var*rho;
    q_Mat(:, ii + 1) = STM(:,:, ii)*q_Mat(:, ii);
end


fig = 0;

fig = fig + 1;


figure(fig)
plot(t, q_Mat(:,2:end));

legend('qw', 'qv1', 'qv2', 'qv3')



fig = fig + 1;

figure(fig)
plot(t, omega_Mat)


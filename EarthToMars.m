clear; clc;
load('Launch_Data.mat', "min_vi", "min_pi", "min_depart","min_arrive", "tof", "min_C3");
addpath 'C:\Users\abina\Downloads\mice\mice';
addpath 'C:\Users\abina\Downloads\mice\mice\src\mice';
addpath 'C:\Users\abina\Downloads\mice\mice\lib';
cspice_furnsh('naif0012.tls');
cspice_furnsh('pck00010.tpc'); 
cspice_furnsh('de440.bsp'); 
days = 24*60*60;
%% Initial Conditions
x0 = min_pi;
v0 = min_vi;
%% Interplanetary Transfer
depart = min_depart;
arrive = min_arrive;
tof = arrive - depart;
t = [0, tof];
sunMu = 1.327 * 10^11; % km^3/s^2
options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
[t,sol1] = ode45(@orbdyn, t, [x0;v0], options, sunMu);
tstep = t(2)- t(1);
tmax = t(end);
plot3(sol1(:, 1), sol1(:, 2), sol1(:, 3), 'MarkerSize',10);
hold on;
%% Arrival at Mars
% MarsPos= cspice_spkgeo(4, min_arrive, 'J2000', 0);
% state0 = sol1(end, 1:6) - MarsPos';
% marsMu = 4.2828 * 10^4;
% t = [0, tof];
% options = odeset('RelTol',1e-12, 'AbsTol',1e-12);
% [t,sol2] = ode45(@orbdyn, t, state0, options, marsMu);
% plot3(sol2(:, 1), sol2(:, 2), sol2(:, 3), 'MarkerSize',10);
%% Orbit Plot

% Finding the State Vectors for Earth's Orbit for a Year starting departure date.
[EarthVec, ~] = cspice_spkezr('EARTH BARYCENTER', [min_depart:tstep: min_arrive], 'J2000', 'NONE', 'SUN');
% Finding State Vectors for Mars' orbit for 2 years since departure date.
[MarsVec, ~] = cspice_spkezr('MARS BARYCENTER', [min_depart:days: min_depart + (2*365*days)], 'J2000', 'NONE', 'SUN');

% Finding the Final Position of Mars at Arrival.
[MarsPos, ~] = cspice_spkezr('MARS BARYCENTER', min_arrive, 'J2000', 'NONE', 'SUN');

figure(1);
% Plotting Tajectory of the Transfer.
plot3(sol1(:, 1), sol1(:, 2), sol1(:, 3), 'MarkerSize',10);
hold on;
% Plotting Earth's Orbit.
plot3(EarthVec(1,:), EarthVec(2,:), EarthVec(3,:), 'MarkerSize',10);
hold on;
% Plotting Mars' orbit.
plot3(MarsVec(1,:), MarsVec(2,:), MarsVec(3,:), 'MarkerSize',10);
axis equal;
title('Minimum C3 trajectory from Earth to Mars');
xlabel('x-axis (km)');
ylabel('y-axis (km)');
zlabel('z-axis (km)');

grid on;
hold on
% Plotting Sun as sphere at the center.
[X,Y,Z] = sphere;
r = 695700*100;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;

% Plotting the Earth as sphere at the departure location.
surf(X2,Y2,Z2);
hold on;
[X,Y,Z] = sphere;
r = 6378*1000;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2+EarthVec(1,1),Y2+EarthVec(2,1),Z2+EarthVec(3,1));

hold on

% Plotting Mars as sphere at the arrival location.
r = 3396.2*1000;
X2 = X * r;
Y2 = Y * r;
Z2 = Z * r;
surf(X2+MarsPos(1),Y2+MarsPos(2),Z2+MarsPos(3));
legend('Min C3 Trajectory', 'Earth Orbit', 'Mars Orbit', 'Sun', 'Earth', 'Mars');
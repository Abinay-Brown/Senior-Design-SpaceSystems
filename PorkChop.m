clear; clc;
%% Load Mice Kernels
addpath 'C:\Users\abina\Downloads\mice\mice';
addpath 'C:\Users\abina\Downloads\mice\mice\src\mice';
addpath 'C:\Users\abina\Downloads\mice\mice\lib';
cspice_furnsh('naif0012.tls');
cspice_furnsh('pck00010.tpc'); 
cspice_furnsh('de440.bsp'); 

%% Looping through different departure and arrival date combinations
days = 24*60*60; % # of seconds in a day
mu = 1.327 * 10^11; % Sun's Gravitational Parameter

depart = cspice_str2et('Jan 01, 2027 00:00:00.0');

t0_arrive = cspice_str2et('Jan 01, 2028 00:00:0.0');
t1_arrive = cspice_str2et('Dec 31, 2028 23:59:59');

st = depart:days: depart + (1*365*days);
fn = t0_arrive:days: t1_arrive;

% Creating a meshgrid of vectors of starting and ending dates
[X, Y] = meshgrid(st, fn);

% Preallocating zero vectors to store the Magnitudes of C3_Calc 
% and Vinf arrival data for the different combinations.
C3_Calc = zeros(length(fn), length(st));
Vinf_Calc = zeros(length(fn), length(st));


i = 0; j = 0; % Looping variables for saving C3 and Vinf data.
for start = st
    i = i + 1;
    for finish = fn
        j = j+1;
        % Compute Earth and Mars State Vector for the given departure
        % and arrival Combination.
        [state1, ~] = cspice_spkezr('EARTH BARYCENTER', start, 'J2000', 'NONE', 'SUN');
        [state2, ~] = cspice_spkezr('MARS BARYCENTER', finish, 'J2000', 'NONE', 'SUN');
        % Time of flight
        tof = (finish - start);
        % Determine Heliocentric Vel at departure and arrival for the
        % trajectory using lamberts problem.
        [vi, vf] = glambert(mu, state1, state2, tof, 0);
        % Calculate C3.
        c3 = norm(vi - state1(4:6))^2;
        % Calculate Vinf
        vinf = norm(vf - state2(4:6));
        % Store C3 and Vinf.
        C3_Calc(j, i) = c3;
        Vinf_Calc(j, i ) = vinf;
    end
    j = 0;
end

%% Specific C3

t_arrive = cspice_str2et('Aug 17, 2028 0:0:0.0');
[~, arrive_idx]= min(abs(fn - t_arrive));

result = zeros(1, length(st));

for i = 1:length(st)
    [state1, ~] = cspice_spkezr('EARTH BARYCENTER', st(i), 'J2000', 'NONE', 'SUN');
    [state2, ~] = cspice_spkezr('MARS BARYCENTER', t_arrive, 'J2000', 'NONE', 'SUN');
    tof = (t_arrive - st(i));
    [vi, vf] = glambert(mu, state1, state2, tof, 0);
    c3_spec = norm(vi - state1(4:6))^2;
    result(i) = c3_spec;
end
[min_C3, depart_idx] = min(result);
min_depart = st(depart_idx);
min_arrive = t_arrive;
tof = (min_arrive - min_depart);
[state1, ~] = cspice_spkezr('EARTH BARYCENTER', min_depart, 'J2000', 'NONE', 'SUN');
[state2, ~] = cspice_spkezr('MARS BARYCENTER', min_arrive, 'J2000', 'NONE', 'SUN');
[vi, vf] = glambert(mu, state1, state2, tof, 0);
min_vi = vi;
min_pi = state1(1:3);
disp(['Min Depart: ', cspice_et2utc(min_depart, 'C', 6)]);
disp(['Min Arrival: ', cspice_et2utc(min_arrive, 'C', 6)]);
sent = sprintf('Min C3 (%0.3f km^2/s^2)', min_C3);
disp(sent);

save('Launch_Data.mat', "min_vi", "min_pi", "min_depart","min_arrive","tof", "min_C3");
%% Plotting
% Vectors for the contour plots.
C3 = 0:5:30;
Vinf = 0:1:10;
% Plot 1 is Porkchop of C3 magnitudes at different departure and arrival
% date combinations.
subplot(1,2,1);
contour(X,Y,C3_Calc, C3, 'ShowText','on');
hold on;
% Plotting the Min C3 Values.
scatter(st(depart_idx), fn(arrive_idx), 'r*');
grid on;
sent = sprintf('Min C3 (%0.3f km^2/s^2)', min_C3);
legend('C3 Contour', sent);
title('PorkChop Plot For Launch Window 1st Jan 2026 to 1st Jan 2027');
xlabel("Launch Date (Days Since Jan 1st 2027)");
ylabel("Arrival Date (Days Since Jan 1st 2028)");
xticks(depart:days*91: depart + (1*365*days));
xticklabels({'0', '91.25', '182.5', '273.75', '365'});
yticks(t0_arrive:days*91: t1_arrive);
yticklabels({'0', '91.25', '182.5', '273.75', '365'});
title('Launch C3 Porkchop Plot');
colormap('jet')
colorbar;

subplot(1,2,2);
% Plot 2 Porkchop of Vinf arrival magnitudes at different
% departure and arrival combinations.
contour(X,Y,Vinf_Calc, Vinf, 'ShowText','on');

hold on;
% Plotting the Min C3 Values.
scatter(st(depart_idx), fn(arrive_idx), 'r*');
sent = sprintf('Selected Arrival V_I_n_f (%0.3f km/s)', Vinf_Calc(arrive_idx, depart_idx));
legend('Arrival V_I_n_f Contour', sent);
xlabel("Launch Date (Days Since Jan 1st 2027)");
ylabel("Arrival Date (Days Since Jan 1st 2028)");
xticks(depart:days*91: depart + (1*365*days));
xticklabels({'0', '91.25', '182.5', '273.75', '365'});
yticks(t0_arrive:days*91: t1_arrive);
yticklabels({'0', '91.25', '182.5', '273.75', '365'});
title('Arrival V_i_n_f Porkchop Plot');
colormap('jet')
colorbar;
grid on;


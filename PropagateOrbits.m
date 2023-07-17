% Format workspace
clear all;
close all;
clc;
format longG;
delete('*~');

%% Script Constants
% Constants
mu = 398600.4418;           % Gravitational parameter for Earth
RE = 6378.137;              % Earth's radius (km)
omega = [0; 0; 7.2921e-5;];     % The Earth's rotation when seen from the ECI frame (rad/s)

% Orbit Variables (COEs)
a = 26562;                  % Semi-major axis (km)
e = 0.74105;                % Eccentricity
M_0 = (32) * pi / 180;      % Inital Mean Anomaly
i = 63.4 * pi / 180;        % Inclination
RAAN = 260 * pi / 180;      % Right Ascension
argPeri = 270 * pi / 180;   % Argument of Perigee 

%% Calculating useful quantities
% Find initial eccentric anomaly, EA, and true anomaly, TA
EA = Kepler(e, M_0, 1e-10);
TA = 2 * atan2( sqrt((1+e)/(1-e)), 1 / (tan(EA / 2)));
%fprintf("True Anomaly, Theta:    %f         (rad)\n", TA);

% Calculate orbital period, T
n = sqrt( mu / (a^3) );     % Mean Motion
T = 2 * pi / n;
%fprintf("Orbit Period, T:        %.4f       (s)\n", T);

%% Propagate Orbits (Analytical Solution)
% Setup
timeResolution = 1000;
timeVector = 0:T/timeResolution:T;
MA = M_0 + n*(timeVector);

% Loop
for nn = 1:length(timeVector)

    % Find EA
    EA(nn) = Kepler(e, MA(nn), 1e-10);

    % Find TA
    TA(nn) =  2 * atan2( sqrt((1+e)/(1-e)), 1 / (tan(EA(nn) / 2)));
    %TA(nn) =  2 * atan2( sqrt(1+e) * tan(EA(nn) / 2), sqrt(1-e));

    coe = [a, e, i, RAAN, argPeri, TA(nn)];
    ECI(:,nn) = COE2RV(coe, mu);
end

%% Propagate Orbits (ODE solution, ECI)
% Setting options and initial conditions.
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
initialState = ECI(:, 1);

% Solving
[t, Y] = ode45(@(t, X) TBP_ECI(t, X, mu), timeVector, initialState, options);
ECISol = transpose(Y); % Transpose to be in same format as ECI matrix.

%% Propagate Orbits (ODE solution, ECEF)
% Setting options and initial conditions.
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
initialState = ECI(:,1);
initialState(4:6) = initialState(4:6) - cross(omega,initialState(1:3));
timeVectorECEF = 0:T/timeResolution:2*T; % New timeVector

% Solving
[t, Z] = ode45(@(t, X) TBP_ECEF(t, X, mu), timeVectorECEF, initialState, options);
ECEF = transpose(Z); % Transpose to be in same format as ECI matrix.

%% Find Error between Analytical and ODE solutions
% Position error
numPos = ECI(1:3, :);
odePos = ECISol(1:3, :);
posError = numPos - odePos;
posError = vecnorm(posError);

% Velocity error
numVel = ECI(4:6, :);
odeVel = ECISol(4:6, :);
velError = numVel - odeVel;
velError = vecnorm(velError);

%% Finding specific angular momentum over time, for TBP
h = cross(ECISol(1:3, :), ECISol(4:6, :));
h = vecnorm(h);
avgH = mean(h);

%% Plot
%% Graph Anomalies
%figure('units','normalized','outerposition',[0 0 1 1]); % Full-screen figure.
%subplot(2, 3, 1);
%hold on;
%grid on;
%plot(timeVector, MA);
%plot(timeVector, EA);
%plot(timeVector, TA);
%yticks(0:pi/2:2*pi);
%yticklabels({'0','\pi/2','2\pi','3\pi/2','2\pi'});
%legend('Mean Anomaly', 'Eccentric Anomaly', 'True Anomaly', 'Location', 'northwest')
%title("Anomalies over Time");
%xlabel("Time (s)");
%ylabel("Anomaly (rad)");

% Plot Earth for reference
%subplot(2, 3, 2);
figure();
[Xe,Ye,Ze] = sphere(50);
surf(RE*Xe, RE*Ye, RE*Ze, 'EdgeColor', 'none', 'FaceColor', '#235761');
hold on;
grid on;
axis equal;
xlabel('X, ECI (km)');
ylabel('Y, ECI (km)');
zlabel('Z, ECI (km)');

% Plot orbit in 3D (from analytical solution)
%plot3(ECI(1,:), ECI(2,:), ECI(3,:), 'Color', 'k', 'LineWidth', 2);
%plot3(ECI(1,1), ECI(2,1), ECI(3,1), 'ok', 'MarkerFaceColor', 'b'); % Start
%plot3(ECI(1,end), ECI(2,end), ECI(3,end), 'ok', 'MarkerFaceColor', 'r'); %End
% Printing first and last ECI state (from analytical solution)
fprintf("--ECI--\n");
fprintf("First Position (X, Y ,Z): %.4f, %.4f, %.4f\nFirst Velocity (X, Y, Z): %.4f, %.4f, %.4f\n", ECI(:, 1));
fprintf("Last Position (X, Y ,Z):  %.4f, %.4f, %.4f\nLast Velocity (X, Y, Z):  %.4f, %.4f, %.4f\n", ECI(:, end));
% Plot orbit in 3D (from ODE solution)
plot3(ECISol(1,:), ECISol(2,:), ECISol(3,:),'Color', 'Magenta', 'LineStyle', '--', 'LineWidth', 2);
title("Orbital Solutions")
legend('Earth', 'Analytical Solution', 'Initial Position', 'Final Position', 'ODE Solution');
view(cross(ECI(1:3, 1), ECI(1:3, timeResolution/4)));

% Plot ECEF orbit:
%subplot(2, 3, 3)
figure();
[Xe,Ye,Ze] = sphere(50);
surf(RE*Xe, RE*Ye, RE*Ze, 'EdgeColor', 'none', 'FaceColor', '#235761');
hold on;
grid on;
axis equal;
xlabel('X, ECEF (km)');
ylabel('Y, ECEF (km)');
zlabel('Z, ECEF (km)');
plot3(ECEF(1,:), ECEF(2,:), ECEF(3,:), 'Color', 'k', 'LineWidth', 2);
plot3(ECEF(1,1), ECEF(2,1), ECEF(3,1), 'ok', 'MarkerFaceColor', 'b'); % Start
plot3(ECEF(1,end), ECEF(2,end), ECEF(3,end), 'ok', 'MarkerFaceColor', 'r'); %End
legend('Earth', 'Orbital Path', 'Initial Position', 'Final Position');
title("Orbital Solution over Two Orbits (ECEF)")
% Printing first and last ECEF state
fprintf("--ECEF--\n");
fprintf("First Position (X, Y ,Z): %.4f, %.4f, %.4f\nFirst Velocity (X, Y, Z): %.4f, %.4f, %.4f\n", ECEF(:, 1));
fprintf("Last Position (X, Y ,Z):  %.4f, %.4f, %.4f\nLast Velocity (X, Y, Z):  %.4f, %.4f, %.4f\n", ECEF(:, end));

%% Plot Analytical vs ODE error
%subplot(2, 3, 4);
%hold on;
%grid on;
%yyaxis left;
%plot(timeVector, posError);
%ylabel("Position Error (kilometers)");
%yyaxis right;
%plot(timeVector, velError);
%ylabel("Velocity Error (kilometers per second)");
%title("Position and Velocity Differences between Analytical and ODE Solutions over Time");
%xlabel("Time (seconds)");

%% Plot specific angular momentum over time
%subplot(2, 3, 5);
%hold on;
%grid on;
%plot(timeVector, h, 'k', 'LineWidth', 2);
%xlabel("Time (seconds)");
%title("Specific Angular Momentum over Time");
%ylabel("Specific Angular Momentum, h (square kilometers per second)");
%ylim([0, h(1) .* 1.5]);
%approximation = sqrt(mu.*a.*(1 - e^2));
%yline(approximation, 'Color', 'cyan', 'LineStyle', '--', 'LineWidth', 2);
%leg1 = legend('Specific Angular Momentum', '$\sqrt{\mu a (1 - e^2)}$');
%set(leg1,'Interpreter','latex');





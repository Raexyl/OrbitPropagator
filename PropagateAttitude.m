% Format workspace
clear all;
close all;
clc;
format longG;
delete('*~');

%% Script Variables
theta = [    74.0792 * pi/180;  % Radians
             28.7436 * pi/180;
            -56.0137 * pi/180;
];
omega = [   -4.3198e-5;         % Radians / second
             9.2422e-5;
             0.0001;
];
I = [       2500, 0, 0;         % kg m^2
            0, 5000, 0;
            0, 0, 6500;
];
timeVector = 0:10:3600;         % seconds

%% Propagate Attitude (ODE solution)
% Setting options and initial conditions.
options = odeset('RelTol', 1e-12, 'AbsTol', 1e-12);
initialState = [theta; omega;];

% Solving
[t, Y] = ode45(@(t, X) AttitudeDynamics(t, X, I), timeVector, initialState, options);
sol = transpose(Y); % Transpose for convenience.

%% Plotting
figure('units','normalized','outerposition',[0 0 1 1])
subplot(2, 2, 1);
plot(t, sol(1, :));
hold on;
plot(t, sol(2, :));
plot(t, sol(3, :));
title("Euler Angles over Time");
legend("Roll", "Pitch", "Yaw");
xlabel("Time (s)");
ylabel("Angle (rad)");

subplot(2, 2, 2);
plot(t, sol(4, :));
hold on;
plot(t, sol(5, :));
plot(t, sol(6, :));
title("Rotational Velocity Components over Time");
legend("Roll", "Pitch", "Yaw");
xlabel("Time (s)");
ylabel("Angular Speed (rad / s)");

%% More Calculations
% Angular Momentum:
am = I * sol(4:6, :);
for nn = 1:length(am)
    nam(nn) = norm(am(:, nn));
end

% Rotational Kinetic Energy:
rke =  dot(sol(4:6, :), am) .* 0.5;

%% More Plots
subplot(2, 2, 3)
plot(t, nam);
title("Specific Angular Momentum over Time");
xlabel("Time (s)");
ylabel("Specific Angular Momentum (m^2 / s)");
ylim([0 max(nam) * 1.5]);

subplot(2, 2, 4)
plot(t, rke);
title("Rotational Kinetic Energy over Time");
xlabel("Time (s)");
ylabel("Rotational Kinetic Energy (joules)");
ylim([0 max(rke) * 1.5]);


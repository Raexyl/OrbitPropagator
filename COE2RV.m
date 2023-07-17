function [X] = COE2RV(coe, mu)
%COE2RV Converts from classical Keplerian orbital elements to the ECI frame of reference.
%
%   This function is of the form COE2RV(coe, mu)
%   where mu is the gravitational parameter and 
%   coe is an input vector of classical orbital elements as shown:
%       [semi-major axis, eccentricity, inclincation, right-ascension, argument of periapsis, true anomaly]
%
%   A vector [X] is returned, which first contains the 
%   3D position vector, then the 3D velocity vector:
%       [x, y, z, xdot, ydot, zdot]

% Re-naming coe elements for convenience
a = coe(1);
e = coe(2);
i = coe(3);
RAAN = coe(4);
argPeri = coe(5);
TA = coe(6);

%% Finding transformation matrix IP (perifocal to inertial ECI)
MatNI = [
     cos(-RAAN), sin(-RAAN),  0;
    -sin(-RAAN), cos(-RAAN),  0;
    0,          0,          1;
];
MatHN = [
        1,         0,          0;
        0,          cos(-i),     sin(-i);
        0,          -sin(-i),    cos(-i);
];
MatPH = [
     cos(-argPeri), sin(-argPeri),    0;
    -sin(-argPeri), cos(-argPeri),    0;
    0,          0,                  1;
];

MatrixIP = MatNI * MatHN * MatPH;
%MatrixIP = transpose(MatrixIP);

%% Finding r in the perifocal frame
% Using the orbit equation:
r_min = a.*(1-e);
P = r_min .* (1 + e);   % h^2 / mu
h2 = P * mu;            % h^2
h = sqrt(h2);           % h
mu_h = mu / h;          % mu / h
r = P / (1 + e.*cos(TA));

% Calculating position
pos = MatrixIP * [r.*cos(TA); r.*sin(TA); 0;];

% Calculating velocity
vel = MatrixIP * [-mu_h.*sin(TA); mu_h.*(cos(TA) + e); 0;];

[X] = [pos(1); pos(2); pos(3); vel(1); vel(2); vel(3);];

end


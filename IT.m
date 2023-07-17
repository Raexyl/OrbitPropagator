%% Calculating [IT]
% Required variables:
r = [7578.37; 15678.86; 9466.82;]; % km
r = r./norm(r);
v = [-2.3580; -0.7963; -4.3611;];  % km / s
v = v./norm(v);
h = cross(r, v);
h = h./norm(h);

firstAxis = cross(v, h);

DCM = [
    firstAxis(1), v(1), h(1);
    firstAxis(2), v(2), h(2);
    firstAxis(3), v(3), h(3);
]

%% Calculating [BI]
% Calculating [TI]
DCM = transpose(DCM);   % Equal to inverse, as is an orthogonal matrix.

% Calculating [BT] from (1-2-1)(30, 20, 10) rotation sequence
a1 = 30 * pi / 180; % Angle 1
a2 = 20 * pi / 180; % Angle 2
a3 = 10 * pi / 180; % Angle 3

R1 = [
    1,      0,          0;
    0,      cos(a1),    sin(a1);
    0,     -sin(a1),    cos(a1);
];
R2 = [
    cos(a2),    0,     -sin(a2); 
    0,          1,      0;
    sin(a2),    0,      cos(a2);
];
R3 = [
    1,      0,          0;
    0,      cos(a3),    sin(a3);
    0,     -sin(a3),    cos(a3);
];

BT = R3 * R2 * R1;
format short
BI = BT * DCM   

%% Calculating Principal Axes and Angle
phi = acos((trace(BI) - 1) / 2 );    % Angle
doubleSinPhi = sin(phi) .* 2;
e1 = (BI(2, 3) - BI(3, 2)) / doubleSinPhi;
e2 = (BI(3, 1) - BI(1, 3)) / doubleSinPhi;
e3 = (BI(1, 2) - BI(2, 1)) / doubleSinPhi;

PrincipalAxes = [phi; e1; e2; e3;]

%% Calculating Equivalent Quaternion of [BI]
Q = [
    cos(phi / 2);
    e1 .* sin(phi / 2);
    e2 .* sin(phi / 2);
    e3 .* sin(phi / 2);
]
Qnorm = norm(Q)

%% Calculating Equivalent Euler Angles of [BI]
theta2 = asin(-BI(1,3));
theta1 = asin(BI(1, 2) / cos(theta2));
theta3 = asin(BI(2, 3) / sin(theta2));

% Euler Angles in degrees
EulerAngles = [theta1 * 180 / pi, theta2 * 180 /pi, theta3 * 180 / pi]

%% Calculating Inverse Matrix for Kinematic Differential Equation
mat = inv(BI)
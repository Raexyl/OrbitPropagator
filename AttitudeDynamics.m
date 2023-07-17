function [dXdT] = AttitudeDynamics(t, X, I)
%AttitudeDynamics A representation of the time-evolution of satellite attitude, using euler angles in the principal axes frame of reference.
%   Two first order ODEs re-create a second-order equation of motion.
%
%   X should be of the cartesian form [posX, posY, posZ, velX, velY, velZ]
%   where the prefix "pos-" represents angular position values, "vel-"
%   represents angular velocity values.
%
%   I refers to the matrix of inertia of the satellite.

    % Separate position and velocity vectors.
    theta = X(1:3);
    omega = X(4:6);

    % Matrix B from the kinematic differential equation for 3-2-1 euler angles.
    B = [1,     0,               -sin(theta(2));
         0,     cos(theta(1)),    sin(theta(1)).*cos(theta(2));
         0,    -sin(theta(1)),    cos(theta(1)).*cos(theta(2));
    ];
    B = inv(B);
    
    % Omega tilde:
    omegatilde = [  0,       -omega(3), omega(2);
                    omega(3), 0,       -omega(1);
                   -omega(2), omega(1), 0;
        ];
    
    % Apply equations.
    thetadot = B * omega;
    omegadot = inv(I) * -(omegatilde * I * omega);

    % Return output of the form
    dXdT = [thetadot; omegadot;];
end


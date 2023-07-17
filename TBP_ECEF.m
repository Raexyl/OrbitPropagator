function [dXdT] = TBP_ECEF(t, X, mu)
%TBP_ECEF A representation of the two-body problem in the ECEF reference frame.
%   Two first order ODEs re-create the second-order equation of motion.
%
%   X should be of the cartesian form [posX, posY, posZ, velX, velY, velZ]
%   where the prefix "pos-" represents position values, and "vel-"
%   represents velocity values.
%
%   Mu refers to the gravitational parameter, equal to GM.

    omega = [0; 0; 7.2921e-5;];     % The Earth's rotation when seen from the ECI frame (rad/s)


    % Separate position and velocity vectors. X.
    r = X(1:3);
    rdot = X(4:6);

    % Apply equations of motion
    rMag = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    rdotdot = - cross(omega .* 2, rdot) - cross(omega, cross(omega, r))  - r .* (mu / (rMag^3));

    % Return output of the form
    dXdT = [rdot; rdotdot;];
end


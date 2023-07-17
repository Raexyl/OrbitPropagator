function [dXdT] = TBP_ECI(t, X, mu)
%TBP_ECI A representation of the two-body problem in the ECI reference frame.
%   Two first order ODEs re-create the second-order equation of motion.
%
%   X should be of the cartesian form [posX, posY, posZ, velX, velY, velZ]
%   where the prefix "pos-" represents position values, "vel-"
%   represents velocity values.
%
%   Mu refers to the gravitational parameter, equal to GM.

    % Separate position and velocity vectors.
    r = X(1:3);
    rdot = X(4:6);
    
    % Calculate acceleration
    rMag = sqrt(r(1)^2 + r(2)^2 + r(3)^2);
    rdotdot =  r .* (-mu / (rMag^3));

    % Return output of the form
    dXdT = [rdot; rdotdot;];
end


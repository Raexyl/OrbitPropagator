function [E] = Kepler(e, M, tol)
%KEPLER This function takes an (e)ccentricity, an initial guess for M(ean Anomaly) and an error tolerance.
% It will use the Netwon-Raphson method to solve
% Kepler's Equation (M = E - e sin(E)) and plot the solution.

    %% Solve Equation

    % Initalising variables
    err = tol * 10; % initial error
    E = M;          % initial guess for solution
    n = 0;          % to store number of iterations

    % Iterating on that guess:
    while err > tol
        f = E - e.*sin(E) - M;
        df = 1 - (e.*cos(E));

        En = E - (f/df);

        err = En - e.*sin(En) - M;
        E = En;
        n = n + 1;
    end
    %fprintf("Newton-Raphson Iterations: %d\n", n);
    %fprintf("Eccentric Anomaly, E:   %f         (rad)\n", E);
end
function X = refineProjector(X)
% Refines a projector using a Newton-Schulz iteration
%
% Args:
%   X (double(\*,\*)): Matrix with $|| X^2-X || \le 1/4$
%
% Returns:
%   double(\*,\*): Refined projector
    X = (-2*X + 3*eye(size(X)))*X^2;
end

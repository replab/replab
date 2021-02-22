function [P, omega] = biorthoStepP(I, P)
% Biorthogonalization step for the given projection map
%
% The biorthogonality condition is ``P * I == eye(d)``.
%
% This function first scales ``P`` so that ``trace(P*I) == d``. It then computes
% a measure of orthogonality ``omega = norm(eye(d) - P*I, 'fro')``.
%
% If ``omega < 0.9``, it performs a orthogonalization step (Newton iteration); otherwise
% it applies the proper correction.
%
% Args:
%   I (double(D,d)): Injection map
%   P (double(d,D)): Projection map
%
% Returns
% -------
%   P: double(d,D)
%     Approximately biorthogonal projection map
%   omega: double
%     Orthogonality measure after rescaling but before other corrections
    d = size(I, 2);
    PI = P*I;
    f = trace(PI)/d;
    P = P/f;
    PI = PI/f;
    Delta = eye(d) - PI;
    omega = norm(Delta, 'fro');
    % the spectral radius of PI is bounded by its Frobenius norm
    % https://math.stackexchange.com/questions/2965360/frobenius-norm-inequality-spectral-radius-is-smaller-than-frobenius-norm
    if omega < 0.9
        P = P + Delta*P;
    else
        replab.msg(2, 'Newton biortho, deviation=%6.2E => taking full inversion step', omega);
        P = PI\P;
        return
    end
end

function [P, ortho] = biorthoStepP(I, P)
% Biorthogonalizes the given injection/projection maps
%
% The biorthogonality condition is ``P * I == eye(d)``.
%
% Args:
%   I (double(D,d)): Injection map
%   P (double(d,D)): Projection map
%
% Returns
% -------
%   P: double(d,D)
%     Approximately biorthogonal projection map
%   ortho: double
%     Orthogonality measure before correction
    d = size(I, 2);
    %ortho = norm(P*I - eye(d), 'fro');
    %P = (P*I)\P;
    %return
    PI = P*I;
    f = trace(PI)/d;
    P = P/f;
    PI = PI/f;
    Delta = eye(d) - PI;
    ortho = norm(Delta, 'fro');
    % the spectral radius of PI is bounded by its Frobenius norm
    % https://math.stackexchange.com/questions/2965360/frobenius-norm-inequality-spectral-radius-is-smaller-than-frobenius-norm
    if ortho < 0.9
        P = P + Delta*P;
    else
        replab.msg(2, 'Newton biortho, deviation=%6.2E => taking full inversion step', ortho);
        P = PI\P;
        return
    end
end

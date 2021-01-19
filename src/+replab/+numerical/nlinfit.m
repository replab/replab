function [beta, R, J, covB, mse] = nlinfit(x, y, modelfun, beta0)
% Nonlinear least-squares fit
%
% Input/output parameters are mostly compatible with the ``nlinfit`` MATLAB
%
% Args:
%   x (double(n,\*)): Matrix of independent variables
%   y (double(n, 1)): Observed output
%   modelfun (function_handle): Function handle that computes ``y = modelfun(beta, x)``
%   beta0 (double(p, 1)): Initial model parameters
%
% Returns
% -------
%   beta: double(p, 1)
%     Estimated parameters
%   R: double(n, 1)
%     Residuals
%   J: double(n, p)
%     Estimated Jacobian
%   covB: double(p, p)
%     Estimated variance-covariance matrix
%   mse: double
%     Estimated mean squared error
    p = length(beta0);
    n = length(y);
    beta = replab.numerical.LMFnlsq2(@(b) modelfun(b, x) - y, beta0);
    yy = modelfun(beta, x);
    R = y - yy;
    % Compute Jacobian (again)
    J = zeros(n, p);
    for k = 1:p
        delta = eps^(1/3);
        beta1 = beta;
        beta1(k) = beta1(k) + delta;
        yy1 = modelfun(beta1, x);
        J(:,k) = (yy1 - yy)/delta;
    end
    mse = sum(abs(R).^2)/(n-p);
    covB = mse*inv(J'*J);
end

function J = finiteDifferenceJacobian(fun, r, x, fdss, typX)
% Numerical approximation of the Jacobian matrix
%
% Args:
%   fun (function_handle): Function to evaluate the Jacobian of
%   r (double(m, 1)): Value of the function at ``x``
%   x (double(n, 1)): Point at which to evaluate the function
%   fdss (double or double(n,1)): Scalar of vector step size factor
%   typX (double or double(n,1)): Typical x values
%
% Returns:
%   double(m,n): Evaluated Jacobian matrix
    r = r(:);
    m = length(r);
    n = length(x);
    J  = zeros(m, n);
    if isscalar(fdss)
        fdss = ones(n, 1) * fdss;
    end
    if isscalar(typX)
        typX = ones(n, 1) * typX;
    end

    % Compute finite step
    sg = sign(x);
    sg(sg == 0) = 1;
    delta = fdss.*sg.*max(abs(x), typX);

    for k = 1:n % for each parameter
        dx = delta(k);
        x1 = x;
        x1(k) = x1(k) + dx;
        r1 = fun(x1);
        J(:,k) = (r1(:)-r)/dx;
    end
end

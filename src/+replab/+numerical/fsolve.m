function [x, fval, exitFlag, J] = fsolve(fun, x0, varargin)
% Solve one or more of (over)determined nonlinear equations in the least squares sense
%
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/39564-lmfnlsq2
%
% Adapted some calling conventions according to the ``fsolve`` function of the Optimization Toolbox
%
% Not sure that the stopping criteria are exactly the ones advertised, a problem that may be present
% in the original LMFnlsq2 code; this needs to be checked if this function gets used outside RepLAB.
%
% A solution is obtained by a Fletcher's version of the Levenberg-Maquardt algorithm for minimization
% of a sum of squares of equation residuals.
%
% A main application is in curve fitting during processing of experimental or noisy data (see `.nlinfit`).
%
% Reference: Fletcher, R., (1971): A Modified Marquardt Subroutine for Nonlinear Least Squares. Rpt. AERE-R 6799, Harwell
%
% Original author: M. Balda, Institute of Thermomechanics, Academy of Sciences of The Czech Republic, balda AT cdm DOT cas DOT cz, miroslav AT balda DOT cz
%
% Args:
%   fun (function_handle): Function handle that evaluates a column vector of equation residuals.
%                          Residuals are defined as ``res = fun(x) - y``, where y is a column vector of length ``m``
%   x0 (double(n,1)): Vector of initial guesses of the unknown parameters, cannot be empty
%
% Keyword Args:
%   initDamping (double): Set the initial Levenberg-Marquardt parameter, default: 1e-2
%   finiteDifferenceStepSize (double or double(n,1)): Scalar or vector step size factor for finite differences
%   scaleProblem (logical): Whether to scale the problem using the Jacobian (sometimes improves the convergence if poorly scaled)
%   functionTolerance (double): Tolerance for the function value, default: 1e-6
%   maxIterations (integer): Maximum number of iterations, default: 500
%   jacobian (function_handle): Function to compute the Jacobian, returns a matrix of size ``m x n``
%   stepTolerance (double): Tolerance on the estimated parameters, default: 1e-6
%   typicalX (double(n,1)): Typical x values, default: ``ones(1,n)``
%
% Returns:
%   x: double(n,1)
%     Solution
%   fval: double(m,1)
%     Function value at solution point
%   exitFlag: integer
%     Flag corresponding to the termination criterion triggered. Solved (1), maximum iterations reached (0), error in evaluation (-1)

    n = length(x0);

    args = struct('initDamping', 1e-2, 'finiteDifferenceStepSize', sqrt(eps), ...
                  'scaleProblem', false, 'functionTolerance', 1e-6, 'maxIterations', 500, ...
                  'jacobian', [], 'stepTolerance', 1e-6, 'typicalX', ones(1, n));
    args = replab.util.populateStruct(args, varargin);

    assert(isa(fun, 'function_handle'), 'fun must be a function handle');
    x0 = x0(:);
    assert(isa(x0, 'double'));

    % Initialization

    t0 = clock;
    x = x0;
    maxit = args.maxIterations;
    epsf = args.functionTolerance;
    epsx = args.stepTolerance;
    if isscalar(epsx)
        epsx = epsx*ones(n,1);
    end
    jacobian = args.jacobian;

    l  = args.initDamping;
    lc = 1;
    r  = fun(x); % Initial residuals
    if any(isnan(r))
        exitFlag = -1;
        return
    end
    SS = r'*r;
    dx = zeros(n, 1);
    res = 1;
    cnt = 0;

    if isempty(jacobian)
        J = replab.numerical.finiteDifferenceJacobian(fun, r, x, args.finiteDifferenceStepSize, args.typicalX);
    else
        J = jacobian(x);
    end
    xJ = x;
    if any(isnan(J(:)))
        exitFlag = -1;
        return
    end

    A = J'*J;
    v = J'*r;

    if args.scaleProblem
        D = diag(A);
        D(D<=0) = 1;
        T = sqrt(D);
    else
        D = ones(n, 1);
        T = ones(n, 1);
    end

    Rlo = 0.25;
    Rhi = 0.75;

    % Solution

    while 1 % Main iteration cycle
        cnt = cnt + 1;
        d = diag(A);
        s = zeros(n,1);
        while 1 % Internal cycle (Denis: Broyden's update?)
            while 1 % Verify that A is positive definite?
                UA = triu(A,1);
                A = UA' + UA + diag(d+l*D);
                [U, p] = chol(A); %   Choleski decomposition
                if p == 0
                    break
                end
                l = 2*l;
                if l == 0
                    l = 1;
                end
            end
            dx = U\(U'\v); % Vector of x increments
            vw = dx'*v;
            fin = -1;
            if vw <= 0
                break % The END
            end

            for i = 1:n
                z = d(i)*dx(i);
                if i > 1
                    z = A(i,1:i-1)*dx(1:i-1)+z;
                end
                if i < n
                    z = A(i+1:n,i)'*dx(i+1:n)+z;
                end
                s(i) = 2*v(i)-z;
            end

            dq = s'*dx;
            s  = x - dx;
            rd = fun(s);
            if any(isnan(rd))
                exitFlag = -1;
                return
            end

            % ~~~~~~~~~~~~

            res = res + 1;
            SSP = rd'*rd;
            dS  = SS - SSP;
            fin = 1;
            if all(abs(dx) <= epsx) || res >= maxit || abs(dS) <= epsf
                break % The END
            end
            fin = 0;
            if dS >= Rlo*dq
                break
            end
            A = U;
            y = .5;
            z = 2*vw - dS;
            if z > 0
                y = vw/z;
            end
            if y > .5
                y = .5;
            end
            if y < .1
                y = .1;
            end
            if l == 0
                y = 2*y;
                for i = 1:n
                    A(i,i) = 1/A(i,i);
                end
                for i = 2:n
                    ii = i-1;
                    for j = 1:ii
                        A(j,i) = -A(j,j:ii)*A(j:ii,i).*A(i,i);
                    end
                end
                for i = 1:n
                    for j = i:n
                        A(i,j) = abs(A(i,j:n)*A(j,j:n)');
                    end
                end
                l = 0;
                tr = diag(A)'*D;
                for i = 1:n
                    z = A(1:i,i)'*T(1:i)+z;
                    if i < n
                        ii = i+1;
                        z = A(i,ii:n)*T(ii:n)+z;
                    end
                    z = z*T(i);
                    if z>l
                        l=z;
                    end
                end
                if tr<l
                    l=tr;
                end
                l = 1/l;
                lc = l;
            end
            l = l/y;
            if dS>0
                dS=-1e300;
                break
            end
        end % Internal cycle

        if fin
            break
        end
        if dS > Rhi*dq
            l = l/2;
            if l < lc
                l=0;
            end
        end
        SS = SSP;

        x = s;
        r = rd;
        if isempty(jacobian)
            J = replab.numerical.finiteDifferenceJacobian(fun, r, x, args.finiteDifferenceStepSize, args.typicalX);
        else
            J = jacobian(x);
        end
        xJ = x;
        if any(isnan(J(:)))
            exitFlag = -1;
            return
        end

        A = J'*J;
        v = J'*r;
    end % Main iteration

    if fin > 0
        if dS > 0
            x = s;
            r = rd;
        end
    end

    fval = r;

    if nargout > 3 && any(xJ ~= x)
        if isempty(jacobian)
            J = replab.numerical.finiteDifferenceJacobian(fun, r, x, args.finiteDifferenceStepSize, args.typicalX);
        else
            J = jacobian(x);
        end
    end

    if res >= maxit
        exitFlag = 0;
    else
        exitFlag = 1;
    end
end

classdef Tolerances < replab.Str
% Tolerance criteria used for the termination of iterative algorithms

    properties (SetAccess = protected)
        minIterations % (integer): Minimum number of iterations, default: 5
        maxIterations % (integer): Maximum number of iterations, default: 200
        windowSize % (integer): Number of steps to include in the monotonicity check, default: 8
        optimalityTol % (double): Optimality below which the iterations can stop, default: 1e-6
        relativeStepSizeTol % (double): Relative step size below which the iterations can stop, default: 1e-6
        comb % (integer(2,\*)): All pairs of indices in the window
    end

    methods

        function self = Tolerances(varargin)
            args = struct('minIterations', 5, 'maxIterations', 200, 'windowSize', 8, 'optimalityTol', 1e-6, 'relativeStepSizeTol', 1e-6);
            args = replab.util.populateStruct(args, varargin);
            comb = zeros(2, 0);
            for i = 1:args.windowSize
                for j = i+1:args.windowSize
                    comb(:, end+1) = [i;j];
                end
            end
            self.minIterations = args.minIterations;
            self.maxIterations = args.maxIterations;
            self.windowSize = args.windowSize;
            self.optimalityTol = args.optimalityTol;
            self.relativeStepSizeTol = args.relativeStepSizeTol;
            self.comb = comb;
        end

        function logHeader(self)
        % Logs a header at the verbosity level used during the iterations
            replab.msg(2, ' #iter   optim.   (slp)   stepSize (slp)');
            replab.msg(2, '-----------------------------------------');
        end

        function exitFlag = test(self, omega, delta, k)
        % Tests the termination criteria
        %
        % Args:
        %   omega (double(1,k)): Optimality values
        %   delta (double(1,k)): Relative step sizes
        %   k (integer): Current iteration number
            exitFlag = 0;
            R = self.windowSize;
            if k < R
                replab.msg(2, '%6d   %6.2E         %6.2E', k, omega(k), delta(k));
            elseif k >= R
                omegaWin = omega(k-R+1:k);
                deltaWin = delta(k-R+1:k);
                omegaMax = max(omegaWin);
                deltaMax = max(deltaWin);
                omegaLog = log10(max(omegaWin, 1e-200));
                deltaLog = log10(max(deltaWin, 1e-200));
                % Theil–Sen estimator
                iterDiff = self.comb(2,:) - self.comb(1,:);
                omegaDiff = omegaLog(self.comb(2,:)) - omegaLog(self.comb(1,:));
                omegaSlp = median(omegaDiff./iterDiff);
                deltaDiff = deltaLog(self.comb(2,:)) - deltaLog(self.comb(1,:));
                deltaSlp = median(deltaDiff./iterDiff);
                replab.msg(2, '%6d   %6.2E (%+1.1f)  %6.2E (%+1.1f)', k, omega(k), omegaSlp, delta(k), deltaSlp);
                if k < self.minIterations
                    return
                end
                if omegaMax <= self.optimalityTol && deltaMax <= self.relativeStepSizeTol && omegaSlp >= -1/self.maxIterations && deltaSlp >= -1/self.maxIterations
                    replab.msg(2, 'Stop: converged');
                    exitFlag = 1;
                end
            end
            if k >= self.maxIterations && exitFlag == 0
                replab.msg(2, 'Stop: number of iterations exceeded');
                exitFlag = -1;
            end
        end

    end

end

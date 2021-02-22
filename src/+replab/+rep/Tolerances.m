classdef Tolerances < replab.Str
% Tolerance criteria used for the termination of iterative algorithms

    properties (SetAccess = protected)
        minIterations % (integer): Minimum number of iterations, default: 5
        maxIterations % (integer): Maximum number of iterations, default: 200
        regularizationWindowSize % (integer): Number of steps to include in the monotonicity check, default: 8
        optimalityStandardTol % (double): Tolerance below which the problem is considered solved succesfully, default: 1e-10
        optimalityFinalTol % (double) : Tolerance at which the iterations stop, default: 1e-15
        relativeStepSizeTol % (double): Relative step size below which the iterations stop, default: 0

    end

    methods

        function self = Tolerances(varargin)
            args = struct('minIterations', 5, 'maxIterations', 200, 'regularizationWindowSize', 5, ...
                          'optimalityStandardTol', 1e-10, 'optimalityFinalTol', 1e-15, ...
                          'relativeStepSizeTol', 0);
            args = replab.util.populateStruct(args, varargin);
            self.minIterations = args.minIterations;
            self.maxIterations = args.maxIterations;
            self.regularizationWindowSize = args.regularizationWindowSize;
            self.optimalityStandardTol = args.optimalityStandardTol;
            self.optimalityFinalTol = args.optimalityFinalTol;
            self.relativeStepSizeTol = args.relativeStepSizeTol;
        end

        function logHeader(self)
        % Logs a header at the verbosity level used during the iterations
            replab.msg(2, ' #iter   optim.   (reg)      stepSize (reg)');
            replab.msg(2, '------------------------------------------------');
        end

        function exitFlag = test(self, omega, delta, k)
        % Tests the termination criteria
        %
        % Args:
        %   omega (double(1,k)): Optimality values
        %   delta (double(1,k)): Relative step sizes
        %   k (integer): Current iteration number
            exitFlag = 0;
            R = self.regularizationWindowSize;
            if k < R
                replab.msg(2, '%6d   %6.2E            %6.2E', k, omega(k), delta(k));
            elseif k >= R
                omegaMax = max(omega(k-R+1:k));
                deltaMax = max(delta(k-R+1:k));
                replab.msg(2, '%6d   %6.2E (%6.2E) %6.2E (%6.2E)', k, omega(k), omegaMax, delta(k), deltaMax);
                if k < self.minIterations
                    return
                end
                if omegaMax <= self.optimalityFinalTol
                    replab.msg(2, 'Stop: reached final optimality');
                    exitFlag = 1;
                elseif omegaMax <= self.optimalityStandardTol
                    if deltaMax <= self.relativeStepSizeTol
                        replab.msg(2, 'Stop: step size too small, and acceptable optimality');
                        exitFlag = 2;
                    elseif k > R
                        omegaPol = polyfit(1:R, log10(omega(k-R+1:k)), 1);
                        deltaPol = polyfit(1:R, log10(delta(k-R+1:k)), 1);
                        if omegaPol(1) >= 0 && deltaPol(1) >= 0
                            replab.msg(2, 'Stop: stalled, and acceptable optimality');
                            exitFlag = 3;
                        end
                    end
                end
            end
            if k >= self.maxIterations && exitFlag == 0
                replab.msg(2, 'Stop: number of iterations exceeded');
                exitFlag = -1;
            end
        end

    end

end

classdef Equivariant_forCompactGroup < replab.Equivariant

    methods

        function self = Equivariant_forCompactGroup(repR, repC, special)
            self@replab.Equivariant(repR, repC, special);
        end

    end

    methods (Access = protected)

        function [X, err] = project_double_sparse(self, X)
            X = full(X);
            % Parameters
            minIterations = 10;
            maxIterations = 200;
            relTol = 1e-6;
            relZero = 1e-15;
            windowSize = 8;
            nSamples = 3; % number of samples per iteration
            useInverses = true; % use inverses of samples as well
            dR = self.repR.dimension;
            dC = self.repC.dimension;
            replab.msg(1, 'Equivariant projection, %d x %d', dR, dC);
            replab.msg(1, '');
            replab.msg(2, ' #iter   delta    (slp)  norm');
            replab.msg(2, '---------------------------------');
            delta = zeros(1, maxIterations);
            exitFlag = 0;
            k = 1;
            comb = zeros(2, 0);
            for i = 1:windowSize
                for j = i+1:windowSize
                    comb(:, end+1) = [i;j];
                end
            end
            nX = norm(X, 'fro');
            while exitFlag == 0
                X1 = zeros(dR, dC);
                for j = 1:nSamples
                    g = self.group.sample;
                    S1 = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X));
                    if useInverses
                        ginv = self.group.inverse(g);
                        S2 = self.repR.matrixRowAction(ginv, self.repC.matrixColAction(ginv, X));
                        X1 = X1 + S1 + S2;
                    else
                        X1 = X1 + S1;
                    end
                end
                X1 = X1/nSamples;
                if useInverses
                    X1 = X1/2;
                end
                nX1 = norm(X1, 'fro');
                delta(k) = norm(X1 - X, 'fro')/nX;
                if k >= windowSize && k >= minIterations
                    window = delta(k-windowSize+1:k);
                    deltaMax = max(window);
                    deltaLog = log10(max(window, 1e-200));
                    % Theil–Sen estimator
                    iterDiff = comb(2,:) - comb(1,:);
                    deltaDiff = deltaLog(comb(2,:)) - deltaLog(comb(1,:));
                    slope = median(deltaDiff./iterDiff);
                    replab.msg(2, '%6d   %6.2E (%+1.1f) %6.2E', k, delta(k), slope, nX1/nX);
                    if nX1/nX <= relZero
                        exitFlag = 3;
                        replab.msg(1, 'Stop: relative norm of iterate is below %6.2E', relZero);
                    end
                    if all(window == 1e-100)
                        exitFlag = 2;
                        replab.msg(1, 'Stop: complete stall over the regularization window');
                    end
                    if slope >= -1/maxIterations && deltaMax <= relTol
                        exitFlag = 1;
                        replab.msg(1, 'Stop: estimated slope=%6.2E, and max delta over window=%6.2E', slope, deltaMax);
                    end
                else
                    replab.msg(2, '%6d   %6.2E        %6.2E', k, delta(k), nX1/nX);
                end
                if k >= maxIterations
                    exitFlag = -1;
                end
                X = X1;
                k = k + 1;
            end
            if exitFlag < 0
                error('replab:equivariantConvergence', 'Cannot satisfy convergence criterion in equivariant projection after %d iterations', k);
            end
            err = max(delta(k-windowSize+1:k));
        end

    end

end

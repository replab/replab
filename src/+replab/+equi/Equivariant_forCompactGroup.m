classdef Equivariant_forCompactGroup < replab.Equivariant

    methods

        function self = Equivariant_forCompactGroup(repR, repC, special)
            self@replab.Equivariant(repR, repC, special);
        end

        function B = blocks(self)
            B = self.cached('blocks', @() self.computeBlocks);
        end

    end

    methods (Access = protected)

        function B = computeBlocks(self)
            if ~self.repR.hasMaximalTorusExponents || ~self.repC.hasMaximalTorusExponents
                B = {1:self.repR.dimension; 1:self.repC.dimension};
                return
            end
            powersR = self.repR.maximalTorusExponents;
            powersC = self.repC.maximalTorusExponents;
            bR = find(all(powersR == 0, 2));
            bC = find(all(powersC == 0, 2));
            if ~isempty(bR) && ~isempty(bC)
                B = {bR; bC};
            else
                B = cell(2, 0);
            end
            r = 1;
            while r <= self.repR.dimension
                row = powersR(r, :);
                if any(row ~= 0)
                    bR = find(all(bsxfun(@eq, powersR, row), 2));
                    bC = find(all(bsxfun(@eq, powersC, row), 2));
                    B = horzcat(B, {bR; bC});
                    powersR(bR, :) = 0;
                    powersC(bC, :) = 0;
                end
                r = r + 1;
            end
        end

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
            delta = zeros(1, maxIterations);
            exitFlag = 0;
            k = 1;
            comb = zeros(2, 0);
            for i = 1:windowSize
                for j = i+1:windowSize
                    comb(:, end+1) = [i;j];
                end
            end
            group = self.group;
            repR = self.repR;
            repC = self.repC;
            if replab.globals.useReconstruction && repR.overC && repC.overC && repR.hasMaximalTorusExponents && repC.hasMaximalTorusExponents
                blocks = self.blocks;
                [~, R] = group.reconstruction;
                useTorus = true;
            else
                useTorus = false;
            end
            replab.msg(1, 'Equivariant projection, %d x %d', dR, dC);
            if useTorus
                replab.msg(1, 'Using maximal torus block projection');
            else
                replab.msg(1, 'Using standard group averaging');
            end
            replab.msg(1, '');
            replab.msg(2, ' #iter   delta    (slp)  norm');
            replab.msg(2, '---------------------------------');
            nX = norm(X, 'fro');
            if nX == 0
                err = 0;
                return
            end
            while exitFlag == 0
                if useTorus
                    X0 = zeros(dR, dC);
                    for i = 1:size(blocks, 2)
                        bR = blocks{1, i};
                        bC = blocks{2, i};
                        X0(bR, bC) = X(bR, bC);
                    end
                else
                    X0 = X;
                end
                X1 = zeros(dR, dC);
                for j = 1:nSamples
                    g = group.sample;
                    S1 = repR.matrixRowAction(g, repC.matrixColAction(g, X0));
                    if useInverses
                        ginv = group.inverse(g);
                        S2 = repR.matrixRowAction(ginv, repC.matrixColAction(ginv, X0));
                        X1 = X1 + S1 + S2;
                    else
                        X1 = X1 + S1;
                    end
                end
                X1 = X1/nSamples;
                if useInverses
                    X1 = X1/2;
                end
                if useTorus
                    for i = length(R.sets):-1:1
                        S = R.sets{i};
                        X2 = X1;
                        for j = 2:length(S)
                            g = S{j};
                            X2 = X2 + repR.matrixRowAction(g, repC.matrixColAction(g, X1));
                        end
                        X1 = X2/length(S);
                    end
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

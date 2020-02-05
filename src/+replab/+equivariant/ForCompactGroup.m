classdef ForCompactGroup < replab.Equivariant

    methods

        function self = ForCompactGroup(repR, repC, special)
            assert(ismember(exist('nlinfit'), [2 6]), 'Computations on compact groups require nonlinear curve fit.');
            self = self@replab.Equivariant(repR, repC, special);
        end

        function [X err] = project(self, X)
            hasNewFigure = false;

            % Parameters
            nSamplesPerIter = 3; % number of samples per iteration
            useInverses = true; % use inverses of samples as well
            maxIters = 5000; % maximum number of iterations before erroring
            nWarmUpIters = 20; % warm up iterations before first fit
            nCheckIters = 10; % iterations to verify the noise floor
            maxNoiseFloorError = 1; % maximum error on the noise floor (in log10)
            noiseMargin = 1; % noise floor margin (in log10)
            iter = 1;
            errs = [];
            % we fit 10^logfloor + exp(-slope*x + offset)
            modelfun = @(b,x) log10(10^b(1)+10.^(-b(2)*x+b(3)));

            for iter = 1:maxIters
                acc = zeros(self.repR.dimension, self.repC.dimension);
                for j = 1:nSamplesPerIter
                    g = self.group.sample;
                    S1 = self.repR.matrixRowAction(g, self.repC.matrixColAction(g, X));
                    if useInverses
                        ginv = self.group.inverse(g);
                        S2 = self.repR.matrixRowAction(ginv, self.repC.matrixColAction(ginv, X));
                        acc = acc + S1 + S2;
                    else
                        acc = acc + S1;
                    end
                end
                if useInverses
                    X1 = acc/(2*nSamplesPerIter);
                    d = norm(X1 - X, 'fro');
                    X = X1;
                else
                    X1 = acc/nSamplesPerIter;
                    d = norm(X1 - X, 'fro');
                    X = X1;
                end
                errs(1, iter) = d + 1e-100;
                if iter > nWarmUpIters
                    logfloor0 = log10(min(errs));
                    slope0 = (log10(errs(1)) - log10(errs(nWarmUpIters)))/nWarmUpIters;
                    offset0 = log10(errs(1));
                    beta0 = [logfloor0 slope0 offset0];
                    w = warning('off');
                    errored = false;
                    try
                        [beta, R, J, CovB, MSE] = nlinfit(1:iter, log10(errs(1:iter)), modelfun, beta0);
                    catch
                        errored = true;
                    end
                    warning(w);
                    if errored || rank(J) < 3
                        % the rank is deficient, usually because the curve didn't flatten
                        % yet, and thus the noise floor cannot be estimated
                    else
                        logfloor = beta(1);
                        slope = beta(2);
                        offset = beta(3);
                        % compute error on fit parameters
                        delta = sqrt(diag(CovB)) * tinv(1-0.05/2, iter - 3);
                        crossing = [];
                        if delta < maxNoiseFloorError/2
                            % -slope*x + offset == logfloor
                            % x = - (logfloor - offset)/slope
                            crossing = round(-(logfloor - offset)/slope);
                            if iter > crossing + nCheckIters
                                check = log10(errs((iter-nCheckIters):iter));
                                if all(check < logfloor + noiseMargin)
                                    err = logfloor + noiseMargin;
                                    return
                                end
                            end
                        end
                        if replab.equivariant.plotConvergence || i == maxIters
                            if ~hasNewFigure
                                figure
                                hasNewFigure = true;
                            end
                            clf
                            hold on
                            plot(1:iter, log10(errs(1:iter)), 'x');
                            plot(1:0.1:iter, modelfun(beta, 1:0.1:iter), '-');
                            if ~isempty(crossing)
                                plot([crossing crossing], [log10(min(errs)) log10(max(errs))], 'r-');
                            end
                            xlabel('Iteration #');
                            ylabel('Log10 of approximation error');
                            drawnow
                        end
                    end
                end
            end
            error('replab:equivariantConvergence', 'Cannot satisfy convergence criterion in equivariant projection after %d iterations', iter);
        end

    end

end

classdef Equivariant_forCompactGroup < replab.Equivariant

    methods

        function self = Equivariant_forCompactGroup(repR, repC, special)
            self@replab.Equivariant(repR, repC, special);
        end

    end

    methods (Access = protected)

        function [X, err] = project_double_sparse(self, X)
            X = full(X);
            hasNewFigure = false;
            % Parameters
            maxIters = 5000; % maximum number of iterations before erroring
            nSamplesPerIter = 3; % number of samples per iteration
            useInverses = true; % use inverses of samples as well
            nWarmUpIters = 20; % warm up iterations before first fit
            nCheckIters = 10; % number of the last iterations on which to verify the noise floor and the parameter errors
            nbSigmasUncertainties = 3; % Statistical confidence at which we estimate the uncertainty of the model parameters
            uncertaintiesConfidenceLevel = 1-(1-erf(nbSigmasUncertainties/sqrt(2)))/2; % Corresponding threshold probability
            maxParameterRelativeError = 0.1; % maximum relative error on the model parameters (defines if the parameters of the fit can be considered as meaningful)
            nbSigmasCrossing = 5; % statistical confidence that the crossing point was overcome
            maxNoiseFloorAbsoluteError = 1; % maximum absolute error on the noise floor (in log10)
            noiseMargin = 1; % noise floor margin (in log10)

            % we fit 10^logfloor + exp(-slope*x + offset)
            modelfun = @(b,x) log10(10^min(b(1),100)+10.^(-b(2)*x+b(3)));
            % For comparison purpose, we also model a purely exponential
            % fit
            modelfunExp = @(c,x) -c(1)*x+c(2);

            errs = [];
            modelRelativeErrors = [];
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
                    y = log10(errs(1:iter));
                    logfloor0 = log10(min(errs));
                    slope0 = (log10(errs(1)) - log10(errs(nWarmUpIters)))/nWarmUpIters;
                    offset0 = log10(errs(1));
                    beta0 = [logfloor0;slope0;offset0];
                    w = warning('off');
                    errored = false;
                    try
                        [beta, R, J, CovB, MSE] = replab.numerical.nlinfit((1:iter)', y(:), modelfun, beta0);
                        beta = beta(:).';
                    catch
                        errored = true;
                    end
                    warning(w);

                    if errored || (rank(J) < 3)
                        % the rank is deficient, usually because the curve didn't flatten
                        % yet, and thus the noise floor cannot be estimated
                    else
                        % extract the fit parameters and crossing point
                        logfloor = beta(1);
                        slope = beta(2);
                        offset = beta(3);
                        crossing = round(-(logfloor - offset)/slope);

                        % compute approximate uncertainty on these
                        % parameters due to finite sampling
                        delta = sqrt(diag(CovB)) * tinv(1-uncertaintiesConfidenceLevel, iter - 3);
                        jacobian = [-1/slope; -(logfloor - offset)/slope^2; 1/slope];
                        variance = jacobian'*CovB*jacobian;
                        crossingDeltaCov = sqrt(variance)*tinv(1-uncertaintiesConfidenceLevel, iter - 3);

                        % We monitor the errors
                        modelRelativeErrors(iter,:) = abs([delta(1:2)'./beta(1:2) crossingDeltaCov/crossing]);

                        % Depending on the fit parameters, we check if we
                        % performed enough iterations. For this we have the
                        % following conditions:
                        %  1. The noise floor level must be estimated
                        %     with sufficient accuracy
                        preciseNoiseFloor = delta(1) < maxNoiseFloorAbsoluteError/2;
                        %  2. The last few errors must be close enough to
                        %     the logfloor
                        lastErrors = log10(errs((iter-nCheckIters):iter));
                        lastErrorsOk = all(lastErrors < logfloor + noiseMargin);
                        %  3. The crossing point must be sufficiently well
                        %     defined
                        wellDefinedCrossingPoint = max(max(modelRelativeErrors((iter-nCheckIters):iter,:))) < maxParameterRelativeError;
                        %  4. We must be 'clearly' beyond the crossing
                        %     point
                        offset2 = offset + nbSigmasCrossing*std(R);
                        crossingThreshold = round(-(logfloor - offset2)/slope)-crossing;
                        minIter = crossing + max(nCheckIters, crossingThreshold);
                        minIterOk = (iter > minIter);
                        % All in all:
                        exitCondition = preciseNoiseFloor && lastErrorsOk && wellDefinedCrossingPoint && minIterOk;

                        % Before possibly exiting we eventually produce
                        % some plots
                        if (replab.globals.equivariantPlotConvergence && ((mod(iter, 20) == 1)  || exitCondition)) || (iter == maxIters)
                            if ~hasNewFigure
                                figure
                                hasNewFigure = true;
                            end
                            clf
                            subplot(2,1,1);
                            hold on
                            plot(1:iter, log10(errs(1:iter)), 'x');
                            plot(1:0.1:iter, modelfun(beta, 1:0.1:iter), '-');
                            plot(1:0.1:iter, modelfunExp([slope, offset2], 1:0.1:iter), '-');
                            plot(crossing*[1 1], [log10(min(errs)) log10(max(errs))], 'r-');
                            plot((crossing + crossingThreshold)*[1 1], [log10(min(errs)) log10(max(errs))], 'b-');
                            xlabel('Iteration #');
                            ylabel('Log10 of approximation error');
                            subplot(2,1,2);
                            hold on
                            plot(1:iter, modelRelativeErrors)
                            xlabel('Iteration #');
                            ylabel('Parameter relative uncertainty')
                            axis([1 iter -1 1]);
                            drawnow;
                            pause(0.01);
                        end

                        if exitCondition
                            err = max(lastErrors);
                            return
                        end
                    end
                end
            end
            error('replab:equivariantConvergence', 'Cannot satisfy convergence criterion in equivariant projection after %d iterations', iter);
        end

    end

end

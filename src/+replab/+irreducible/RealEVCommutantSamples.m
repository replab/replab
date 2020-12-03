classdef RealEVCommutantSamples < replab.Str

    properties (SetAccess = protected)
        rep % (`+replab.Rep`): Representation from which to take samples
    end

    properties (Access = protected)
        samples % (cell(1,\*) of double(\*,\*)): Samples
        errors % (double(1,\*)): Errors
    end

    methods

        function self = RealEVCommutantSamples(rep)
            self.rep = rep;
            self.samples = cell(1, 0);
            self.errors = zeros(1, 0);
        end

        function [X, err] = sample(self, ind)
            for i = length(self.samples)+1:ind
                rep = self.rep;
                dom = replab.domain.SelfAdjointMatrices(rep.field, rep.dimension);
                X = dom.sample;
                if ~rep.knownUnitary
                    U = rep.unitarize;
                    X = U.Ainv * X * U.A;
                end
                [X, err] = rep.commutant.project(X);
                self.samples{i} = X;
                self.errors(i) = err;
            end
            X = self.samples{ind};
            err = self.errors(ind);
        end

    end

end

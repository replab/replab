classdef Samples < replab.Str
% Lazy linked list of samples from an equivariant space

    properties (SetAccess = protected)
        E % (`+replab.Equivariant`): Equivariant space
    end

    properties (Access = protected)
        sample_ % (double(*,*)): Sampled matrix
        error_ % (double): Estimated error
        next_ % (`+replab.+equivariant.Samples`): Next samples in chain
    end

    methods (Access = protected)

        function computeValue(self)
            [sample error] = self.E.sampleWithError;
            self.sample_ = sample;
            self.error_ = error;
        end

        function computeNext(self)
            self.next_ = replab.equivariant.Samples(self.E);
        end

    end

    methods

        function self = Samples(E)
            assert(isa(E, 'replab.Equivariant'));
            self.E = E;
        end

        function [X err] = value(self)
        % Returns the current sample
            if isequal(self.sample_, [])
                self.computeValue;
            end
            X = self.sample_;
            err = self.error_;
        end

        function s = next(self)
        % Returns the next sample in chain
            if isequal(self.next_, [])
                self.computeNext;
            end
            s = self.next_;
        end

    end

end

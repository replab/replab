classdef Samples < replab.Str
% Lazy linked list of samples from an equivariant space

    properties (SetAccess = protected)
        E % (`+replab.Equivariant`): Equivariant space
    end

    properties (Access = protected)
        headSample_ % (double(*,*)): Sampled matrix
        headError_ % (double): Estimated error
        tail_ % (`+replab.+equivariant.Samples`): Next sample in chain
    end

    methods (Access = protected)

        function computeHead(self)
            [headSample headError] = self.E.sampleWithError;
            self.headSample_ = headSample;
            self.headError_ = headError;
        end

        function computeTail(self)
            self.tail_ = replab.equivariant.Samples(self.E);
        end

    end

    methods

        function self = Samples(E)
            assert(isa(E, 'replab.Equivariant'));
            self.E = E;
        end

        function [X err] = head(self)
        % Returns the current sample
            if isequal(self.headSample_, [])
                self.computeHead;
            end
            X = self.headSample_;
            err = self.headError_;
        end

        function s = tail(self)
        % Returns the next sample in chain
            if isequal(self.tail_, [])
                self.computeTail;
            end
            s = self.tail_;
        end

    end

end

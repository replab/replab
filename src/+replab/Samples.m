classdef Samples < replab.Obj
% Sequence of random samples from a domain

    properties (SetAccess = protected)
        domain % (`.Domain`): Domain to sample
    end

    properties (Access = protected)
        samples % (cell(1,\*) of samples): Already computed samples
    end

    methods

        function self = Samples(domain)
            self.domain = domain;
            self.samples = cell(1, 0);
        end

        function I = iterator(self)
        % Returns a mutable iterator into this sequence of samples
        %
        % Returns:
        %   `.SamplesIterator`: Iterator starting at the first sample
            I = replab.SamplesIterator(self);
        end

        function s = at(self, ind)
        % Returns the sample in the sequence at a particular index, and computes the samples up to that index if necessary
        %
        % Args:
        %   ind (integer): Positive integer indexing the samples
        %
        % Returns:
        %   Sample
            if length(self.samples) < ind
                for i = length(self.samples)+1:ind
                    self.samples{1,i} = self.domain.sample;
                end
            end
            s = self.samples{ind};
        end

    end

end

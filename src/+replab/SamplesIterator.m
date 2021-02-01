classdef SamplesIterator < replab.Str
% (Mutable) iterator into a sequence of samples

    properties (SetAccess = protected)
        samples % (`.Samples`): Samples
        last % (integer): Index of the last returned sample
    end

    methods

        function self = SamplesIterator(samples)
            self.samples = samples;
            self.last = 0;
        end

        function s = next(self)
        % Returns the next sample in the sequence and advances the iterator
        %
        % Returns:
        %   Sample
            ind = self.last + 1;
            s = self.samples.at(ind);
            self.last = ind;
        end

    end

end

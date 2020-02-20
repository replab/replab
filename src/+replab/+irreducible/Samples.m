classdef Samples < replab.Str
% Computes on-demand and caches samples from the commutant or trivial space of a representation
%
% Abstract base class

    properties
        rep % (`+replab.Rep`): Representation to compute things from
    end

    methods

        function s = forSubRep(self, sub)
            s = replab.irreducible.SubRepSamples(sub, self.rep, self);
        end

    end

    methods % Abstract

        function X = commutantSample(self, i)
        % Returns or computes the i-th sample from the commutant in the sample sequence
        %
        % Args:
        %   i (integer): 1-based index of the sample
        %
        % Returns:
        %   double(*,*): Sample
            error('Abstract');
        end

        function X = trivialSample(self, i)
        % Returns or computes the i-th samples from the trivial space in the sample sequence
        %
        % Args:
        %   i (integer): 1-based index of the sample
        %
        % Returns:
        %   double(*,*): Matrix whose columns are elements of the trivial space
            error('Abstract');
        end

    end

end

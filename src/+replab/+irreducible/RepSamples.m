classdef RepSamples < replab.Str
% Computes and caches samples from the commutant or trivial space of a representation
%
% Uses instances of `+replab.Equivariant` to compute the samples.

    properties (SetAccess = protected)
        commutantSamples % (cell(1,*) of double(*,*)): Samples from the commutant
        trivialSamples % (cell(1,*) of double(*,*)): Samples from the trivial space
    end

    properties (Access = protected)
        trivial_ % (`+replab.Equivariant`): Equivariant space from rep to the trivial representation of size ``rep.dimension``
    end

    methods

        function self = RepSamples(rep)
            self.rep = rep;
            self.commutantSamples = {};
            self.trivialSamples = {};
        end

        function T = trivial(self)
        % Returns a "window" into the trivial space of a representation
        %
        % It is defined as the space of equivariant maps from the underlying representation
        % to the trivial representation of same dimension.
        %
        % Returns:
        %   `+replab.Equivariant`: Trivial space as an equivariant space
            if isempty(self.trivial_)
                trivialRep = self.rep.group.trivialRep(self.rep.field, self.rep.dimension);
                self.trivial_ = trivialRep.equivariant(self.rep);
            end
            T = self.trivial_;
        end

        function X = commutantSample(self, i)
        % Returns or computes the i-th sample from the commutant in the sample sequence
        %
        % Args:
        %   i (integer): 1-based index of the sample
        %
        % Returns:
        %   double(*,*): Sample
            while length(self.commutantSamples) < i
                self.commutantSamples{1, end+1} = self.rep.commutant.sample;
            end
            X = self.commutantSamples{i};
        end

        function X = commutantSampleForSubRep(self, i, sub)
        % Returns or computes the i-th sample for a subrepresentation from the commutant
        %
        % Args:
        %   i (integer): 1-based index of the sample
        %   sub (`+replab.SubRep`): Subrepresentation of `.rep`
        %
        % Returns:
        %   double(*,*): Sample
            assert(sub.parent == self.rep);
            X = full(sub.F_internal * self.commutantSample(i) * self.H_internal);
        end

        function X = trivialSample(self, i)
        % Returns or computes the i-th samples from the trivial sapce in the sample sequence
        %
        % Args:
        %   i (integer): 1-based index of the sample
        %
        % Returns:
        %   double(*,*): Sample
            while length(self.trivialSamples) < i
                self.trivialSamples{1, end+1} = self.trivial.sample;
            end
            X = self.trivialSamples{i};
        end

    end

end

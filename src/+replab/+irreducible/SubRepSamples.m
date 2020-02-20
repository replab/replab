classdef SubRepSamples < replab.irreducible < Samples
% Computes and caches samples from the commutant or trivial space of a representation
%
% Reuses the on-demand samples of a parent representation

    properties
        parent % (`+replab.Rep`): Parent representation
        parentSamples % (`+replab.+irreducible.Samples`): Samples for the parent representation
    end

    methods

        function self = SubRepSamples(rep, parent, parentSamples)
        %
        % Args:
        %   rep (`+replab.SubRep`): Subrepresentation of the parent representation
        %   parent (`+replab.Rep`): Parent representation
        %   parentSamples (`+replab.+irreducible.Samples`): Samples for the parent representation
            assert(isa(rep, 'replab.SubRep'));
            assert(rep.parent == parent);
            assert(parentSamples.rep == parent);
            self.rep = rep;
            self.parent = parent;
            self.parentSamples = parentSamples;
        end

        function X = commutantSample(self, i)
            X = full(self.rep.F_internal * self.parentSamples.commutantSample(i) * self.rep.H_internal);
        end

        function X = trivialSample(self, i)
            d = self.rep.dimension;
            dP = self.parent.dimension;
            cut = sparse(1:d, 1:d, ones(1, d), dP, d);
            X = full(self.rep.F_internal * self.parentSamples.trivialSample(i) * cut);
        end

    end

end

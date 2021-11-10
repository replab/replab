classdef NiceIsomorphism < replab.gen.NiceIsomorphism

    properties (SetAccess = protected)
        domainSize % (integer): Size of the signed permutations
    end

    methods

        function self = NiceIsomorphism(domainSize, sourceType)
            self.domainSize = domainSize;
            sourceGenerators = {[-1 2:domainSize]};
            if domainSize > 1
                sourceGenerators{1,end+1} = [2:domainSize 1];
            end
            if domainSize > 2
                sourceGenerators{1,end+1} = [2 1 3:domainSize];
            end
            targetType = replab.PermutationGroupType.make(2*domainSize);
            self.finishConstruction(@(niceIso) replab.SignedPermutationGroup(domainSize, sourceGenerators, 'type', sourceType, 'niceIsomorphism', niceIso), sourceGenerators, targetType);
        end

    end

    methods % Implementations

        % Morphism

        function t = imageElement(self, s)
            t = replab.SignedPermutation.toPermutation(s);
        end

        % Isomorphism

        function s = preimageElement(self, t)
            s = replab.SignedPermutation.fromPermutation(t);
        end

        % NiceIsomorphism

        function l = sourceContains(self, s)
            l = true;
        end

    end

end

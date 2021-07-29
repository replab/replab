classdef SignedPermutationIsomorphism < replab.GeneralIsomorphism

    properties (SetAccess = true)
        domainSize % (integer): Size of the signed permutations
    end

    methods

        function self = SignedPermutationIsomorphism(signedSymmetricGroup)
            self.domainSize = domainSize;

            targetGenerators = {[2 1 3:2*domainSize]};
            if domainSize > 1
                targetGenerators{1,end+1} = [3:2*domainSize 1 2];
            end
            if domainSize > 2
                targetGenerators{1,end+1} = [3 4 1 2 5:2*domainSize];
            end
            self.source =
            self.target =
            self.torusMap = [];
        end

    end

end

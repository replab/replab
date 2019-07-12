classdef GHZ < replab.SemidirectProductGroup
% An approximation of the symmetry group of GHZ states    
    properties
        nParties;
        nLevels;
        rootOrder;
    end
    
    methods
        
        function self = GHZ(nParties, nLevels, rootOrder)
            base = replab.quantum.GHZBase(nParties, nLevels, rootOrder);
            SParties = replab.S(nParties);
            SLevels = replab.S(nLevels);
            quotient = replab.DirectProductGroup({SParties SLevels});
            f = @(q, b) base.permuteParties(q{1}, base.permuteLevels(q{2}, b));
            phi = replab.ActionFun('Permutation of parties/levels', quotient, base, f);
            G = replab.SemidirectProductGroup(phi);
            self = self@replab.SemidirectProductGroup(phi);
            self.nParties = nParties;
            self.nLevels = nLevels;
            self.rootOrder = rootOrder;
        end
        
        function rho = toMatrix(self, g)
            h = g{1};
            h1 = h{1};
            h2 = h{2};
            H1 = self.H.factor(1);
            H2 = self.H.factor(2);
            nL = self.nLevels;
            partyRho = replab.Permutations.toMatrix(H1.indexRelabelingPermutation(h1, nL));
            levelMat = replab.Permutations.toMatrix(h2);
            levelRho = 1;
            for i = 1:self.nParties
                levelRho = kron(levelRho, levelMat);
            end
            phaseRho = self.N.toMatrix(g{2});
            rho = partyRho*levelRho*phaseRho;
        end
        
        function rep = rep(self)
            rep = replab.RepFun(self, 'C', self.N.naturalRep.dimension, @(g) self.toMatrix(g));
        end
        
    end
    
end
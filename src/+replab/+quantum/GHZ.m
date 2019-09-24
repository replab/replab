classdef GHZ < replab.semidirectproduct.OfCompactGroups
% Symmetry group of the GHZ states
    
    properties
        nParties % integer: Number of parties
        nLevels % integer: Number of levels for each party (=2 for qubits)
    end
    
    methods
        
        function self = GHZ(nParties, nLevels)
        % Constructs the GHZ group for a given number of parties and levels
        %
        % Args:
        %   nParties (integer): Number of parties
        %   nLevels (integer): Number of levels
            base = replab.quantum.GHZBase(nParties, nLevels);
            SParties = replab.S(nParties);
            SLevels = replab.S(nLevels);
            discrete = directProduct(SParties, SLevels);
            f = @(q, b) base.permuteParties(q{1}, base.permuteLevels(q{2}, b));
            phi = replab.Action.lambda('Permutation of parties/levels', discrete, base, f);
            self@replab.semidirectproduct.OfCompactGroups(phi);
            self.nParties = nParties;
            self.nLevels = nLevels;
        end
        
        function rho = toMatrix(self, g)
        % Returns the natural matrix action of a group element
        %
        % It represents the action on a tensor space (C^d)^n, where
        % d = self.nLevels and n = self.nParties.
        %
        % For the action of the phase part (connected group part), 
        % see `replab.quantum.GHZBase.toMatrix`. This action is complemented
        % by the permutation of parties, and the correlated permutation of subsystem levels.
        %
        % Args:
        %   g (element): Group element
        %
        % Returns:
        %   double matrix: Unitary matrix representing `g`
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
        
        function rep = naturalRep(self)
        % Returns the natural representation of this group
        %
        % Returns:
        %   replab.Rep: The unitary natural representation
            rep = replab.Rep.lambda(self, 'C', self.N.naturalRep.dimension, true, @(g) self.toMatrix(g));
        end
        
    end
    
end

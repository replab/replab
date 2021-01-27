classdef LexMin < replab.Obj
% Finds the minimal lexicographic representative of a list under the action of a permutation group

    properties (SetAccess = protected)
        degree % (integer): Domain size, length of `.vec`
        lexChain % (`+replab.+bsgs.Chain`): BSGS chain with base in lexicographic order
        G % (`.PermutationGroup`): Group to search in
        vec % (double(1,\*)): Vector to find the minimal representative of
        K % (`.PermutationGroup`): Subgroup of `.G` that stabilizes `.vec`
    end

    methods

        function self = LexMin(G, vec, K)
            assert(G.domainSize == length(vec));
            self.G = G;
            self.lexChain = G.lexChain;
            self.vec = vec;
            self.K = K;
            self.degree = G.domainSize;
        end

        function [minimal, g] = search(self)
            self.minimal = zeros(1, self.degree);
            self.minimalCorrectBefore = 1;
            self.minimalG = 1:self.degree;
            for level = 1:length(self.lexChain.B)
                self.rec(1, level, 1:self.degree);
            end
            g = self.G.inverse(self.minimalG);
            minimal = self.minimal;
        end

    end

    properties (Access = protected)
        minimal
        minimalCorrectBefore
        minimalG
    end

    methods (Access = protected)

        function rec(self, level, toLevel, curG)
            if level <= toLevel
                chain = self.lexChain;
                candidates = [];
                beta = chain.B(level);
                if level < length(chain.B)
                    nextBeta = chain.B(level+1);
                else
                    nextBeta = self.degree+1;
                end
                if nextBeta > self.minimalCorrectBefore
                    for k = self.minimalCorrectBefore:nextBeta-1
                        self.minimal(k) = self.vec(self.minimalG(k));
                    end
                    self.minimalCorrectBefore = nextBeta;
                end
                orbit = chain.Delta{level};
                for oi = 1:length(orbit)
                    b = orbit(oi);
                    bg = curG(b);
                    comp = self.vec(bg) - self.minimal(beta);
                    k = beta + 1;
                    nextG = curG(chain.U{level}(:,oi));
                    while k < nextBeta && comp == 0
                        comp = self.vec(nextG(k)) - self.minimal(k);
                        k = k + 1;
                    end
                    if comp <= 0
                        if comp < 0
                            for k = beta:self.minimalCorrectBefore-1
                                self.minimal(k) = self.vec(nextG(k));
                            end
                            self.minimalG = nextG;
                            candidates = [];
                        end
                        candidates(end+1) = oi;
                    end
                end
                for oi = candidates
                    b = orbit(oi);
                    bg = curG(b);
                    nextG = curG(chain.U{level}(:,oi));
                    self.rec(level + 1, toLevel, nextG);
                end
            end
        end

    end

end

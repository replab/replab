classdef LexMin < replab.Obj
% Finds the minimal lexicographic representative of a list under the action of a permutation group

    properties (SetAccess = protected)
        degree % (integer): Domain size, length of `.vec`
        lexChain % (`+replab.+bsgs.Chain`): BSGS chain with base in lexicographic order
        G % (`.PermutationGroup`): Group to search in
        vec % (double(1,\*)): Vector to find the minimal representative of
        vecStabilizer % (`.PermutationGroup`): Subgroup of `.G` that stabilizes `.vec`
    end

    methods

        function self = LexMin(G, vec, vecStabilizer)
            assert(G.domainSize == length(vec));
            self.G = G;
            self.lexChain = G.lexChain;
            self.vec = vec;
            self.vecStabilizer = vecStabilizer;
            self.degree = G.domainSize;
        end

        function [minimal, g] = search(self)
        % Performs the backtrack search
            self.minimal = zeros(1, self.degree); % initialize the variables
            self.minimalCorrectBefore = 1;
            self.minimalG = 1:self.degree;
            for level = 1:length(self.lexChain.B)
                % find lex-minimal prefix of length "level"
                self.rec(1, level, 1:self.degree, self.vecStabilizer.chain.mutableCopy, 1);
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

        function rec(self, level, toLevel, curG, stab, stabLevel)
        % Implements breath-first in the cosets ``G / vecStabilizer``
        %
        % It filters elements that do not lead to a minimal lexicographic representative at each step in the stabilizer chain.
        % It also searchs only one element per coset.
        %
        % Args:
        %   level (integer): Current level
        %   toLevel (integer): Maximum level to explore
        %   curG (permutation): Current product of transversals
        %   stab (`.Chain`): Mutable chain of stabilizers, with changed base
        %   stabLevel (integer): Current level in the stabilizer chain
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
                if isempty(stab)
                    mask = true(1, self.degree);
                else
                    mask = replab.bsgs.minimalMaskInOrbit(self.degree, stab.S(:,stab.Sind(stabLevel):end));
                end
                for oi = candidates
                    b = orbit(oi);
                    bg = curG(b);
                    if mask(bg)
                        if isempty(stab)
                            stab1 = stab;
                            stabLevel1 = 0;
                        elseif stab.stabilizes(stabLevel, bg)
                            % if it already stabilizes, pass it on unchanged
                            stab1 = stab;
                            stabLevel1 = stabLevel;
                        else
                            % otherwise stabilize
                            stab.changeBasePointAt(stabLevel, bg);
                            stabLevel1 = stabLevel + 1;
                            os = stab.orbitSizes;
                            if all(os(stabLevel1:end) == 1)
                                stab1 = [];
                            else
                                stab1 = stab;
                            end
                        end
                        nextG = curG(chain.U{level}(:,oi));
                        self.rec(level + 1, toLevel, nextG, stab1, stabLevel1);
                    end
                end
            end
        end

    end

end

classdef Node < replab.bsgs.Chain
% A node with non-trivial orbit in the BSGS chain
    
    properties
        % Invariants:
        % we have leftAction(u(i), beta) = orbit(i)

        beta; % base point
        orbit; % orbit of beta stored as 1 x orbitSize cell array
        u; % transversal elements stored as 1 x orbitSize cell array
        uInv; % inverse of transversal elements stored as 1 x orbitSize cell array
        uInvW;
        newWord;
        emptyWords;
        ownSG; % strong generators found at this node, stored as 1 x nStrongGens cell array
        next; % next node in chain, or [] for the terminal node
    end
        
    methods
                
        function self = Node(A, beta)
            self.A = A;
            self.beta = beta;
            self.orbit = {beta};
            self.u = {A.G.identity};
            self.uInv = {A.G.identity};
            self.ownSG = {};
            self.uInvW = {replab.Word.identity};
            self.newWord = [false];
            self.emptyWords = 0;
        end
        
        function l = orbitSize(self)
            l = length(self.orbit);
        end
        
        function r = randomU(self)
            r = self.u{randi(self.orbitSize)};
        end
                
        function ind = orbitIndex(self, b)
        % Returns index of orbit element, or 0 if not found
            ind = 0;
            for i = 1:self.orbitSize
                if self.P.eqv(self.orbit{i}, b)
                    ind = i;
                    return
                end
            end
        end
       
        function addOwnStrongGenerator(self, g)
            self.ownSG{end + 1} = g;
        end
        
        function update(self)
            S = self.strongGeneratingSet;
            nS = length(S);
            i = 1;
            while i <= self.orbitSize
                b = self.orbit{i};
                for j = 1:length(S)
                    g = S{j};
                    newB = self.A.leftAction(g, b);
                    if self.orbitIndex(newB) == 0
                        newU = self.G.compose(g, self.u{i});
                        self.orbit{end+1} = newB;
                        self.u{end+1} = newU;
                        self.uInv{end+1} = self.G.inverse(newU);
                        self.uInvW{end+1} = [];
                        self.newWord(end+1) = false;
                        self.emptyWords = self.emptyWords + 1;
                    end
                end
                i = i + 1;
            end
        end
        
    end
    
end

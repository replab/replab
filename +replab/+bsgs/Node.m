classdef Node < replab.bsgs.Chain
% A node with non-trivial orbit in the BSGS chain
    
    properties
        % Invariants:
        % we have leftAction(u(i), beta) = orbit(i)

        beta; % base point
        orbit; % orbit of beta stored as 1 x orbitSize cell array
        u; % transversal elements stored as 1 x orbitSize cell array
        uInv; % inverse of transversal elements stored as 1 x orbitSize cell array
        uWords;
        uInvWords;
        ownSG; % strong generators found at this node, stored as 1 x nStrongGens cell array
        ownSGWords;
        next; % next node in chain, or [] for the terminal node
    end
        
    methods
        
        function check(self)
            beta = self.beta;
            for i = 1:length(self.orbit)
                u = self.u{i};
                uInv = self.uInv{i};
                b = self.orbit{i};
                self.P.assertEqv(self.A.leftAction(u, beta), b);
                self.P.assertEqv(self.A.leftAction(uInv, b), beta);
            end
        end
        
        function self = Node(A, beta)
            self.A = A;
            self.beta = beta;
            self.orbit = {beta};
            self.u = {A.G.identity};
            self.uInv = {A.G.identity};
            self.ownSG = {};
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
                    end
                end
                i = i + 1;
            end
        end
        
    end
    
end

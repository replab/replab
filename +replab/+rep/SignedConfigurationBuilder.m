classdef SignedConfigurationBuilder < handle
    
    properties
        n; % this configuration has size n x n
        nGenerators; % number of generators in the group
        generators;  % generators in a matrix, one per row
        invGenerators; % inverses of generators in a matrix, one per row
        fibers; % fibers of this configuration, a partition with F blocks
        fiberOrbits; % FxF cell array of orbit indices present in a block
        zeroStartRow; % if a zero orbit is known, its starting row,col, otherwise 0,0
        zeroStartCol;
        zeroLastRow; % last cell added to the zero orbitx 
        zeroLastCol;
        nOrbits; % number of orbits known so far
        orbitStartRow; % start cell of the i-th orbit, i = 1,...,nOrbits
        orbitStartCol;
        orbitIndex; % orbit index for the cell (r, c)
        orbitNextRow; % next element in the current orbit
        orbitNextCol;
        flip; % whether there is a sign flip relative to orbitStartRow/Col
    end
    
    methods
        
        function self = SignedConfigurationBuilder(generators)
        % Builds a SignedConfigurationBuilder from generators given
        % as signed permutations, one per row, in the matrix generators
            nGenerators = size(generators, 1);
            n = size(generators, 2);
            self.n = n;
            self.nGenerators = nGenerators;
            self.generators = int32(generators);
            self.invGenerators = zeros(nGenerators, n, 'int32');
            G = replab.SignedPermutations(n);
            for i = 1:nGenerators
                self.invGenerators(i, :) = G.inverse(self.generators(i, :));
            end
            self.fibers = replab.Partition.permutationsOrbits(abs(generators));
            self.fiberOrbits = cell(self.fibers.nBlocks, self.fibers.nBlocks);
            self.zeroStartRow = 0;
            self.zeroStartCol = 0;
            self.zeroLastRow = 0;
            self.zeroLastCol = 0;
            self.nOrbits = 0;
            self.orbitStartRow = [];
            self.orbitStartCol = [];
            self.orbitNextRow = zeros(self.n, self.n, 'int32');
            self.orbitNextCol = zeros(self.n, self.n, 'int32');
            self.orbitIndex = -ones(self.n, self.n, 'int32');
            self.flip = false(self.n, self.n);
            self.iterateFibers;
        end
        
        function C = toPhaseConfiguration(self)
            orbitStart = zeros(1, 0, 'int32');
            orbitRow = zeros(1, 0, 'int32');
            orbitCol = zeros(1, 0, 'int32');
            for o = 1:self.nOrbits
                r = self.orbitStartRow(o);
                c = self.orbitStartCol(o);
                orbitStart(1, o) = length(orbitRow) + 1;
                while r ~= 0
                    orbitRow(1, end+1) = r;
                    orbitCol(1, end+1) = c;
                    lr = r;
                    lc = c;
                    r = self.orbitNextRow(lr, lc);
                    c = self.orbitNextCol(lr, lc);
                end
            end
            orbitStart(1, self.nOrbits + 1) = length(orbitRow) + 1;
            C = replab.rep.PhaseConfiguration(self.n, self.fibers, self.fiberOrbits, 1-self.flip*2, orbitStart, orbitRow, orbitCol, self.orbitIndex);
        end
        
        function orbitOpt = iterate(self, r, c)
        % Iterates the orbit starting from the cell (r,c)
        % Returns [] if the orbit is zero, or [orbitInd] if nonzero
            T = 1;
            testR = [r];
            testC = [c];
            self.nOrbits = self.nOrbits + 1;
            o = self.nOrbits;
            self.orbitIndex(r, c) = o;
            self.orbitStartRow(o) = r;
            self.orbitStartCol(o) = c;
            self.flip(r, c) = false;
            lastR = r;
            lastC = c;
            isZero = false;
            while T > 0
                tr = testR(T);
                tc = testC(T);
                T = T - 1;
                for i = 1:self.nGenerators
                    ir = self.generators(i, tr);
                    ic = self.generators(i, tc);
                    s = sign(ir)*sign(ic) == -1;
                    ir = abs(ir);
                    ic = abs(ic);
                    oi = self.orbitIndex(ir, ic);
                    iflip = xor(s, self.flip(tr, tc));
                    switch oi
                      case -1 % unknown index, add to orbit
                        self.orbitIndex(ir, ic) = o;
                        self.orbitNextRow(lastR, lastC) = ir;
                        self.orbitNextCol(lastR, lastC) = ic;
                        self.flip(ir, ic) = iflip;
                        lastR = ir;
                        lastC = ic;
                        T = T + 1;
                        testR(T) = ir;
                        testC(T) = ic;
                      case 0
                        error('Should not happen');
                      case o
                        if self.flip(ir, ic) ~= iflip
                            isZero = true;
                        end
                      otherwise
                        error('Should not happen');
                    end
                end
            end
            % post-process zero orbit
            if isZero
                r = self.orbitStartRow(o);
                c = self.orbitStartCol(o);
                if self.zeroStartRow == 0
                    assert(self.zeroStartCol == 0); % no preexisting linked list for the zero orbit
                    self.zeroStartRow = r; % start of the orbit
                    self.zeroStartCol = c;
                    self.orbitIndex(r, c) = 0; % apply already one iteration
                    lr = r; lc = c;
                    r = self.orbitNextRow(lr, lc);
                    c = self.orbitNextCol(lr, lc);
                else
                    % merge the current orbit with the last zero orbit
                    lr = self.zeroLastRow;
                    lc = self.zeroLastCol;
                end
                while r ~= 0
                    self.orbitNextRow(lr, lc) = r;
                    self.orbitNextCol(lr, lc) = c;
                    self.orbitIndex(r, c) = 0;
                    lr = r; lc = c;
                    r = self.orbitNextRow(lr, lc);
                    c = self.orbitNextCol(lr, lc);
                end
                self.zeroLastRow = lr;
                self.zeroLastCol = lc;
                self.nOrbits = self.nOrbits - 1;
                orbitOpt = [];
            else
                orbitOpt = o;
            end
        end
        
        function iterateFibers(self)
            nFibers = self.fibers.nBlocks;
            for fr = 1:nFibers % iterate over row fibers
                for fc = 1:nFibers % iterate over col fibers
                    Fr = self.fibers.block(fr);
                    Fc = self.fibers.block(fc);
                    fo = zeros(1, 0);
                    for r = Fr
                        for c = Fc
                            if self.orbitIndex(r, c) == -1
                                orbitOpt = self.iterate(r, c);
                                if ~isempty(orbitOpt)
                                    fo(end+1) = orbitOpt;
                                end
                            end
                        end
                    end
                    self.fiberOrbits{fr, fc} = fo;
                end
            end
        end
        
    end
end

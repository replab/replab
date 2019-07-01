classdef SignedPhaseConfigurationBuilder < handle
    
    properties (SetAccess = protected)
        nR; % number of rows in phase configuration
        nC; % number of columns in phase configuration
        nGenerators; % number of generators in the group
        rowAction; % generator image, row action
                   % in a matrix, one per matrix row
        colAction; % generator image, col action
                   % in a matrix, one per matrix row
        zeroStartR; % if a zero orbit is known, its starting row,col, otherwise 0,0
        zeroStartC;
        zeroLastR; % last cell added to the zero orbit
        zeroLastC;
        nOrbits; % number of orbits known so far
        orbitStartR; % start cell of the i-th orbit, i = 1,...,nOrbits
        orbitStartC;
        orbitIndex; % orbit index for the cell (r, c)
        orbitNextR; % next element in the current orbit
        orbitNextC;
        flip; % whether there is a sign flip relative to orbitStartRow/Col
    end
    
    methods
        
        function self = SignedPhaseConfigurationBuilder(rowAction, colAction)
        % Builds a SignedConfigurationBuilder from action of
        % generators given as signed permutations, one per row, 
        % in the matrix rowAction and colAction
            nGenerators = size(rowAction, 1);
            assert(nGenerators == size(colAction, 1));
            self.nGenerators = nGenerators;
            
            nRows = size(rowAction, 2);
            nCols = size(colAction, 2);
            self.nRows = nRows;
            self.nCols = nCols;
            
            self.rowAction = int32(rowAction); % optimize access
            self.colAction = int32(colActino);
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
            for r = 1:nR
                for c = 1:nC
                    if self.orbitIndex(r, c) == -1
                        self.iterate(r, c);
                    end
                end
            end
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
                
    end
end

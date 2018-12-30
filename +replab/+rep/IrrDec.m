classdef IrrDec
% Describes a decomposition of an associative algebra in irreducible representations
% over the reals.
%
% TODO: it identifies but does not handle quaternionic representations.
    
    properties (SetAccess = protected)
        algebra;      % Associative algebra which we decompose
        fromFiber;    % fromFiber(i) is the index of fiber in from which
                      % the basis vector U(:,i) comes from        
        U;            % Orthonormal change of basis matrix
        nComponents;  % Number of isotypic components
        compDims;     % Isotypic component dimensions
        repDims;      % Representation dimensions
        repMuls;      % Representation multiplicities
                      % all types are real
        repTypes;     % 1 - real, 2 - complex, 3 - quaternionic
    end
    
    methods
               
        function self = IrrDec(algebra, fromFiber, U, repDims, repMuls, repTypes)
            assert(isreal(U));
            self.algebra = algebra;
            self.fromFiber = fromFiber;
            self.U = U;
            self.nComponents = length(repDims);
            self.repDims = repDims(:)';
            self.repMuls = repMuls(:)';
            self.compDims = repDims(:)' .* repMuls(:)';
            self.repTypes = repTypes(:)';
        end
        
        function R = compRange(self, r)
        % Indices corresponding to the r-th isotypic component
        % Correspond to columns of U, and to indices in fromOrbit
            from = sum(self.compDims(1:r-1)) + 1;
            to = sum(self.compDims(1:r));
            R = from:to;
        end
        
        function Urep = compBasis(self, r)
        % Returns the basis of the r-th isotypic component
            Urep = self.U(:, self.compRange(r));
        end

        function perm = swap(self)
        % Returns the permutation that switches between the kron(A, id) and kron(id, B) forms
            shift = 0;
            perm = [];
            for r = 1:self.nComponents
                d = self.repDims(r);
                m = self.repMuls(r);
                s = reshape(reshape(1:d*m, [d m])', [1 d*m]);
                perm = [perm, shift + s];
                shift = shift + d*m;
            end
        end
        
% $$$         function check(self)
% $$$         % Checks the validity of this irreducible decomposition
% $$$             import qdimsum.*
% $$$             % Perform isotypic checks
% $$$             self.toIsoDec.check;
% $$$             tol = self.settings.blockDiagMatTol;
% $$$             % Checks that the isotypic components are correct by considering
% $$$             % a sample from matrices that commute with the group
% $$$             sampleI = self.group.phaseConfiguration.sampleRealGaussian;
% $$$             collected = false;
% $$$             blocks = self.projectInIrrBasis(sampleI, collected, true);
% $$$             testI = self.U*blkdiag(blocks{:})*self.U';
% $$$             assert(~isNonZeroMatrix(testI - sampleI, tol), 'Matrix that commutes with the group has not the proper form');
% $$$             % Checks correctness by sampling from the group algebra
% $$$             sampleG = randn*GenPerm.orthogonalMatrix(self.group.randomElement) + randn*GenPerm.orthogonalMatrix(self.group.randomElement);
% $$$             collected = true;
% $$$             blocks = self.projectInIrrBasis(sampleG, collected, true);
% $$$             testG = self.U*blkdiag(blocks{:})*self.U';
% $$$             assert(~isNonZeroMatrix(testG - sampleG, tol), 'Group algebra matrix has not the proper form');
% $$$         end
        
    end

    methods (Static)
        
        function refinedBasis = orderedComplexBasis(iso, r)
        % Reorder the components of a complex representation so that the complex structure is visible
        % or returns [] if the representation is quaternionic
        % Optimized for precision
            tol = replab.Settings.doubleEigTol;
            range = iso.compRange(r);
            d = iso.repDims(r);
            m = iso.repMuls(r);
            newU = iso.U;
            fibersForRange = iso.fromFiber(range);
            refinedBasis = zeros(iso.algebra.n, length(range));
            % iterate over all fibers present in this rep
            for f = unique(fibersForRange) 
                basisInd = range(fibersForRange == f); % fiber indices for rep
                realRank = length(basisInd);
                % Euclidean coordinates of the f-th fiber elements (regardless of representation)
                fFiber = iso.algebra.fibers.block(f);
                n = length(fFiber);
                % find restriction of algebra to the f-th fiber
                resAlgebra = iso.algebra.fiber(f);
                % basis for the r-th representation in the f-th orbit
                basis = iso.U(fFiber, basisInd);
                % compute a generic invariant sample (non-symmetric matrix), restricted to fFiber x fFiber
                sample = basis*replab.rep.sampleRealMatrix(realRank, realRank)*basis';
                sample = resAlgebra.project(sample); % project in the invariant subspace
                assert(isreal(sample));
                % Perform the Schur decomposition
                [U T] = schur(sample, 'real');
                % Get diagonal part = real part of eigenvalues
                D = diag(T);
                % Find non zero eigenvalues => part of the isotypic component
                nzInd = find(abs(D) > tol);
                % Find the zero eigenvalues => part of the null space
                zInd = find(abs(D) <= tol);
                % Cluster eigenvalues by reordering the Schur blocks
                clusters = zeros(1, n);
                % Put null space last
                clusters(zInd) = 1;
                % Find similar eigenvalues
                nzD = D(nzInd(1:2:end));
                distD = abs(bsxfun(@minus, nzD, nzD'));
                maskD = distD <= tol;
                conD = replab.Partition.connectedComponents(maskD).blocks;
                % If the rank is not saturated, we have a quaternionic representation
                if length(conD) * 2 ~= realRank
                    % wrong rank found, it is a quaternionic representation
                    refinedBasis = [];
                    return
                end
                % Group similar eigenvalues
                for i = 1:length(conD)
                    compD1 = nzInd(conD{i}*2-1);
                    compD2 = nzInd(conD{i}*2);
                    assert(length(compD1) == d/2);
                    clusters(compD1) = i + 1;
                    clusters(compD2) = i + 1;
                end
                % Reorder
                [US TS] = ordschur(U, T, clusters);
                % Force blocks corresponding to the same complex eigenvalue to express it the same way
                % by removing the degeneracy due to conjugation
                start = 1;
                b = 1;
                for i = 1:m
                    b = b + 2;
                    for j = 2:n/m/2 % number of blocks
                        if sign(TS(b, b + 1)) == -sign(TS(start, start+1))
                            TS(:,[b b+1]) = TS(:,[b+1 b]);
                            TS([b b+1],:) = TS([b+1 b],:);
                            US(:,[b b+1]) = US(:,[b+1 b]);
                        end
                        b = b + 2;
                    end
                    start = b;
                end
                % Refined reordered basis is found
                refinedBasis(fFiber, fibersForRange == f) = US(:, 1:realRank);
            end
        end
        
    end
    
    methods (Static)
        
        function I = fromIsoDec(iso)
        % Constructs an irreducible decomposition from an isotypic decomposition
            sample = iso.algebra.sampleGeneric;
            U = iso.U;
            repTypes = zeros(1, iso.nComponents);
            for r = 1:iso.nComponents
                d = iso.repDims(r);
                m = iso.repMuls(r);
                range = iso.compRange(r);
                if iso.repIsReal(r)
                    % use a second sample to put all irreducible components in the same basis
                    if iso.ordered
                        Urep = U(:, range);
                    else
                        Urep = iso.refinedBasis(r);
                    end
                    % only need the first row of blocks
                    repSample = Urep(:, 1:d)'*(sample + sample')*Urep;
                    % constructing the change of basis matrices
                    P = cell(1, m);
                    P{1} = eye(d);
                    for j = 2:m
                        P{j} = repSample(:, (j-1)*d+(1:d))';
                        P{j} = P{j} * sqrt(d/trace(P{j}*P{j}'));
                    end
                    Urep = Urep * blkdiag(P{:});
                    U(:, range) = Urep;
                    repTypes(r) = 1;
                else
                    Urep = replab.rep.IrrDec.orderedComplexBasis(iso, r);
                    if isequal(Urep, [])
                        % Quaternionic type
                        U(:, range) = NaN;
                        repTypes(r) = 3;
                    else
                        % Complex type
                        repSample = Urep(:, 1:d)'*sample*Urep;
                        P = cell(1, m);
                        P{1} = eye(d);
                        for j = 2:m
                            P{j} = repSample(:, (j-1)*d+(1:d))';
                            P{j} = P{j} * sqrt(d/trace(P{j}*P{j}'));
                        end
                        Urep = Urep * blkdiag(P{:});
                        U(:, range) = Urep;
                        repTypes(r) = 2;
                    end
                end
            end
            I = replab.rep.IrrDec(iso.algebra, iso.fromFiber, U, iso.repDims, iso.repMuls, repTypes);
        end
                
    end
    
end

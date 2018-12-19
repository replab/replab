classdef PhaseConfiguration
    
    properties
        n;          % matrix size
        phase;      % phase array, containing roots of unity
                    %
                    % The o-th orbit has indices between position
                    % orbitStart(o) and orbitStart(o+1)-1
                    % in the flat arrays orbitRow, orbitCol and orbitInd.
                    % orbitInd corresponds to Matlab linear indexing of matrices.
                    %
        orbitStart; % start position of orbit indices for the o-th orbit (int32)
        orbitRow;   % row indices of orbits, flat array (int32)
        orbitCol;   % col indices of orbits, flat array (int32)
        orbitInd;   % flat indices of orbits, flat array (int32)
        index;      % index of the cell for each element, or 0 if element is zero
    end
    
    methods
        
        function self = PhaseConfiguration(n, phase, orbitStart, orbitRow, orbitCol, index)
            self.n = double(n);
            self.phase = double(phase);
            self.orbitStart = int32(orbitStart(:)');
            self.orbitRow = int32(orbitRow(:)');
            self.orbitCol = int32(orbitCol(:)');
            self.orbitInd = self.orbitRow + n * (self.orbitCol - 1);
            self.index = double(index);
        end
        
        function pc = restrict(self, block)
        % Takes a slice of this phase configuration for the
        % indices in block; this is used for example to restrict
        % the representation decomposition algorithm to the orbits
        % of a permutation group
            newN = length(block);
            newPhase = self.phase(block, block);
            oldIndex = self.index(block, block);
            oldOrbits = unique(self.index(:));
            newIndex = zeros(newN, newN);
            newOrbitStart = [];
            newOrbitRow = [];
            newOrbitCol = [];
            newO = 1;
            for oldO = 1:length(oldOrbits)
                [R C] = find(oldIndex == oldOrbits(oldO));
                if ~isempty(R)
                    newOrbitStart = [newOrbitStart length(newOrbitRow)+1];
                    newOrbitRow = [newOrbitRow R(:)'];
                    newOrbitCol = [newOrbitCol C(:)'];
                    newIndex(oldIndex == oldOrbits(newO)) = newO;
                    newO = newO + 1;
                end
            end
            newOrbitStart = [newOrbitStart length(newOrbitRow)+1];
            pc = replab.rep.PhaseConfiguration(newN, newPhase, ...
                                               newOrbitStart, ...
                                               newOrbitRow, ...
                                               newOrbitCol, ...
                                               newIndex);
        end

        function m = nOrbits(self)
        % Returns the number of orbits
            m = length(self.orbitStart) - 1;
        end

        function m = orbitSize(self, o)
        % Returns the size of the o-th orbit
            m = double(self.orbitStart(o + 1) - self.orbitStart(o));
        end

        function t = orbitTranspose(self, o)
        % Returns the index of the transpose of this orbit, i.e. the orbit
        % corresponding to {(j,i) for all (i,j) in the o-th orbit}
            i = self.orbitStart(o);
            t = self.index(self.orbitCol(i), self.orbitRow(i));
        end

        function t = orbitTransposeScalar(self, o, s)
        % Returns the scalar t corresponding to the transpose of the o-th orbit
        % when the o-th orbit is assigned the scalar s
        %
        % Enables 
        % setOrbit(M, o, s)
        % setOrbit(M, orbitTranspose(o), orbitTransposeScalar(o, s))
            i = self.orbitStart(o);
            r = self.orbitRow(i);
            c = self.orbitCol(i);
            t = s / self.phase(r, c) * self.phase(c, r);
        end

        function M1 = project(self, M)
            n = self.n;
            M1 = zeros(n, n);
            for o = 1:self.nOrbits
                M1 = self.setOrbit(M1, o, self.getOrbitAverage(M, o));
            end
        end

        function s = getOrbitAverage(self, M, o)
            linIndices = self.orbitInd(self.orbitStart(o):self.orbitStart(o+1)-1);
            s = dot(M(linIndices), conj(self.phase(linIndices))/self.orbitSize(o));
        end

        function M = setOrbit(self, M, o, s)
        % Sets all elements of the o-th orbit in the matrix M to the value s
            linIndices = self.orbitInd(self.orbitStart(o):self.orbitStart(o+1)-1);
            M(linIndices) = self.phase(linIndices) * s;
        end

        function b = isOrbitDiagonal(self, o)
            i = self.orbitStart(o);
            b = (self.orbitRow(i) == self.orbitCol(i)); % if the first element is diagonal, all elements are diagonal 
        end

        function b = isOrbitSelfTranspose(self, o)
        % Returns whether the o-th orbit contains its transpose, i.e.
        % for every (i,j) in the orbit, (j,i) is also contained
            b = self.orbitTranspose(o) == o;
        end

        function M = sampleRealGaussian(self)
        % Samples from a matrix with real entries distributed
        % according to the standard normal distribution, and
        % then symmetrized under the current phase configuration
        %
        % Similar to Random.realGaussian followed by self.project
        %
        % Standard deviation of average of n random normal
        % variables is origStdDev/sqrt(n)
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                s = randn / sqrt(self.orbitSize(o));
                M = self.setOrbit(M, o, s);
            end
        end

        function M = sampleComplexGaussian(self)
        % Similar to sampleRealGaussian except it samples
        % from complex normal entries
        %
        % Standard deviation of real/imag part is 1/sqrt(2)
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                s = (randn + randn*1i) / sqrt(self.orbitSize(o) * 2);
                M = self.setOrbit(M, o, s);
            end
        end

        function M = sampleSymmetricGaussian(self)
        % Works as if it samples from a matrix from the 
        % Gaussian Orthogonal Ensemble
        % and then symmetrizes it under the current phase configuration
        %
        % Similar to Random.symmetricGaussian followed by self.project
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                if self.isOrbitDiagonal(o)
                    % diagonal elements are scaled up by sqrt(2)
                    s = randn * sqrt(2/self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                elseif self.isOrbitSelfTranspose(o)
                    % elements are averaged over standard normals
                    s = randn / sqrt(self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                else
                    ot = self.orbitTranspose(o);
                    if o < ot
                        % the orbit index o is minimal under transpose, so do it
                        % factor 2 due to the transpose part
                        s = randn / sqrt(self.orbitSize(o)*2);
                        st = self.orbitTransposeScalar(o, s);
                        M = self.setOrbit(M, o, s);
                        M = self.setOrbit(M, ot, st);
                    end % else do nothing, will be handled by the orbit transpose
                end
            end
        end

        function M = sampleHermitianGaussian(self)
        % Similar to sampleSymmetricGaussian for matrices
        % from the Gaussian Unitary Ensemble
        %
        % Similar to Random.hermitianGaussian followed by self.project
            n = self.n;
            M = zeros(n, n);
            for o = 1:self.nOrbits
                if self.isOrbitDiagonal(o)
                    % diagonal elements are real
                    s = randn / sqrt(self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                elseif self.isOrbitSelfTranspose(o)
                    % self transpose so real
                    % off-diagonal elements are scaled down by sqrt(2)
                    s = randn / sqrt(2*self.orbitSize(o));
                    M = self.setOrbit(M, o, s);
                else
                    ot = self.orbitTranspose(o);
                    if o < ot
                        % factor 2 for off-diagonal, and factor 2 for
                        % averaging over the orbit transpose
                        s = (randn + randn*1i) / sqrt(self.orbitSize(o)*4);
                        st = self.orbitTransposeScalar(o, s);
                        M = self.setOrbit(M, o, s);
                        M = self.setOrbit(M, ot, conj(st));
                    end
                end
            end
        end


    end
    
    methods (Static)

        function phase = computePhase(shift, order)
            if order == 2
                phase = 1 - 2*shift;
            else
                phase = zeros(n, n);
                for r = 1:n
                    for c = 1:n
                        phase(r, c) = replab.rep.PhaseConfiguration.rootOfUnity(shift(r, c), order);
                    end
                end
            end
        end
        
        function C = fromSignedPermMatlab(generators)
            n = size(generators, 2);
            B = replab.rep.PhaseConfigurationBuilder.fromSignedPerm(generators);
            phase = replab.rep.PhaseConfiguration.computePhase(B.shift, B.order);
            C = replab.rep.PhaseConfiguration(n, phase, B.orbitStart, B.orbitRow, B.orbitCol, B.index);
        end
        
        function C = fromSignedPermJava(generators)
            n = size(generators, 2);
            assert(n > 0, 'Java code does not work with empty array, as Matlab will pass a null pointer.');
            % convert 1-based indices to 0-based indices
            generators(generators > 0) = generators(generators > 0) - 1;
            B = com.faacets.qdimsum.PhaseConfigurationBuilder.fromSignedPerm(n, generators);
            phase = replab.rep.PhaseConfiguration.computePhase(double(B.shift), double(B.order));
            C = replab.rep.PhaseConfiguration(n, phase, B.orbitStart + 1, B.orbitRow + 1, B.orbitCol + 1, B.index + 1);
        end
        
        function C = fromSignedPerm(generators)
            if exist('com.faacets.qdimsum.PhaseConfigurationBuilder')
                if size(generators, 1) > 0
                    C = replab.rep.PhaseConfiguration.fromSignedPermJava(generators);
                else
                    C = replab.rep.PhaseConfiguration.fromSignedPermMatlab(generators);
                end
            else
                warning('Warning: qdimsum.jar not in the Java path, using Matlab code as fallback.');
                C = replab.rep.PhaseConfiguration.fromSignedPermMatlab(generators);
            end
        end
        
        function u = rootOfUnity(k, n)
            switch n
              case 1
                u = 1;
              case 2
                switch k
                  case 0
                    u = 1;
                  case 1
                    u = -1;
                end
              case 4
                switch k
                  case 0
                    u = 1;
                  case 1
                    u = 1i;
                  case 2
                    u = -1;
                  case 3
                    u = -1i;
                end
              otherwise
                u = exp(2*pi*1i*k/n);
            end
        end
    end


end

classdef Matrices < replab.cat.Domain
    
    properties (SetAccess = protected)
        nRows;
        nCols;
        complex;
        selfAdjoint;
        canEqv = false;
        canHash = false;
        canSample = true;
    end
    
    methods
        
        function self = Matrices(nRows, nCols, complex, selfAdjoint)
            if self.nRows ~= self.nCols
                assert(~selfAdjoint);
            end
            self.nRows = nRows;
            self.nCols = nCols;
            self.complex = complex;
            self.selfAdjoint = selfAdjoint;
        end
        
        function s = str(self)
            s = sprintf('%d x %d', self.nRows, self.nCols);
            if self.complex
                if self.selfAdjoint
                    s = [s ' hermitian '];
                else
                    s = [s ' complex '];
                end
            else
                if self.selfAdjoint
                    s = [s ' symmetric '];
                else
                    s = [s ' real '];
                end
            end
            s = [s 'double matrices'];
        end
        
        function b = eqv(self, x, y)
            error('Cannot test floating point matrices for equality');
        end
        
        function h = hash(self, x)
            error('Cannot hash floating point matrices');
        end
        
        function M = sample(self)
            if self.selfAdjoint
                n = self.nRows;
                if self.complex
                    % Generates a Hermitian matrix with measure invariant
                    % under unitary transformations,
                    % sampled from the Gaussian Unitary Ensemble, see
                    % http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
                    M = zeros(n, n);
                    for r = 1:n
                        M(r, r) = randn;
                        M(r, r+1:end) = (randn(1, n-r) + randn(1, n-r)*1i)/sqrt(2);
                        M(r+1:end, r) = conj(M(r, r+1:end));
                    end
                else
                    % Generates a symmetric matrix with measure invariant 
                    % under orthogonal transformations,
                    % sampled from the Gaussian Orthogonal Ensemble, see
                    % http://staff.math.su.se/shapiro/UIUC/random_matrices.pdf
                    M = zeros(n, n);
                    for r = 1:n
                        % diagonal elements are scaled up by sqrt(2)
                        M(r, r) = randn * sqrt(2);
                        % while other elements are standard normals
                        M(r, r+1:end) = randn(1, n-r);
                        M(r+1:end, r) = M(r, r+1:end);
                    end
                end
            else % not self adjoint
                nRows = self.nRows;
                nCols = self.nCols;
                if self.complex
            M = (randn(nRows, nCols) + randn(nRows, nCols) * 1i)/sqrt(2);
                else
                    M = randn(nRows, nCols);
                end
                    
            end
        end
        
    end
    
end

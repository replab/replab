classdef CommutantConfiguration < replab.Commutant

    properties (SetAccess = protected)
        field; % vector space field, 'R' (real) or 'C' (complex)
        n;     % matrix size
        group; % group
        rep;   % representation
        phaseConfiguration;
    end
    
    properties (Access = protected)
        parent_; % parent domain: real or complex matrices
    end
    
    methods
        
        function self = CommutantConfiguration(rep, phaseConfiguration)
            self.rep = rep;
            self.n = rep.dimension;
            self.field = rep.field;
            self.group = rep.group;
            self.phaseConfiguration = phaseConfiguration;
            switch self.field
              case 'R'
                self.parent_ = replab.domain.RealMatrices(self.n, self.n);
              case 'C'
                self.parent_ = replab.domain.ComplexMatrices(self.n, self.n);
              case 'H'
                self.parent_ = replab.domain.QuaternionMatrices(self.n, self.n);
              otherwise
                error('Unknown field');
            end
        end
        
        function X = project(self, X)
        % Projects any n x n matrix in the equivariant subspace
            X = self.phaseConfiguration.project(X);
        end

        % Str
        
        function s = headerStr(self)
            s = sprintf('%d x %d %s commutant matrices', ...
                        self.n, self.n, ...
                        replab.str.field(self.field));
        end
        
        % Domain
        
        function b = eqv(self, X, Y)
            b = self.parent_.eqv(X, Y);
        end
        
        function X = sample(self)
            X = self.project(self.parent_.sample);
        end
        
        function X = sampleHermitian(self)
            X = self.parent_.sample;
            X = self.project(X);
            X = (X + X')/2;
        end
        
    end
    
end

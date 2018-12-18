classdef Settings
    
    methods (Static)
        
        function value = bsgsFailureProbability(newValue)
            persistent BsgsFailureProbability;
            if nargin == 1
                BsgsFailureProbability = newValue;
            elseif isequal(BsgsFailureProbability, [])
                BsgsFailureProbability = 2^-100;
            end
            value = BsgsFailureProbability;
        end
        
        function value = matrixNormTol(m, n, field)
            
        end
        
        function value = vectorNormTol(d, field)
        end
        
        function value = eigTol(field)
            switch field
              case {'R7', 'C7'}
                value = sqrt(replab.doubleEigTol)
              case {'R15', 'C15'}
                value = replab.doubleEigTol;
              otherwise
                error('Not implemented');
            end
        end
        
        function value = doubleEigTol(newValue)
            persistent DoubleEigTol;
            if nargin == 1
                DoubleEigTol = newValue;
            elseif isequal(DoubleEigTol, [])
                DoubleEigTol = 1e-9;
            end
            value = DoubleEigTol;
        end
            
    end
    
end

classdef IrrepInfo < replab.Str

    properties (SetAccess = protected)
        label % char or []: Short description of the irreducible representation, or [] if no description is known
              %             By convention, the string '1' represents the trivial representation
        divisionAlgebra % {'R', 'C', 'H', []} : For representations over the reals, the division algebra type if known
                        %                       If it is unknown, or the subrepresentation is complex, must be [].
        isDivisionAlgebraCanonical % logical: Whether the division algebra is known to be in the canonical form,
                                   %          as in `replab.domain.ComplexTypeMatrices` and
                                   %          `replab.domain.QuaternionTypeMatrices`. False is means either not true
                                   %          or status unknown.
    end
   
    methods
        
        function self = IrrepInfo(label, divisionAlgebra, isDivisionAlgebraCanonical)
            if nargin < 2
                divisionAlgebra = [];
                isDivisionAlgebraCanonical = [];
            end
            self.label = label;
            self.divisionAlgebra = divisionAlgebra;
            self.isDivisionAlgebraCanonical = isDivisionAlgebraCanonical;
        end
        
        %% Str methods
        
        function s = shortStr(self, maxColumns)
            s = {};
            if ~isempty(self.label)
                s{end+1} = ['''' self.label ''''];
            end
            if ~isempty(self.divisionAlgebra)
                switch self.divisionAlgebra
                  case 'R'
                    s{end+1} = 'real of real type';
                  case 'C'
                    s{end+1} = 'real of complex type';
                  case 'H'
                    s{end+1} = 'real of quaternionic type';
                end
            end
            if isequal(self.isDivisionAlgebraCanonical, true)
                s{end+1} = 'canonical basis';
            end
            s = strjoin(s, ', ');
        end
        
    end
    
end

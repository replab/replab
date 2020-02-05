classdef Info < replab.Str

    properties (SetAccess = protected)
        divisionAlgebra % ({'R', 'C', 'H', []}): For representations over the reals, the division algebra type if known
                        %
                        %                        If the division algebra is unknown, the value is ``[]``. If the
                        %                        representation is complex, the value is ``[]`` as well.
        isDivisionAlgebraCanonical % (logical or []): Whether the division algebra is known to be in the canonical form
                                   %
                                   %                  The canonical form is described in `+replab.+domain.ComplexTypeMatrices` and
                                   %                  `+replab.+domain.QuaternionTypeMatrices`. The value must be []
                                   %                  if the representation is complex, or if the division algebra is ``R``.
    end

    methods

        function self = Info(divisionAlgebra, isDivisionAlgebraCanonical)
            self.divisionAlgebra = divisionAlgebra;
            if isequal(self.divisionAlgebra, 'R')
                assert(isequal(self.isDivisionAlgebraCanonical, []));
            end
            self.isDivisionAlgebraCanonical = isDivisionAlgebraCanonical;
        end

        %% Str methods

        function s = shortStr(self, maxColumns)
            s = {'irreducible'};
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

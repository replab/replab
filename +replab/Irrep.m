classdef Irrep < replab.SubRep
% Describes a subrepresentation of a unitary finite representation
%
% 
    properties (SetAccess = protected)
        realDivisionAlgebra;
    end
    
    methods
        
        function self = Irrep(parent, U0, realDivisionAlgebra)
            self = self@replab.SubRep(parent, U0);
            if isequal(self.field, 'R')
                assert(isa(realDivisionAlgebra, 'replab.DivisionAlgebra'));
            else
                assert(isequal(realDivisionAlgebra, []));
            end
            self.realDivisionAlgebra = realDivisionAlgebra;
        end
        
        function s = headerStr(self)
            if isequal(self.realDivisionAlgebra, [])
                s = 'Irreducible subrepresentation';
            else
                switch self.realDivisionAlgebra.shortName
                  case 'R'
                    s = 'Real-type real irreducible subrepresentation';
                  case 'C'
                    s = 'Complex-type real irreducible subrepresentation';
                  case 'H'
                    s = 'Quaternion-type real irreducible subrepresentation';
                end
            end
        end

        
        function sub1 = recoverRational(self)
            sub1 = self;
            % TODO: handle various real division algebras
        end
        
    end

end

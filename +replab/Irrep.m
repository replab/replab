classdef Irrep < replab.SubRep
% Describes an irreducible subrepresentation
%
% The additional constraints compared to replab.SubRep are:
% - this subrepresentation is irreducible over its field,
% - if the field is R, then the type (real, complex, quaternion)
%   of this irreducible representation is known, and the
%   encoding of complex/quaternion types is canonical
    properties (SetAccess = protected)
        realDivisionAlgebra; % = [] for representation over 'C'
                             % or a replab.DivisionAlgebra
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
            if self.overC || self.realDivisionAlgebra.isReal
                U1 = replab.rep.recoverRational(self);
                if isequal(U1, [])
                    sub1 = self;
                else
                    sub1 = replab.Irrep(self.parent, U1, self.realDivisionAlgebra);
                end
            else
                sub1 = self;
            end
        end
        
    end

end

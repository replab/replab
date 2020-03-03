classdef HarmonizedIsotypic < replab.Isotypic
% Describes an isotypic component in the decomposition of a representation
%
% Improves on `.Isotypic` by having all irreducible representations in the same basis.

    methods

        function self = HarmonizedIsotypic(parent, irreps)
            self = self@replab.Isotypic(parent, irreps);
        end

        %% Rep methods

% $$$         function rho = image(self, g)
% $$$             rho = image@replab.SubRep(self, g);
% $$$             return
% $$$             p = self.parent.image(g);
% $$$             U = self.irrep(1).U;
% $$$             rho = U*p*U';
% $$$             for i = 2:self.nIrreps
% $$$                 U = self.irrep(i).U;
% $$$                 rho = rho + U*p*U';
% $$$             end
% $$$             rho = rho / self.nIrreps;
% $$$             rho = kron(eye(self.nIrreps), rho);
% $$$         end
% $$$
% $$$         function c = commutant(self)
% $$$             c = commutant@replab.SubRep(self);
% $$$             return
% $$$             if isempty(self.commutant_)
% $$$                 if self.overC
% $$$                     self.commutant_ = replab.IsotypicSimpleCommutant(self);
% $$$                 else
% $$$                     switch self.irrep(1).irrepInfo.divisionAlgebra
% $$$                       case 'R'
% $$$                         self.commutant_ = replab.IsotypicSimpleCommutant(self);
% $$$                       case 'C'
% $$$                         self.commutant_ = replab.IsotypicComplexCommutant(self);
% $$$                       case 'H'
% $$$                         self.commutant_ = replab.IsotypicQuaternionCommutant(self);
% $$$                     end
% $$$                 end
% $$$             end
% $$$             c = self.commutant_;
% $$$         end

    end

end

classdef HarmonizedIsotypic < replab.Isotypic
% Describes an isotypic component in the decomposition of a representation
%
% Improves on `.Isotypic` by having all irreducible representations in the same basis.

    methods

        function self = HarmonizedIsotypic(parent, irreps, E_internal)
            self = self@replab.Isotypic(parent, irreps, E_internal);
        end

        %% Rep methods

        function rho = image_internal(self, g)
            p = self.parent.image_internal(g);
            E = self.irrep(1).E_internal;
            B = self.irrep(1).B_internal;
            rho = E*p*B;
            for i = 2:self.nIrreps
                E = self.irrep(i).E_internal;
                B = self.irrep(i).B_internal;
                rho = rho + E*p*B;
            end
            rho = rho / self.nIrreps;
            rho = kron(eye(self.nIrreps), rho);
        end

        function c = commutant(self)
            if isempty(self.commutant_)
                if self.overC
                    self.commutant_ = replab.IsotypicSimpleCommutant(self);
                else
                    switch self.irrep(1).frobeniusSchurIndicator
                      case 1
                        self.commutant_ = replab.IsotypicSimpleCommutant(self);
                      case 0
                        self.commutant_ = replab.IsotypicComplexCommutant(self);
                      case -1
                        self.commutant_ = replab.IsotypicQuaternionCommutant(self);
                      otherwise
                        error('Unknown indicator');
                    end
                end
            end
            c = self.commutant_;
        end

        %% Isotypic methods

        function [A Ainv] = changeOfBasis(self, i, j, context)
            A = eye(self.irrepDimension);
            Ainv = A;
        end

    end

end

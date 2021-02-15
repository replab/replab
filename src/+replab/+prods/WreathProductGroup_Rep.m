classdef WreathProductGroup_Rep < replab.Rep
% Imprimitive or primitive representation of a wreath product group

    properties (SetAccess = protected)
        factorRep % (`+replab.Rep`): Representation of the group whose copies are acted upon
        type % (``'primitive'`` or ``'imprimitive'``): Representation type
    end

    properties (Access = protected)
        baseRep % (`+replab.Rep`): Representation of the base group
        actingRep % (`+replab.Rep`): Representation of the acting group
    end

    methods

        function self = WreathProductGroup_Rep(group, type, factorRep, actingRep, baseRep)
        % Constructs an imprimitive representation of a wreath product group
        %
        % Args:
        %   group (`+replab.WreathProductGroup`): Wreath product group
        %   (``'primitive'`` or ``'imprimitive'``): Representation type
        %   factorRep (`+replab.Rep`): Representation of the base factor group
        %   actingRep (`+replab.Rep`): Permutation representation of the acting group
        %   baseRep (`+replab.Rep`): Representation of the base group
            assert(isa(group, 'replab.WreathProductGroup'));
            assert(ismember(type, {'primitive', 'imprimitive'}));
            assert(strcmp(factorRep.field, actingRep.field));
            assert(strcmp(factorRep.field, baseRep.field));
            assert(isa(factorRep, 'replab.Rep') && factorRep.group == group.A);
            assert(isa(actingRep, 'replab.Rep') && actingRep.group == group.H);
            assert(actingRep.isExact && actingRep.isUnitary);
            assert(actingRep.dimension == baseRep.dimension);
            assert(isa(baseRep, 'replab.Rep') && baseRep.group == group.N);
            self@replab.Rep(group, baseRep.field, baseRep.dimension, 'isUnitary', baseRep.isUnitary);
            self.type = type;
            self.factorRep = factorRep;
            self.actingRep = actingRep;
            self.baseRep = baseRep;
        end

    end

    methods (Access = protected) % Implementations

        % Rep

        function rho = image_exact(self, g)
            rho = self.actingRep.image(g{1}, 'exact') * self.baseRep.image(g{2}, 'exact');
        end

        function rho = image_double_sparse(self, g)
            rho = self.actingRep.image(g{1}, 'double/sparse') * self.baseRep.image(g{2}, 'double/sparse');
        end

        function e = computeErrorBound(self)
            e = self.baseRep.errorBound;
        end

        function c = computeConditionNumberEstimate(self)
            c = self.baseRep.conditionNumberEstimate;
        end

        function b = computeIsUnitary(self)
            b = self.dimension <= 1 || self.factorRep.isUnitary;
        end

    end

    methods % Implementations

        % Rep

        function b = isExact(self)
            b = self.factorRep.isExact;
        end

    end

end

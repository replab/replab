classdef SemigroupFun < replab.DomainFun & replab.Semigroup
    properties (SetAccess = protected)
        composeFun;
    end
    methods
        function self = SemigroupFun(description, eqvFun, sampleFun, ... % Domain
                                     composeFun) % Semigroup
            self = self@replab.DomainFun(description, eqvFun, sampleFun);
            self.composeFun = composeFun;
        end
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.DomainFun(self), ...
                hiddenFields@replab.Semigroup(self), ...
                {'composeFun'} ...
                );
        end
    end
end

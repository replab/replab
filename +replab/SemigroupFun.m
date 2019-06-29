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
            names1 = hiddenFields@replab.DomainFun(self);
            names2 = hiddenFields@replab.Semigroup(self);
            names = vertcat(names1(:), names2(:));
            names{end+1, 1} = 'composeFun';
            names = unique(names);
        end
    end    
end

classdef MonoidFun < replab.SemigroupFun & replab.Monoid
    methods
        function self = MonoidFun(description, eqvFun, sampleFun, ... % Domain
                                     composeFun, ... % Semigroup
                                     identity) % Monoid
            self = self@replab.SemigroupFun(description, eqvFun, sampleFun, composeFun);
            self.identity = identity;
        end
        function names = hiddenFields(self)
            names1 = hiddenFields@replab.SemigroupFun(self);
            names2 = hiddenFields@replab.Monoid(self);
            names = vertcat(names1(:), names2(:));
            names{end+1, 1} = 'inverseFun';
            names = unique(names);
        end
    end
end

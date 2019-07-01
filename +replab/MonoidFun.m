classdef MonoidFun < replab.SemigroupFun & replab.Monoid
    methods
        function self = MonoidFun(description, eqvFun, sampleFun, ... % Domain
                                     composeFun, ... % Semigroup
                                     identity) % Monoid
            self = self@replab.SemigroupFun(description, eqvFun, sampleFun, composeFun);
            self.identity = identity;
        end
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.SemigroupFun(self), ...
                hiddenFields@replab.Monoid(self), ...
                {'inverseFun'} ...
                );
        end
    end
end

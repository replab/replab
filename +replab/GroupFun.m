classdef GroupFun < replab.MonoidFun & replab.Group
    properties (SetAccess = protected)
        inverseFun;
    end
    methods
        function self = GroupFun(description, eqvFun, sampleFun, ... % Domain
                                 composeFun, ... % Semigroup
                                 identity, ... % Monoid
                                 inverseFun) % Group
            self = self@replab.MonoidFun(description, eqvFun, sampleFun, composeFun, identity);
            self.inverseFun = inverseFun;
        end
        function names = hiddenFields(self)
            names1 = hiddenFields@replab.MonoidFun(self);
            names2 = hiddenFields@replab.Group(self);
            names = horczat(names1(:), names2(:));
            names{end+1, 1} = 'inverseFun';
            names = unique(names);
        end
    end    
end

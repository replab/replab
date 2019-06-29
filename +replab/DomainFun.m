classdef DomainFun < replab.Domain & replab.StrFun
    properties (SetAccess = protected)
        eqvFun;
        sampleFun;
    end
    methods
        function self = DomainFun(description, eqvFun, sampleFun)
            self = self@replab.StrFun(description, description);
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
        end
        
        function names = hiddenFields(self)
            names1 = hiddenFields@replab.Domain(self);
            names2 = hiddenFields@replab.StrFun(self);
            names = vertcat(names1(:), names2(:));
            names{end+1, 1} = 'eqvFun';
            names{end+1, 1} = 'sampleFun';
            names = unique(names);
        end

    end
end

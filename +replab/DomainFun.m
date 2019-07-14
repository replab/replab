classdef DomainFun < replab.Domain & replab.StrFun
    properties (SetAccess = protected)
        eqvFun;
        sampleFun;
    end
    methods
        function self = DomainFun(description, eqvFun, sampleFun)
            self@replab.StrFun(description, description);
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
        end
        
        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Domain(self), ...
                hiddenFields@replab.StrFun(self), ...
                {'eqvFun' 'sampleFun'} ...
                );
        end

    end
end

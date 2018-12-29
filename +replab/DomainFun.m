classdef DomainFun < replab.Domain & replab.Str
    properties (SetAccess = protected)
        eqvFun;
        sampleFun;
    end
    methods
        function self = DomainFun(description, eqvFun, sampleFun)
            self = self@replab.Str(description);
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
        end
    end
end

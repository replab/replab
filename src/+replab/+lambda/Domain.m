classdef Domain < replab.Domain
% An implementation of a domain defined by image functions

    properties (SetAccess = protected)
        header
        eqvFun
        sampleFun
    end

    methods

        function self = Domain(header, eqvFun, sampleFun)
            self.header = header;
            self.eqvFun = eqvFun;
            self.sampleFun = sampleFun;
        end

        function str = headerStr(self)
            str = self.header;
        end

        function names = hiddenFields(self)
            names = replab.str.uniqueNames( ...
                hiddenFields@replab.Domain(self), ...
                {'header'} ...
                );
        end

        % Domain methods

        function b = eqv(self, t, u)
            f = self.eqvFun;
            b = f(t, u);
        end

        function t = sample(self)
            f = self.sampleFun;
            t = f();
        end

    end

end

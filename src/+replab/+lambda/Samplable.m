classdef Samplable < replab.Samplable
% A function_handle based implementation of `+replab.Samplable`

    properties (SetAccess = protected)
        header
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

        % Samplable methods

        function t = sample(self)
            f = self.sampleFun;
            t = f();
        end

    end

end

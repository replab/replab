classdef DispatchNext < replab.Str

    properties (SetAccess = protected)
        message % (charstring): Explanation of the failure of the particular implementation
    end

    methods

        function self = DispatchNext(message)
            if nargin >= 1
                self.message = message;
            else
                self.message = '';
            end
        end

        function s = shortStr(self)
            if isempty(self.message)
                s = 'DispatchNext';
            else
                s = sprintf('DispatchNext: %s', self.message);
            end
        end

    end

end

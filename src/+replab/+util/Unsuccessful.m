classdef Unsuccessful < replab.Str
% Unsuccessful event
%
% This class implements a value signifying that an attempt to do something was
% unsuccessful. It is not an error to be unsuccessful, but an allowed behavior.

    properties (SetAccess = protected)
        message % (charstring): Explanation of the reason why a particular procedure was not successful
    end

    methods

        function self = Unsuccessful(message)
            if nargin >= 1
                self.message = message;
            else
                self.message = '';
            end
        end

        function s = shortStr(self)
            if isempty(self.message)
                s = 'Unsuccessful';
            else
                s = sprintf('Unsuccessful: %s', self.message);
            end
        end

    end

end

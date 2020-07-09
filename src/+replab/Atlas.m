classdef Atlas < replab.Str
% An atlas of finite groups

    properties
        maximalOrder % (integer): Maximal order for which to attempt group recognition
    end

    methods

        function self = Atlas(maximalOrder)
            self.maximalOrder = maximalOrder;
        end

        function setAsDefault(self)
        % Sets this atlas as the default atlas
        %
        % Modifies the global variable `+replab.+globals.atlas`
            replab.globals.atlas(self);
        end

        function result = recognize(self, group)
        % Attempts to identify the given group
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case of positive identification; or ``[]`` if unrecognized.
            error('Abstract');
        end

    end

end

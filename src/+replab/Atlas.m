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
            replab.globals.atlas(self);
        end

        function result = recognize(self, group)
        % Attemps to identify the given group
            error('Abstract');
        end

    end

end

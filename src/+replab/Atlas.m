classdef Atlas
% An atlas of finite groups

    methods (Static)

        function result = recognize(self, group)
        % Attempts to identify the given group
        %
        % Returns:
        %   `+replab.AtlasResult` or []: A result in case of positive identification; or ``[]`` if unrecognized.
            error('Abstract');
        end

    end

end

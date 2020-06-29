classdef Morphism < replab.Str

    properties (SetAccess = protected)
        source % (`.Group`): Source group
        target % (`.Group`): Target group
    end

    methods

        function t = image(s)
            error('Abstract');
        end

    end

end

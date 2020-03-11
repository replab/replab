classdef Root < replab.lobster.Node

    methods

        function str = render(self, context)
            if nargin < 2
                context = struct();
            end

            str = self.render_children(context);
        end

    end

end

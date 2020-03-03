classdef IfNode < replab.lobster.Node

    properties (Constant)
        SUPPORTED_OPERATORS = {'<' '>' '<=' '>=' '==' '~='}
    end

    properties
        expression = ''

        if_branch = cell(0)
        else_branch = cell(0)
    end

    methods
        function self = IfNode(fragment)
            process_fragment(self, fragment);
            self.creates_scope = true;
        end

        function process_fragment(self, fragment)
            self.expression = fragment;
        end

        function str = render(self, context)

            % Problem: the IF block needs to have the others as children
            if replab.lobster.eval_with_context(self.expression, context)
                str = self.render_children(context, self.if_branch);
            else
                str = self.render_children(context, self.else_branch);
            end
        end

        function enter_scope(~)

            % We don't have to do anything here.
        end

        function exit_scope(self)

            % We now have enough information to figure out what to print if the
            % expression evaluates to true or false.

            on_if_branch = true;
            for k = 1:length(self.children)
                if isa(self.children{k}, 'replab.lobster.ElseNode')
                    on_if_branch = false;
                end

                if on_if_branch;
                    self.if_branch = {self.if_branch{:}, self.children{k}};
                else
                    self.else_branch = {self.else_branch{:}, self.children{k}};
                end
            end

        end

    end

end

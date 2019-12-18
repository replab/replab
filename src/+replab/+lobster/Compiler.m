classdef Compiler < handle
   
    properties (Constant)
        VAR_TOKEN_START   = '{{'
        VAR_TOKEN_END     = '}}'
        BLOCK_TOKEN_START = '{%'
        BLOCK_TOKEN_END   = '%}'
        
        TOKEN_REGEX_PATTERN = sprintf('(%s.*?%s|%s.*?%s)', replab.lobster.Compiler.VAR_TOKEN_START, replab.lobster.Compiler.VAR_TOKEN_END, replab.lobster.Compiler.BLOCK_TOKEN_START, replab.lobster.Compiler.BLOCK_TOKEN_END)
    end
    
    properties
       template_string
    end
    
    methods (Access = private)
        
        function fragments = make_fragments(self)
            % Seperate the text sections and the variables/blocks
            [vars, text] = regexp(self.template_string, self.TOKEN_REGEX_PATTERN, 'match', 'split');
            
            % Now we need to use the vars and text to rebuild the original
            % string as a stream of tokens. 
            %
            % To explain how we achieve this consider the following inputs:
            %
            %   >> [vars, text] = regexp('{{ var }} Hello World!', ...
            % vars = 
            %     '{{ var }}'
            % text = 
            %     ''    ' Hello World!'
            %
            % >> [vars, text] = regexp('Hello World! {{ var }}'', ...)
            % vars = 
            %     '{{ var }}'
            % text = 
            %     'Hello World! '    ''
            %
            % >> [vars, text] = regexp('{{ var }}{{ var }}', ...)
            % vars = 
            %     '{{ var }}'    '{{ var }}'
            % text = 
            %     ''    ''    ''
            %
            % To rebuild the string we just have to alternate tokens from 'vars'
            % and 'text'. We will simply have to remove the empty strings from
            % the token string later on.
            
            tokens = cell(0);
            for k = 1:length(vars)
                tokens{k*2} = vars{k};
            end
            for k = 1:length(text)
                tokens{(k-1)*2 + 1} = text{k};
            end
            
            tokens = tokens(cellfun(@(x) ~isempty(x), tokens));
            fragments = cellfun(@(x) replab.lobster.Fragment(x), tokens, 'Uniform', false);
            
            % Post process the fragments to remove text fragments that only
            % contain a newline after a block:
            discard = zeros(1, length(fragments));
            for k = 2:length(fragments)
                if fragments{k-1}.type ~= replab.lobster.FRAGMENT_TYPE.TEXT && ...
                   fragments{k}.type   == replab.lobster.FRAGMENT_TYPE.TEXT && ...
                   strcmp(fragments{k}.raw, sprintf('\n'))
                    discard(k) = 1;
                end
            end
            
            fragments = fragments(~discard);
        end
        
        function new_node = create_node(~, fragment)
           if fragment.type == replab.lobster.FRAGMENT_TYPE.TEXT
               new_node = replab.lobster.TextNode(fragment.clean);
           elseif fragment.type == replab.lobster.FRAGMENT_TYPE.VAR
               new_node = replab.lobster.VarNode(fragment.clean);
           elseif fragment.type == replab.lobster.FRAGMENT_TYPE.BLOCK_START
               block_type = strsplit(fragment.clean);
               
               if strcmp(block_type(1), 'if')
                   new_node = replab.lobster.IfNode(strjoin(block_type(2:end), ' '));
               elseif strcmp(block_type(1), 'else')
                   new_node = replab.lobster.ElseNode();
               elseif strcmp(block_type(1), 'for')
                   new_node = replab.lobster.ForNode(strjoin(block_type(2:end), ' '));
                elseif strcmp(block_type(1), 'call')
                   new_node = replab.lobster.CallNode(strjoin(block_type(2:end), ' '));
               else
                   error('Lobster:TemplateSytaxError', ...
                       '<%s> is not a valid block type', block_type{1});
               end
           else
               error('Lobster:TemplateSyntaxError', ...
                   '<%s> looks like invalid syntax', fragment.raw);
           end
        end
    end
    
    methods
        
        function self = Compiler(template_string)
            self.template_string = template_string;
        end
        
        function root = compile(self)
            root = replab.lobster.Root();
            scope_stack = {root};
            
            fragments = self.make_fragments();
            for k = 1:length(fragments)
                fragment = fragments{k};
                
                if isempty(scope_stack)
                    error('Lobster:NestingError', ...
                        'It look like you have a nesting issue in your template');
                end
                
                parent_scope = scope_stack{end};
                
                % If we are exting a block we need to move down the scope stack
                if fragment.type == replab.lobster.FRAGMENT_TYPE.BLOCK_END
                    parent_scope.exit_scope();
                    scope_stack(end) = [];

                    continue
                end
                
                new_node = self.create_node(fragment);
                parent_scope.add_child(new_node);
                
                if new_node.creates_scope
                    scope_stack{end+1} = new_node;
                    new_node.enter_scope();
                end
            end
        end
    end
    
end

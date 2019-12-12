classdef DocTest < replab.Str
    
    properties
        commands % row cell vector of charstring: commands to be evaluated
        outputs % row cell vector of charstring: expected output
    end

    methods
        
    end
    
    methods (Static)
        
        function dt = parseDocTest(doclines)
            
        end
        
        
        function [command output] = parseCommand(ps)
        end
        
        function [commands outputs] = parse(ps)
            
        end
        
        function dt = parseDoc(doc)
            doctests = {};
            
            for i = 1:doc.nLines
                l = doc.line(i);
                tokens = regexp(l, '^(\s*)(.*)', 'tokens', 'once');
                indent = length(tokens{1});
                content = strtrim(tokens{2});
                if ~isempty(content)
                    if isequal(content, 'Example:')
                        
                    end
                end
                
            end
            
        end
        
        
    end
end
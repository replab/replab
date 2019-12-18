classdef Fragment < handle
    
    properties 
        raw = ''; 
        clean = '';
        type;
    end
    
    methods (Access = private)
        
        function cleaned_text = clean_fragment(~, raw)
        % Strip the token start and end tags from the raw text of the framgent and remove whitespace. 
            token_starts = {replab.lobster.Compiler.VAR_TOKEN_START replab.lobster.Compiler.BLOCK_TOKEN_START};
            
            if length(raw) >= 2 && ~isempty(find(strcmp(raw(1:2), token_starts)))
                cleaned_text = strtrim(raw(3:end-2));
            else 
                cleaned_text = strrep(raw, char(10), '');
                cleaned_text = strrep(cleaned_text, '\n', char(10));
                cleaned_text = strrep(cleaned_text, '\t', char(9));
                cleaned_text = strrep(cleaned_text, '\\', '\');
            end
        end
        
        function compute_type(self)

        % If the length of the raw string is less than the minimum length
        % needed to fit the START_BLOCK and END_BLOCK delimiters then we can
        % be sure that the fragment is a text node.  
            if length(self.raw) < 4
                self.type = replab.lobster.FRAGMENT_TYPE.TEXT;
                return
            end
            
            raw_start = self.raw(1:2);
            if strcmp(raw_start, replab.lobster.Compiler.VAR_TOKEN_START)
                self.type = replab.lobster.FRAGMENT_TYPE.VAR;
            elseif strcmp(raw_start, replab.lobster.Compiler.BLOCK_TOKEN_START)
                if strcmp(self.clean(1:3), 'end')
                    self.type = replab.lobster.FRAGMENT_TYPE.BLOCK_END;
                else
                    self.type = replab.lobster.FRAGMENT_TYPE.BLOCK_START;
                end               
            else
                self.type = replab.lobster.FRAGMENT_TYPE.TEXT;
            end
        end
    end
    
    methods 
        
        function self = Fragment(raw)
            self.raw = raw;
            self.clean = self.clean_fragment(raw);
            self.compute_type();
        end
        
    end
    
end

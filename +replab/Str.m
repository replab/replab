classdef Str < handle
% Defines a 'str' default method and overloads 'disp'
    
    methods
        
        function disp(self)
            disp(self.str);
        end
        
        function s = str(self, short)
            try
                s = self.description;
            catch
                s = sprintf('%s instance', class(self));
            end
        end

    end
    
end

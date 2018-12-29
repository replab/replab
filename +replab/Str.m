classdef Str < handle
% Defines a 'str' default method and overloads 'disp'
%
% Also provides a 'description' property used as default 'str' output
    
    properties (SetAccess = protected)
        description = [];
    end
    
    methods
            
        function self = Str(description)
            if nargin < 1
                self.description = sprintf('%s instance', class(self));
            else
                self.description = description;
            end
        end
        
        function disp(self)
            disp(self.str);
        end
        
        function s = str(self)
            s = self.description;
        end

    end
    
end

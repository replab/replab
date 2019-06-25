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
            disp(strjoin(replab.longStr(self), '\n'));
        end
        
        function s = str(self)
            s = self.description;
        end
        
        function [names values] = additionalFields(self)
        % Returns the name/value pairs corresponding to additional
        % fields to be printed
            names = {};
            values = {};
        end
        
        function [s overLimit] = shortStr(self, maxColumns)
        % Returns a single line description of the current object
        % see replab.str.shortStr for documentation
            [s overLimit] = replab.str.shortStr(self, maxColumns);
        end
        
        function [s overLimit] = longStr(self, maxRows, maxColumns)
        % Returns a multi line description of the current object
        % see replab.str.longStr for documentation
            [s overLimit] = replab.str.longStr(self, maxColumns);
        end

    end
    
end

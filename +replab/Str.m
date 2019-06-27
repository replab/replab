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
        
        function [names values] = additionalFields(self)
        % Returns the name/value pairs corresponding to additional
        % fields to be printed; given as column vectors
        % Classes that override this method should call the
        % superclass method
            names = {};
            values = {};
        end
        
        function names = hiddenFields(self)
        % Returns the names of the fields that are not printed
        % Classes that override this method should call the
        % superclass method
            names = {};
        end
        
        function s = shortStr(self, maxColumns)
        % Returns a single line description of the current object
        % see replab.shortStr for documentation
            s = replab.str.shortStr(self, maxColumns);
        end
        
        function s = longStr(self, maxRows, maxColumns)
        % Returns a multi line description of the current object
        % see replab.longStr for documentation
            s = replab.str.longStr(self, maxRows, maxColumns);
        end

    end
    
end

classdef Str < handle
% Defines a 'str' default method and overloads 'disp'
%
% Also provides methods 'additionalFields' and 'hiddenFields' to
% guide long form object pretty printg
    
    methods
        
        function res = eq(self, rhs)
        % Workaround bug of == not implemented for handles
            if replab.platformIsOctave
                res = true(size(self));
            else
                res = eq@handle(self, rhs);
            end
        end
        
        function disp(self)
            maxRows = replab.Parameters.strMaxRows;
            maxColumns = replab.Parameters.strMaxColumns;
            lines = replab.longStr(self, maxRows, maxColumns);
            lines = replab.str.longFit(lines, maxRows, maxColumns);
            disp(strjoin(lines, '\n'));
        end
        
        function [names values] = additionalFields(self)
        % Returns the name/value pairs corresponding to additional
        % fields to be printed; given as row vectors
        % Classes that override this method should call the
        % superclass method
            names = {};
            values = {};
        end
        
        function names = hiddenFields(self)
        % Returns the names of the fields that are not printed as a row vector
        % Classes that override this method should call the superclass method
            names = {};
        end
        
        function s = shortStr(self, maxColumns)
        % Returns a single line description of the current object
        % see replab.shortStr for documentation
            s = replab.str.shortStr(self, maxColumns);
        end
        
        function s = headerStr(self, maxColumns)
        % Returns a tiny single line description of the current object type
        % see replab.headerStr for documentation
            s = self.shortStr(maxColumns);
        end
        
        function s = longStr(self, maxRows, maxColumns)
        % Returns a multi line description of the current object
        % see replab.longStr for documentation
            s = replab.str.longStr(self, maxRows, maxColumns);
        end

    end
    
end

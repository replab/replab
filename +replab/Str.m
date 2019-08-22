classdef Str < handle
% Provides infrastructure for pretty printing
%
% This class overloads the :meth:`disp` method.
%
% To customize pretty-printing, subclasses can either:
%
% * override the :meth:`tinyStr`, :meth:`shortStr`, :meth:`longStr` methods,
%
% * override the methods :meth:`additionalFields` and :meth:`hiddenFiels` to 
%   control the manner objects are pretty printed.
    
    methods
                
        function [names values] = additionalFields(self)
        % Adds name/value pairs to the printed object properties
        %
        % Classes that override this method should call the
        % superclass method and concatenate the results
        %
        % Returns
        % -------
        %   names: row cell array of strings
        %     Names of the keys of additional properties
        %   values: row cell array of objects
        %     Values of those additional properties
        %
            names = {};
            values = {};
        end
        
        function disp(self)
        % Prints the object on the REPL (overriden standard Matlab method)
            maxRows = replab.Settings.strMaxRows;
            maxColumns = replab.Settings.strMaxColumns;
            lines = replab.longStr(self, maxRows, maxColumns);
            lines = replab.str.longFit(lines, maxRows, maxColumns);
            disp(strjoin(lines, '\n'));
        end

        function s = headerStr(self, maxColumns)
        % Returns a tiny single line description of the current object type
        %
        % See :func:`+replab.headerStr` for documentation
            s = self.shortStr(maxColumns);
        end
        
        function names = hiddenFields(self)
        % Returns the names of the fields that are not printed as a row vector
        %
        % Classes that override this method should call the superclass method
        % and concatenate the results
        %
        % Returns:
        %   row cell array of strings: Names of the hidden keys
            names = {};
        end
        
        function s = longStr(self, maxRows, maxColumns)
        % Returns a multi line description of the current object
        % 
        % See :func:`+replab.longStr` for documentation
            s = replab.str.longStr(self, maxRows, maxColumns);
        end

        function s = shortStr(self, maxColumns)
        % Returns a single line description of the current object
        %
        % See :func:`+replab.shortStr` for documentation
            s = replab.str.shortStr(self, maxColumns);
        end
        
    end
    
end

classdef Str < handle
% Defines a 'str' default method and overloads 'disp'
%
% Also provides methods `.additionalFields` and `.hiddenFields` to guide long form object pretty printing

    methods

        function res = eq(self, rhs)
        % Equality test
        %
        % Workaround bug of == not implemented for handles
        %
        % Args:
        %   self (object): first object
        %   rhs (object): second object to compare to
        %
        % Returns:
        %   boolean: true iff self == rhs
            if replab.compat.isOctave
                res = true(size(self));
            else
                res = eq@handle(self, rhs);
            end
        end

        function disp(self)
            maxRows = replab.settings.strMaxRows;
            maxColumns = replab.settings.strMaxColumns;
            lines = replab.longStr(self, maxRows, maxColumns);
            lines = replab.str.longFit(lines, maxRows, maxColumns);
            disp(strjoin(lines, '\n'));
        end

        function [names, values] = additionalFields(self)
        % Returns the name/value pairs corresponding to additional fields to be printed
        %
        % Classes that override this method should call the superclass method.
        %
        % Returns
        % -------
        %   names: row cell vector of charstring
        %     Additional names to display
        %   values: row cell vector of charstring
        %     Additional values to display
            names = {};
            values = {};
        end

        function names = hiddenFields(self)
        % Returns the names of the fields that are not printed as a row vector
        %
        % Classes that override this method should call the superclass method.
        %
        % Returns:
        %   row cell vector of charstring: Field names to hide
            names = {};
        end

        function s = shortStr(self, maxColumns)
        % Single line text description of the current object
        %
        % See also:
        %   `+replab.shortStr`
            s = replab.str.shortStr(self, maxColumns);
        end

        function s = headerStr(self)
        % Tiny single line description of the current object type
        %
        % See also:
        %   `+replab.headerStr`
            s = replab.str.headerStr(self);
        end

        function s = longStr(self, maxRows, maxColumns)
        % Multi-line description of the current object
        %
        % See also:
        %   `+replab.longStr`
            s = replab.str.longStr(self, maxRows, maxColumns);
        end

    end

end

classdef Str < handle
% Base class that provides sane pretty printing for all RepLAB objects
%
% All classes in RepLAB inherit the `.Str` base class, which provides explicit `~+replab.Str.longStr` and `~+replab.Str.shortStr` methods.
% Those methods take a dimension limit for width (and possibly height).
% In contrast, the functions `+replab.shortStr` and `+replab.longStr` use default values that are
% deduced from the terminal size.
%
% This base class overloads 'disp'
%
% It also provides methods `.additionalFields` and `.hiddenFields` to guide long form object pretty printing.
%
% Compare the two outputs:
%
% Example:
%   >>> P = replab.SymmetricGroup(3)
%     P =
%     Symmetric group acting on 3 elements
%       domainSize: 3
%         identity: [1, 2, 3]
%             type: Symmetric group acting on 3 elements
%     generator(1): [2, 3, 1]
%     generator(2): [2, 1, 3]
%   >>> replab.longStr(P)
%     ans =
%     6x1 cell array
%     {'Symmetric group acting on 3 elements'              }
%     {'  domainSize: 3                                   '}
%     {'    identity: [1, 2, 3]                           '}
%     {'        type: Symmetric group acting on 3 elements'}
%     {'generator(1): [2, 3, 1]                           '}
%     {'generator(2): [2, 1, 3]                           '}
%   >>> replab.shortStr(P)
%     ans =
%     'Symmetric group acting on 3 elements'

    methods

        function disp(self)
        % Standard MATLAB/Octave display method
            maxRows = replab.globals.strMaxRows;
            maxColumns = replab.globals.strMaxColumns;
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
        % The default implementation of 'longStr' is to print a short description of
        % the object on the first line, followed by public properties.
        %
        % See also:
        %   `+replab.longStr`
            s = replab.str.longStr(self, maxRows, maxColumns);
        end

    end

end

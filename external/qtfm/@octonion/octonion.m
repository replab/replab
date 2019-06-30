% Quaternion Toolbox for matlab (QTFM)
%
% Octonion class definition and constructor method/function. The octonion
% class is built on top of the quaternion class and unlike the quaternion
% class, needs a class definition to establish that it is superior to the
% quaternion class. This means that an expression such as q + o which
% attempts to add a quaternion to an octonion (not possible of course) will
% always call the octonion function plus.m. Therefore any error handling
% may be included only in the octonion code and we do not have to rewrite
% existing quaternion code. This also fits with the primary purpose of QTFM
% which is to be a quaternion toolbox - the octonion class is a bonus, but
% likely to remain somewhat experimental for some years to come.

% Copyright (c) 2011 Stephen J. Sangwine and Nicolas Le Bihan.
% See the file : Copyright.m for further details.

classdef (InferiorClasses = {?quaternion}) octonion
    properties (SetAccess = 'private', Hidden = false)
        a = quaternion();
        b = quaternion();
    end
    methods
        function o = octonion(a0, a1, a2, a3, a4, a5, a6, a7)
        % OCTONION   Construct octonions from components. Accepts the following
        % possible arguments, which may be scalars, vectors or matrices:
        %
        % No argument     - returns an empty octonion scalar, vector or matrix.
        % One argument    - An octonion argument returns the argument unmodified.
        %                   A non-octonion argument returns the argument in the
        %                   scalar part and supplies a zero vector part with
        %                   elements of the same type as the scalar part.
        % Two arguments   - returns an octonion, provided the first argument is
        %                   numeric and the second is a pure octonion.
        % Seven arguments - returns a pure octonion scalar, vector or matrix,
        %                   with an empty scalar part.
        % Eight arguments - returns a full octonion scalar, vector or matrix.

        nargoutchk(0, 1) % The number of input arguments is checked
        % as part of the switch statement below.
        switch nargin
            case 0

                % Construct an empty octonion with two empty quaternion components.

                o.a = quaternion(); o.b = o.a; %o = class(o, 'octonion');

            case 1

                if isa(a0, 'octonion')
                    o = a0; % a0 is already an octonion, so return it.
                else
                    % a0 is not a octonion ...
                    if isnumeric(a0)
                        % ... but it is numeric. We can handle this, but we have to
                        % supply the vector components of the octonion, and they must
                        % have the same type as a0 (e.g. double, int8).

                        o.a = quaternion(a0); o.b = zerosq(size(a0), class(a0));
                        %o = class(o, 'octonion');
                    else
                        % ... or it isn't numeric.
                        error(['Cannot construct an octonion with a ',...
                            'non-numeric in the scalar part.']);
                    end
                end

            case 2

                % In this case, the first argument must be the scalar part, and thus
                % numeric, and the second must be the vector part, and thus a pure
                % octonion (with the same component class as the first argument).

                if any(size(a0) ~= size(a1))
                    error('Arguments must have the same dimensions')
                end

                if isnumeric(a0) && isa(a1, 'octonion')
                    c0 = class(a0); c1 = class(x(a1.a));
                    if isempty(s(a1.a))
                        if strcmp(c0, c1)
                            % a1 is already an octonion, so copy it to the output and
                            % insert a0 as the scalar component.
                            o = a1; o.a = quaternion(a0, v(o.a));
                        else
                            error(['Classes of given scalar and vector parts', ...
                                ' differ: ', c0, ' ', c1])
                        end
                    else
                        error('The second argument must be a pure octonion.');
                    end
                else
                    error(['First argument must be numeric and the ',...
                        'second must be a pure octonion.']);
                end

            case 7

                % Construct a pure octonion, since we are given seven arguments and
                % this is the only possibility. They must all be numeric and of the same
                % class.

                s0 = size(a0);
                if any(s0 ~= size(a1)) || any(s0 ~= size(a2)) || any(s0 ~= size(a3)) ...
                        || any(s0 ~= size(a4)) || any(s0 ~= size(a5)) || any(s0 ~= size(a6))
                    error('Arguments must have the same dimensions')
                end

                c0 = class(a0); c1 = class(a1); c2 = class(a2); c3 = class(a3);
                c4 = class(a4); c5 = class(a5); c6 = class(a6);
                % We test only the first argument for numeric, because if the classes
                % match, the others must also be numeric.
                if isnumeric(a0) && strcmp(c0, c1) && strcmp(c0, c2) && strcmp(c0, c3) ...
                        && strcmp(c0, c4) && strcmp(c0, c5) && strcmp(c0, c6)
                    o.a = quaternion(    a0, a1, a2);
                    o.b = quaternion(a3, a4, a5, a6);
                    %o = class(o, 'octonion');
                else
                    error(['All seven arguments must be numeric and of the',  ...
                        ' same class. Given: ', c0, ' ', c1, ' ', c2, ' ', ...
                        c3, ' ', c4, ' ', c5, ' ', c6]);
                end

            case 8 % Return a full octonion. All eight arguments must be numeric and
                % of the same class.

                s0 = size(a0);
                if any(s0 ~= size(a1)) || any(s0 ~= size(a2)) || any(s0 ~= size(a3)) ...
                        || any(s0 ~= size(a4)) || any(s0 ~= size(a5)) || any(s0 ~= size(a6)) ...
                        || any(s0 ~= size(a7))
                    error('Arguments must have the same dimensions')
                end

                c0 = class(a0); c1 = class(a1); c2 = class(a2); c3 = class(a3);
                c4 = class(a4); c5 = class(a5); c6 = class(a6); c7 = class(a7);
                % We test only the first argument for numeric, because if the classes
                % match, the other three must also be numeric.
                if isnumeric(a0) && strcmp(c0, c1) && strcmp(c0, c2) && strcmp(c0, c3) ...
                        && strcmp(c0, c4) && strcmp(c0, c5) && strcmp(c0, c6)
                    o.a = quaternion(a0, a1, a2, a3);
                    o.b = quaternion(a4, a5, a6, a7);
                    % o = class(o, 'octonion');
                else
                    error(['All eight arguments must be numeric and of the',   ...
                        ' same class. Given: ', c0, ' ', c1, ' ', c2, ' ', ...
                        c3, ' ', c4, ' ', c5, ' ', ...
                        c6, ' ', c7]);
                end

            otherwise
                error('Octonion constructor takes 0, 1, 2, 7 or 8 arguments.');
        end
        end
        function n = numArgumentsFromSubscript(~,~,~)
            % Introduction of this function with Matlab R2015b and QTFM 2.4
            % permitted numel to revert to its obvious function of
            % providing the number of elements in an array.
            n = 1;
        end
    end
    methods (Static = true)
        function q = empty(varargin)
            % This function makes it possible to use the dotted notation
            % octonion.empty or octonion.empty(0,1) to create an empty
            % array.
            % TODO Make it possible to write quaternion.empty('double') and
            % specify the class of the empty components.
            d = double.empty(varargin{:});
            q = octonion(d, d, d, d, d, d, d, d);
        end
    end
end

% $Id: octonion.m 1004 2017-11-15 17:14:09Z sangwine $

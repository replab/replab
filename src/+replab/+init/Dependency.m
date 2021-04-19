classdef Dependency < replab.Str
% Describes a RepLAB dependency

    properties (SetAccess = protected)
        name % (charstring): Name of the dependency
    end

    methods % Dependency methods (Abstract)

        function res = inPath(self)
        % Returns true if the library is in the path already
            error('Abstract');
        end

        function res = works(self)
        % Returns true if the library is working (as checked by solving a toy problem)
            error('Abstract');
        end

        function init(self, folderName)
        % Initalizes the library (e.g. by adding it to the path)
        %
        % Args:
        %   folderName (charstring): Base path in which the library is located
            error('Abstract');
        end

    end

    methods % Methods assuming a particular directory structure (external/)

        function require(self)
        % Locates, installs and initializes the library
        %
        % The automatic installation step is only performed is the `+replab.+globals.autoInstall` flag is set
            error('Abstract');
        end

    end

end

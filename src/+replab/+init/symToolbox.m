classdef symToolbox < replab.init.Dependency
% Verifies that symbolic computation is available
%
% This loads the ``symbolic`` Octave package if necessary

    methods

        function self = symToolbox
            self.name = 'symToolbox';
        end

        function res = inPath(self)
            res = any(exist('syms') == [2 6]);
        end

        function init(self, path)
            if replab.compat.isOctave
                try
                    replab.init.log_(1, 'Loading symbolic package for Octave...');
                    pkg load symbolic
                catch
                end
            end
        end

        function res = works(self)
            evalc('sym(''2'');'); % run a first command to trigger loading of the library
            res = double(sym('2') + sym('2')) == 4;
        end

        function require(self)
            if self.inPath
                assert(self.works, 'The symbolic toolbox is in the path but does not work properly.');
            else
                self.init;
                assert(self.inPath, 'Initialization of the symbolic toolbox failed');
                assert(self.works, 'Symbolic toolbox initialized but not working');
            end
        end

    end

end

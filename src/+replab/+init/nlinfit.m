classdef nlinfit < replab.init.Dependency
% Verifies that the nlinfit optimization function is available
%
% This loads the ``optim`` Octave package if necessary

    methods

        function self = nlinfit
            self.name = 'nlinfit';
        end

        function res = inPath(self)
            res = any(exist('nlinfit') == [2 6]);
        end

        function res = works(self)
            res = self.inPath;
        end

        function init(self, folderName)
            if replab.compat.isOctave
                replab.init.log_(1, 'Loading optim package for Octave...');
                pkg load optim
            end
        end

        function require(self)
            if self.inPath
                assert(self.works, 'The nlinfit function is in the path but does not work properly.');
            else
                self.init;
                assert(self.inPath, 'Initialization of nlinfit failed');
                assert(self.works, 'nlinfit initialized but not working');
            end
        end

    end

end

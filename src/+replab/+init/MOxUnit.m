classdef MOxUnit < replab.init.ExternalDependency

    methods

        function self = MOxUnit
            self@replab.init.ExternalDependency('MOxUnit', 'MOxUnit/moxunit_set_path.m');
        end

        function res = works(self)
            res = false;
            try
                assertEqual(2 + 2, 4);
                res = true;
            catch
            end
        end

        function res = inPath(self)
            res = any(exist('moxunit_util_platform_is_octave') == 2);
        end

        function init(self, path)
            run(fullfile(path, 'MOxUnit', 'moxunit_set_path.m'));
        end

    end

end

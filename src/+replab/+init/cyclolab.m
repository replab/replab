classdef cyclolab < replab.init.ExternalDependency
% Adds the cyclolab library to the Java path if it is not yet present

    properties (Constant = true)
            files = {'algebra_2.11-1.0.0.jar' 'cats-kernel_2.11-1.0.1.jar' 'cyclo-core_2.11-0.16.0.1-SNAPSHOT.jar' 'cyclolab_2.11-0.2-SNAPSHOT.jar' 'fastparse_2.11-2.1.2.jar' 'machinist_2.11-0.6.4.jar' 'scala-library-2.11.12.jar' 'scala-reflect-2.11.12.jar' 'scalin-core_2.11-0.16.0.2-SNAPSHOT.jar' 'scalin-macros_2.11-0.16.0.2-SNAPSHOT.jar' 'sourcecode_2.11-0.1.6.jar' 'spire_2.11-0.16.0.jar' 'spire-macros_2.11-0.16.0.jar' 'spire-platform_2.11-0.16.0.jar' 'spire-util_2.11-0.16.0.jar'};
    end

    methods

        function self = cyclolab
            self@replab.init.ExternalDependency('cyclolab', 'cyclolab_2.11-0.2-SNAPSHOT.jar');
        end

        function res = inPath(self)
            res = false;
            if replab.compat.isOctave
                if ~eval('__have_feature__(''JAVA'')')
                    return
                end
            end
            try
                array = javaMethod('parse', 'cyclo.Lab', {'1'});
                res = true;
            catch
            end
        end

        function res = works(self)
            res = self.inPath;
        end

        function init(self, path)
            res = false;
            try
                b = javaMethod('valueOf', 'java.math.BigInteger', 2);
                res = true;
            catch
            end
            if ~res
                error('The Java virtual machine is not available');
            end
            files = replab.init.cyclolab.files;
            for i = 1:length(files)
                f = files{i};
                if ~exist(fullfile(path, f)) == 2
                    error('File %s is missing\n', f);
                end
            end
            for i = 1:length(files)
                f = files{i};
                javaaddpath(fullfile(path, f));
            end
        end

    end

end

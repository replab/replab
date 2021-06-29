function res = initCycloLab(verbose)
% Adds the gluon library to the Java path if it is not yet present
%
% Args:
%   verbose ({0, 1, 2}): Controls the display level
%
% Returns:
%   logical: True if the library is available
    t = cputime;
    basePath = replab.globals.replabPath;
    res = false;
    try
        b = javaMethod('valueOf', 'java.math.BigInteger', 2);
    catch
        if verbose > 0
            disp('The Java virtual machine is not available');
        end
        return
    end
    files = {'algebra_2.11-1.0.0.jar' ...
             'cats-kernel_2.11-1.0.1.jar' ...
             'cyclo-core_2.11-0.16.0.1-SNAPSHOT.jar' ...
             'cyclolab_2.11-0.2-SNAPSHOT.jar' ...
             'fastparse_2.11-2.1.2.jar' ...
             'machinist_2.11-0.6.4.jar' ...
             'scala-library-2.11.12.jar' ...
             'scala-reflect-2.11.12.jar' ...
             'scalin-core_2.11-0.16.0.2-SNAPSHOT.jar' ...
             'scalin-macros_2.11-0.16.0.2-SNAPSHOT.jar' ...
             'sourcecode_2.11-0.1.6.jar' ...
             'spire_2.11-0.16.0.jar' ...
             'spire-macros_2.11-0.16.0.jar' ...
             'spire-platform_2.11-0.16.0.jar' ...
             'spire-util_2.11-0.16.0.jar' ...
            };
    good = true;
    for i = 1:length(files)
        f = files{i};
        if ~exist([basePath '/external/cyclolab/' f]) == 2
            if verbose > 0
                fprintf('File %s is missing\n', f);
            end
            good = false;
        end
    end
    if ~good
        return
    end
    for i = 1:length(files)
        f = files{i};
        javaaddpath([basePath '/external/cyclolab/' f]);
    end
    try
        array = javaMethod('parse', 'cyclo.Lab', {'1' '2' '3'});
        disp('Added CycloLab to the path')
        res = true;
    catch
        if verbose > 0
            disp('Error while testing CycloLab');
        end
    end
end

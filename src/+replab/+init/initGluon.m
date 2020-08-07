function res = initGluon(verbose)
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
    if ~exist('java.math.BigInteger')
        if verbose > 0
            disp('The Java virtual machine is not available');
        end
        return
    end
    files = {'algebra_2.11-1.0.0.jar' ...
             'cats-kernel_2.11-1.0.1.jar' ...
             'cyclo-core_2.11-0.16.0.1-SNAPSHOT.jar' ...
             'gluon_2.11-0.2-SNAPSHOT.jar' ...
             'machinist_2.11-0.6.4.jar' ...
             'scala-compiler-2.11.12.jar' ...
             'scala-library-2.11.12.jar' ...
             'scala-parser-combinators_2.11-1.0.4.jar' ...
             'scala-reflect-2.11.12.jar' ...
             'scala-xml_2.11-1.0.5.jar' ...
             'scalin-core_2.11-0.16.0.2-SNAPSHOT.jar' ...
             'scalin-macros_2.11-0.16.0.2-SNAPSHOT.jar' ...
             'spire_2.11-0.16.0.jar' ...
             'spire-macros_2.11-0.16.0.jar' ...
             'spire-platform_2.11-0.16.0.jar' ...
             'spire-util_2.11-0.16.0.jar' ...
             'fastparse_2.11-2.1.2.jar' ...
             'sourcecode_2.11-0.1.6.jar'
            };
    good = true;
    for i = 1:length(files)
        f = files{i};
        if ~exist([basePath '/external/gluon/' f]) == 2
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
        javaaddpath([basePath '/external/gluon/' f]);
    end
    try
        if verbose > 0
            disp('Gluon: compiling...');
        end
        i = com.faacets.gluon.Interface.compile(com.faacets.gluon.Interface.squareCode);
        if verbose > 0
            disp('Gluon: checking...');
        end
        if i.call(2) ~= 4
            error('Error: wrong result from interface');
        end
        if verbose > 0
            fprintf('Gluon: loaded in %.2f CPU seconds.\n', cputime-t);
        end
        res = true;
    catch
        if verbose > 0
            disp('Error while testing Gluon interface');
        end
    end
end

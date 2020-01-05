classdef PackageData < replab.Str

    properties
        path % row cell vector of charstring: path of the package
        ownFunctions % row cell vector of `+replab.+infra.FunctionLikeData`: Data of functions
        ownClasses % row cell vector of `+replab.+infra.ClassData`: Data of classes
    end

    methods

        function self = PackageData(path, ownFunctions, ownClasses)
            self.path = path;
            self.ownFunctions = ownFunctions;
            self.ownClasses = ownClasses;
        end

    end

end

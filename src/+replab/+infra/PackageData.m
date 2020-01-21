classdef PackageData < replab.Str

    properties
        path % row cell vector of charstring: path of the package
        ownFunctions % row cell vector of `+replab.+infra.FunctionLikeData`: Data of functions
        ownClasses % row cell vector of `+replab.+infra.ClassData`: Data of classes
    end

    methods

        function self = PackageData(path, ownFunctions, ownClasses)
            self.path = path;
            assert(all(cellfun(@(x) all(x ~= '_'), path)), 'Package names cannot contain _');
            % if the restriction above needs to be lifted, the package name mangling needs to
            % be revised, in crawl.m and in the shm code
            self.ownFunctions = ownFunctions;
            self.ownClasses = ownClasses;
        end

    end

end

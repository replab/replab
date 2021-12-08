function ok = replab_checkoop
% Checks the object oriented programming structure of the code
%
% The main reason for this function is the fact that inheritance rules are
% a bit different between Matlab and Octave.
%
% This includes three kinds of checks, performed one after the other one
% upon success:
%
% 1. Check that inherited methods admit only one possible implementation.
% When inheriting from multiple classes, a given method can be provided
% several definitions. Since Matlab and Octave use different rules to
% choose which definition to use (and some of these rules might evolve in
% the future), we check here that no such situation arises. This is needed
% to guarantee that the same code is provided to both interpreters.
%
% 2. When the first test passes, this function also checks that all classes
% which are not inherited by any other class contain no 'Abstract' method,
% i.e. that no top class is abstract.
%
% 3. Finally, we check that all references to replab objects in the code
% point to an existing definition (this is similar to what replab_checkhelp
% does for the documentation).
%
% Returns:
%   logical: True if no problems were detected
    ok = true;
    help('--clear');

    rp = replab.globals.replabPath;
    srcRoot = fullfile(rp, 'src');


    disp('Crawling code base');
    cb = replab.infra.crawl(srcRoot);
    
    disp('Checking for potentially multiple method definitions');
    ac = cb.allClasses;
    % Sort out methods by inheritance number
    [~, I] = sort(cellfun(@(cl) length(cl.allSuperclasses), ac));
    ac = ac(I);
    nbAmbiguousMethods = 0;
    for i = 1:length(ac)
        cl = ac{i};
        if length(cl.ownSuperclasses) >= 2
            disp([replab.infra.repl.linkOpen(cl.fullIdentifier, cl.fullIdentifier, cl.absoluteFilename, 1), ' inherits from ', num2str(length(cl.ownSuperclasses)), ' other classes.']);
            % This class inherits from several classes. We check that all
            % its methods are uniquely defined
            mts = cl.allMethods;
            for j = 1:length(mts)
                decl = mts{j}.declarations.findAll;
                if isa(mts{j}, 'replab.infra.InheritedClassElement') && (length(decl) > 1)
                    % There are several definitions, we check whether a
                    % single definition comes on top of all other ones
                    
                    % Starting from the main class, we explore all
                    % superclasses until finding the actually admissible
                    % definitions
                    nbImplementationsFound = 0;
                    implementations = [];
                    classesToExplore = cl.ownSuperclasses;
                    while ~isempty(classesToExplore)
                        cle = classesToExplore{1};
                        meth = cle.lookup(mts{j}.name);
                        if isempty(meth)
                            % method not present in this branch, nothing to
                            % do here
                        elseif isa(meth, 'replab.infra.InheritedClassElement')
                            % The method comes from further down the
                            % hierarchy, we need to explore further
                            classesToExplore = [classesToExplore, cle.ownSuperclasses];
                        elseif isa(meth, 'replab.infra.ConcreteClassElement')
                            % The method is declared here, we check which
                            % of the initial declarations this corresponds
                            % to
                            nbImplementationsFound = nbImplementationsFound + 1;
                            for k = 1:length(decl)
                                if isequal(meth, decl{k})
                                    implementations(nbImplementationsFound) = k;
                                end
                            end
                        else
                            error('Unknown method type');
                        end
                        classesToExplore = classesToExplore(2:end);
                    end
                    
                    % We check that all reached implementation are
                    % identical
                    assert(min(implementations) > 0, 'Found an implementation that was not registred in the codebase');
                    if min(implementations) ~= max(implementations)
                        nbAmbiguousMethods = nbAmbiguousMethods + 1;
                        disp(['Method ', mts{j}.fullIdentifier, ' admits ', num2str(length(unique(implementations))), ' declarations:']);
                        for k = unique(implementations)
                            disp(['  - in ', replab.infra.repl.linkOpen(decl{k}.fullIdentifier, decl{k}.fullIdentifier, decl{k}.absoluteFilename, decl{k}.startLineNumber)]);
                        end
                        disp(' ');
                    end
                end
            end
        end
    end

    ok = (nbAmbiguousMethods == 0);
    if ok
        disp('No ambiguous method overloading encountered.')
    else
        warning([num2str(nbAmbiguousMethods), ' methods were found who could be ambiguous.']);
        return;
    end
    
    
    
end

function replab_release
% replab_release Release the current develop branch contents as a stable release
%
% replab_release runs the release process to send the current ``develop``
% branch snapshot to the branch ``master``, taking care of version numbers in the process.
%
% There are two types of version numbers in RepLAB:
%
% * Version numbers ending in ``-SNAP`` are snapshot version numbers which are not
%   supposed to be stable: that means different versions of the RepLAB codebase can
%   have the same snapshot number.
% * Version numbers **not** ending in ``-SNAP`` are stable version numbers, which are
%   in one to one correspondence with a version of the codebase. It may correspond at
%   most to two commits: the stable commit on the ``develop`` branch, and the merge commit
%   of it on ``master``.
%
% The master branch has only strictly increasing, stable version numbers. The develop
% branch has snapshot version numbers, except during the release process when we "stabilize"
% the version number for one commit.
%
% The release process is as follows.
%
% 0. (Outside of the script) The user must run ``git fetch origin master develop``. 
%
% 1. We verify that the repository does not have uncommited changes.
%
% 2. We verify that both the ``develop`` and ``master`` branches are in sync with 
%    the ``origin`` remote. If not we abort.
%
% 2. We display the current snapshot version number, which is going to become master, the new snapshot version 
% 
% Move to the right folder
    
    assert(exist('replab_version.txt') == 2, 'The current directory must be the RepLAB root folder.');
    [status, cmdout] = system('git rev-parse --abbrev-ref HEAD');
    assert(isequal(strtrim(cmdout), 'develop'), 'The current repository must be on the develop branch to continue.');
    
    input('Step 0: Press ENTER to confirm that you ran "git fetch origin develop master"');
    
    disp('Step 1: Verifying that the current working tree and index are clean');
    [status, cmdout] = system('git diff --exit-code');
    assert(status == 0, 'The repository has local unstaged changes.');
    [status, cmdout] = system('git diff --cached --exit-code');
    assert(status == 0, 'The repository has staged but uncommitted changes.');

    disp('Step 2: Verifying that master and develop branches are in sync with remote origin.');    
    [status, masterRef] = system('git rev-parse master');
    assert(status == 0, 'Git command failed');
    [status, originMasterRef] = system('git rev-parse origin/master');
    assert(status == 0, 'Git command failed');
    [status, developRef] = system('git rev-parse develop');
    assert(status == 0, 'Git command failed');
    [status, originDevelopRef] = system('git rev-parse origin/develop');
    assert(status == 0, 'Git command failed');
    assert(isequal(masterRef, originMasterRef), 'Please synchronize master with origin/master.');
    assert(isequal(developRef, originDevelopRef), 'Please synchronize develop with origin/develop.');
    
    
% $$$     initialPath = pwd;
% $$$     [pathStr, name, extension] = fileparts(which(mfilename));
% $$$     pathStr = strrep(pathStr, '\', '/');
% $$$     cd(pathStr)
% $$$     cd ..
% $$$     
% $$$     try
% $$$         % Verifies that repository is clean
% $$$         [status, cmdout] = system('git diff --exit-code');
% $$$         assert(status == 0, 'The repository has local unstaged changes.');
% $$$         [status, cmdout] = system('git diff --cached --exit-code');
% $$$         assert(status == 0, 'The repository has staged but uncommitted changes.');
% $$$         [status, masterRef] = system('git rev-parse --verify master');
% $$$         [status, originMasterRef] = system('git rev-parse --verify origin/master');
% $$$         assert(isequal(masterRef, originMasterRef), 'The repository has committed changes that have not been pushed.');
% $$$         [major minor patch snapshot txt] = replab_version;
% $$$         snapshot = false;
% $$$         txt = textVersion(major, minor, patch, snapshot);
% $$$         userTxt = input(sprintf('Release version [%s]:\n', txt), 's');
% $$$         if ~isempty(userTxt)
% $$$             [major minor patch snapshot txt] = replab_version(userTxt);
% $$$         end
% $$$         major1 = major;
% $$$         minor1 = minor;
% $$$         patch1 = patch + 1;
% $$$         snapshot1 = true;
% $$$         txt1 = textVersion(major1, minor1, patch1, snapshot1);
% $$$         userTxt1 = input(sprintf('Next version [%s]:\n', txt1), 's');
% $$$         if ~isempty(userTxt1)
% $$$             [major1 minor1 patch1 snapshot1 txt1] = replab_version(userTxt1);
% $$$         end
% $$$         disp('Preparing release');
% $$$         disp(sprintf('Setting version to %s', txt));
% $$$         writeVersion(txt);
% $$$         disp('Committing version change.');
% $$$         status = system('git add replab_version.txt');
% $$$         assert(status == 0);
% $$$         status = system('git add docs/_config.yml');
% $$$         assert(status == 0);
% $$$         status = system(sprintf('git commit -m "Setting version to %s"', txt));
% $$$         assert(status == 0);
% $$$         status = system('git push origin', '-echo');
% $$$         assert(status == 0);
% $$$         status = system(sprintf('git tag v%s', txt));
% $$$         assert(status == 0);
% $$$         status = system(sprintf('git push origin v%s', txt), '-echo');
% $$$         assert(status == 0);
% $$$         disp(sprintf('Advancing to next version number %s', txt1));
% $$$         writeVersion(txt1);
% $$$         disp('Committing version change.');
% $$$         status = system('git add replab_version.txt');
% $$$         assert(status == 0);
% $$$         status = system('git add docs/_config.yml');
% $$$         assert(status == 0);
% $$$         status = system(sprintf('git commit -m "Setting version to %s"', txt1));
% $$$         assert(status == 0);
% $$$         status = system('git push origin', '-echo');
% $$$         assert(status == 0);
% $$$     catch me
% $$$         % go back to the initial path
% $$$         cd(initialPath);
% $$$ 
% $$$         error(me.message);
% $$$     end
% $$$     
% $$$     % go back to the initial path
% $$$     cd(initialPath);
% $$$     
% $$$     
% $$$     return;
% $$$     
% $$$ 
% $$$     
    function writeVersion(txtV)
        fid = fopen('replab_version.txt', 'w');
        assert(fid ~= -1);
        fprintf(fid, '%s\n', txt);
        status = fclose(fid);
        assert(status == 0);
        contents = fileread('docs/_config.yml');
        lines = strsplit(contents, '\n');
        tag = 'replabVersion:';
        for i = 1:length(lines)
            L = lines{i};
            if length(L) > length(tag) && isequal(L(1:length(tag)), tag)
                lines{i} = sprintf('replabVersion: %s', txtV);
            end
        end
        contents = strjoin(lines, '\n');
        fid = fopen('docs/_config.yml', 'w');
        assert(fid ~= -1);
        fprintf(fid, '%s', contents);
        status = fclose(fid);
        assert(status == 0);        
    end
    function txtV = textVersion(majorV, minorV, patchV, snapshotV)
        if snapshotV
            suffixV = '-SNAPSHOT';
        else
            suffixV = '';
        end
        txtV = sprintf('%d.%d.%d%s', majorV, minorV, patchV, suffixV);
    end

end

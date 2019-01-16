function replab_release
% Automates the release of a RepLAB version
    [pathStr, name, extension] = fileparts(which(mfilename)); % Verifies current directory
    assert(isequal(pathStr, pwd), 'replab_release must be run from the RepLAB folder');
    % Verifies that repository is clean
    [status, cmdout] = system('git diff --exit-code');
    assert(status == 0, 'The repository has local unstaged changes.');
    [status, cmdout] = system('git diff --cached --exit-code');
    assert(status == 0, 'The repository has staged but uncommitted changes.');
    [status, masterRef] = system('git rev-parse --verify master');
    [status, originMasterRef] = system('git rev-parse --verify origin/master');
    assert(isequal(masterRef, originMasterRef), 'The repository has committed changes that have not been pushed.');
    [major minor patch snapshot txt] = replab_version;
    snapshot = false;
    txt = textVersion(major, minor, patch, snapshot);
    userTxt = input(sprintf('Release version [%s]:\n', txt), 's');
    if ~isempty(userTxt)
        [major minor patch snapshot txt] = replab_version(userTxt);
    end
    major1 = major;
    minor1 = minor;
    patch1 = patch + 1;
    snapshot1 = true;
    txt1 = textVersion(major1, minor1, patch1, snapshot1);
    userTxt1 = input(sprintf('Next version [%s]:\n', txt1), 's');
    if ~isempty(userTxt1)
        [major1 minor1 patch1 snapshot1 txt1] = replab_version(userTxt1);
    end
    disp('Preparing release');
    disp(sprintf('Setting version to %s', txt));
    writeVersion(txt);
    disp('Committing version change.');
    status = system('git add replab_version.txt');
    assert(status == 0);
    status = system(sprintf('git commit -m "Setting version to %s"', txt));
    assert(status == 0);
    status = system('git push origin', '-echo');
    assert(status == 0);
    status = system(sprintf('git tag v%s', txt));
    assert(status == 0);
    status = system(sprintf('git push origin v%s', txt), '-echo');
    assert(status == 0);
    disp(sprintf('Advancing to next version number %s', txt1));
    writeVersion(txt1);
    disp('Committing version change.');
    status = system('git add replab_version.txt');
    assert(status == 0);
    status = system(sprintf('git commit -m "Setting version to %s"', txt1));
    assert(status == 0);
    status = system('git push origin', '-echo');
    assert(status == 0);
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

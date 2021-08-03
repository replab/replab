function replab_release
% Releases the current develop branch contents as a stable release
%
% This runs the release process to send the current ``develop``
% branch snapshot to the branch ``master``, taking care of version numbers in the process.
%
% The documentation is no longer generated as we migrated this part to GitHub actions.
%
% This script works fully offline, interactions with the remote repository are done manually
% by the user.
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
% We assume the following:
%
% * The remote ``origin`` corresponds to the repository `<git@github.com:replab/replab.git>`_
% * The user has run ``git fetch origin master develop`` to synchronize the
%   remote content.
% * The ``develop`` branch is checked out in the main worktree, that is in the root RepLAB folder
% * All of the branches ``master``, ``develop`` have their local copy current with the ``origin`` remote.
%
% All the points above are checked by the release process, apart from whether the ``origin`` remote
% has been fetched, as this would require a connection to the remote.
%
% The release steps are as follows.
%
% 0. We ask the user to confirm they ran ``git fetch origin develop master gh-pages``.
%
% 1. We verify that all working trees do not have uncommited changes.
%
% 2. We verify that the branches ``develop`` and ``master`` are in sync with
%    the ``origin`` remote. We verify that the submodule commits are in sync with the data
%    in the ``external/modules.ini`` file.
%
% 3. We ask the user to confirm the version number of the stable release (by default,
%    the develop ``-SNAP`` version with the ``-SNAP`` suffix removed), and the number
%    of the next develop version (by default, the current version number with the
%    minor release number incremented, and the ``-SNAP`` suffix added).
%
%    We set the version number to the stable number in the develop branch.
%
% 4. We run the `replab_runtests` script. We abort in case of errors.
%
% 5. We run the `replab_checkhelp` script.
%
% 6. We set the version number of the develop branch to a stable number. We commit the release
%    on the develop branch.
%
% 7. We checkout the master branch, merge the develop branch by fast-forward.
%
% 8. We tag the master HEAD with a version tag of the form ``vMAJOR.MINOR.PATH``, as in
%    ``v0.5.0``, which is the format that GitHub recognizes for releases.
%
% 9. We checkout the develop branch, set the version number to the next snapshot number,
%    and commit.
%
% 10. We output the command that the user should run to push the changes to the remote.
%     In particular, it involves pushing both the master and develop branches, and the
%     release tag.

    input('Step 0: Press ENTER to confirm that you ran "git fetch origin develop master"');

    path = replab.globals.replabPath;
    gitDir = fullfile(path, '.git');
    mainWT = path;

    wts = replab.infra.Git.listWorktrees(gitDir);
    mainInd = find(cellfun(@(x) isequal(x, mainWT), {wts.worktree}));
    if length(mainInd) ~= 1
        error('Cannot find main working tree %s in worktree list', mainWT);
    end
    assert(isequal(wts(mainInd).branch, 'refs/heads/develop'), 'Main worktree should have the develop branch checked out');
    mainHead = wts(mainInd).head;

    disp(' ');
    disp('Step 1: Verifying that the current working tree and index are clean');
    assert(~replab.infra.Git.hasUnstagedChanges(gitDir, mainWT), ...
           'The main worktree has local unstaged changes. Verify you followed the installation instructions on the website.');
    assert(~replab.infra.Git.hasStagedUncommittedChanges(gitDir, mainWT), ...
           'The main worktree has staged but uncommitted changes. Verify you followed the installation instructions on the website.');

    disp(' ');
    disp('Step 2: Verifying that master and develop branches are in sync with remote origin.');
    assert(isequal(replab.infra.Git.showExactRef(gitDir, 'refs/heads/master'), ...
                   replab.infra.Git.showExactRef(gitDir, 'refs/remotes/origin/master')), ...
           'Please synchronize master with origin/master');
    assert(isequal(replab.infra.Git.showExactRef(gitDir, 'refs/heads/develop'), ...
                   replab.infra.Git.showExactRef(gitDir, 'refs/remotes/origin/develop')), ...
           'Please synchronize develop with origin/develop');

    % Verify submodules' commits
    submodules = replab.infra.Git.executeAssertExitCodeZero('git submodule status', gitDir);
    submodules = cellfun(@strtrim, strsplit(submodules, {'\r' '\n'}), 'uniform', 0);
    mask = cellfun(@isempty, submodules);
    submodules = submodules(~mask);
    for i = 1:length(submodules)
        parts = strsplit(submodules{i}, ' ');
        commit = parts{1};
        module = parts{2};
        assert(isequal(module(1:9), 'external/'), 'There exists a Git submodule not in external/');
        module = module(10:end);
        data = replab.init.ExternalDependency.loadConfig(module);
        assert(isequal(data.commit, commit), 'Git commit differs for %s', module);
    end

    disp(' ');
    disp('Step 3: New version numbers');
    currentVersion = replab_Version.current;
    assert(currentVersion.snapshot, 'Current develop version must be a snapshot');
    releaseVersion = currentVersion.asRelease.prompt('Release version');
    assert(~releaseVersion.snapshot, 'Updated release version cannot be a snapshot');
    newDevelopVersion = releaseVersion.incrementedPatch.asSnapshot.prompt('Develop version').asSnapshot;

    disp(' ');
    disp('Step 4: Run "replab_runtests"');
    assert(replab_runtests, 'Tests failed');

    disp(' ');
    disp('Step 5: Run "replab_checkhelp"');
    assert(replab_checkhelp, 'Help check failed');

    disp(' ');
    disp('Step 6: Set version number to stable release number, commit to the develop branch');
    releaseVersion.updateVersionFile;
    replab.infra.Git.addAll(gitDir, mainWT);
    replab.infra.Git.commit(gitDir, mainWT, sprintf('Version %s', releaseVersion.toText));

    disp(' ');
    disp('Step 7: Merge the stable release from the develop branch unto the master branch');
    replab.infra.Git.checkout(gitDir, mainWT, 'master')
    replab.infra.Git.merge(gitDir, mainWT, 'develop')

    disp(' ');
    disp('Step 8: Tag the stable release');
    replab.infra.Git.tag(gitDir, mainWT, releaseVersion.tag);
    status = system(sprintf('git tag %s', releaseVersion.tag));

    disp(' ');
    disp('Step 9: Checkout the develop branch, set the version number to the next snapshot');
    replab.infra.Git.checkout(gitDir, mainWT, 'develop');
    newDevelopVersion.updateVersionFile;
    replab.infra.Git.addAll(gitDir, mainWT);
    replab.infra.Git.commit(gitDir, mainWT, sprintf('Version %s', newDevelopVersion.toText));

    disp(' ');
    disp('Step 10: Code to copy/paste');
    disp(' ');
    fprintf('git push origin develop master %s\n', releaseVersion.tag);
end

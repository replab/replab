classdef Git < handle
% Interface to call Git from MATLAB
%
% Method names are based on the Git command names, sometimes with name parts shuffled around for readability
% i.e. ``worktree list`` becomes ``listWorktrees``.
%
% Every static method takes a ``gitDir`` argument so that the current directory does not need to match the Git
% repository being queried.

    methods (Static)

        function [cmdout, exitCode] = execute(cmd, varargin)
            sub = cmd;
            for i = 1:length(varargin)
                arg = ['"' varargin{i} '"'];
                sub = strrep(sub, sprintf('$%d', i), arg);
            end
            [exitCode, cmdout] = system(sub);
        end

        function cmdout = executeAssertExitCodeZero(cmd, varargin)
            sub = cmd;
            for i = 1:length(varargin)
                arg = ['"' varargin{i} '"'];
                sub = strrep(sub, sprintf('$%d', i), arg);
            end
            [exitCode, cmdout] = system(sub);
            if length(varargin) > 0
                assert(exitCode == 0, 'Error while executing "%s" with arguments %s', cmd, strjoin(varargin, ', '));
            else
                assert(exitCode == 0, 'Error while executing "%s"', cmd);
            end
        end

        function l = hasUnstagedChanges(gitDir, worktree)
        % Checks whether the worktree has unstaged changes
        %
        % Args:
        %   gitDir (charstring): Full path to the ``.git`` directory
        %   worktree (charstring): Full path to the worktree directory
        %
        % Returns:
        %   logical: Whether the worktree has unstaged changes
            [cmdout, exitCode] = replab.infra.Git.execute('git --git-dir=$1 --work-tree=$2 diff --exit-code', gitDir, worktree);
            l = (exitCode ~= 0);
        end

        function l = hasStagedUncommittedChanges(gitDir, worktree)
        % Checks whether the worktree has staged, uncommitted changes
        %
        % Args:
        %   gitDir (charstring): Full path to the ``.git`` directory
        %   worktree (charstring): Full path to the worktree directory
        %
        % Returns:
        %   logical: Whether the worktree has unstaged changes
            [cmdout, exitCode] = replab.infra.Git.execute('git --git-dir=$1 --work-tree=$2 diff --cached --exit-code', gitDir, worktree);
            l = (exitCode ~= 0);
        end

        function addAll(gitDir, worktree)
        % Runs ``git add -A``
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 --work-tree=$2 add -A');
        end

        function commit(gitDir, worktree, message)
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 --work-tree=$2 commit -m $3', ...
                                                              gitDir, worktree, message);
        end

        function checkout(gitDir, worktree, ref)
        % Checks out the given reference in the given worktree
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 --work-tree=$2 checkout $3', ...
                                                              gitDir, worktree, ref);
        end

        function merge(gitDir, worktree, ref)
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 --work-tree=$2 merge $3', ...
                                                              gitDir, worktree, ref);
        end

        function tag(gitDir, worktree, tag)
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 --work-tree=$2 tag $3', ...
                                                              gitDir, worktree, tag);
        end

        function ref = showExactRef(gitDir, ref)
        % Returns the commit ID associated with the given exact reference
            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 show-ref --verify $2', gitDir, ref);
            parts = strsplit(cmdout, ' ');
            assert(length(parts) == 2);
            ref = parts{1};
        end

        function wt = listWorktrees(gitDir)
        % List all worktrees of a Git repository
        %
        % Args:
        %   gitDir (charstring): Path to the ``.git`` directory
        %
        % Returns:
        %   struct: Struct with fields ``worktree`` (charstring), ``head`` (charstring), ``branch`` (charstring)


            cmdout = replab.infra.Git.executeAssertExitCodeZero('git --git-dir=$1 worktree list --porcelain', gitDir);
            cmdout = strrep(cmdout, '\r\n', '\n');
            lines = strsplit(cmdout, '\n', 'CollapseDelimiters', false);
            i = 1;
            worktrees = cell(1, 0);
            heads = cell(1, 0);
            branches = cell(1, 0);
            while i <= length(lines)
                l = lines{i};
                if isempty(l)
                    break
                end
                assert(replab.compat.startsWith(l, 'worktree '), 'Cannot parse worktree output %s', cmdout);
                worktrees{1, end+1} = l(10:end);
                i = i + 1;

                l = lines{i};
                assert(replab.compat.startsWith(l, 'HEAD '), 'Cannot parse worktree output %s', cmdout);
                heads{1, end+1} = l(6:end);
                i = i + 1;

                l = lines{i};
                if replab.compat.startsWith(l, 'detached')
                    branches{1, end+1} = '';
                elseif replab.compat.startsWith(l, 'branch ')
                    branches{1, end+1} = l(8:end);
                else
                    error('Cannot parse worktree output %s', cmdout);
                end
                i = i + 1;

                l = lines{i};
                assert(isempty(l), 'Cannot parse worktree output %s', cmdout);
                i = i + 1;
            end
            wt = struct('worktree', worktrees, 'head', heads, 'branch', branches);
        end

    end

end

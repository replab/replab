# How to release RepLAB

1) Change the current directory to the replab base directory.

2) Commit all changes, push to the `origin` remote (which is the remote the release process uses)

3) Run tests under MATLAB (if possible) and Octave locally: `replab_runtests`

4) Run `replab_release`

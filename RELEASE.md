# How to release RepLAB

1) Run tests under MATLAB (if possible) and Octave locally: `replab_runtests`

2) Update `replab_version.txt` and `docs/_config.yml` with the new version number.

3) Commit the changes, and tag the commit (`git tag -a v0.X.0`)

4) Push the changes and the tag to GitHub

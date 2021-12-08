% ---
% jupyter:
%   jupytext:
%     text_representation:
%       extension: .m
%       format_name: light
%       format_version: '1.5'
%       jupytext_version: 1.11.2
%   kernelspec:
%     display_name: Octave
%     language: octave
%     name: octave
% ---

% .. module:: +replab
%
% # How to load/save information to MAT files
%
% As Octave does not support loading/saving of objects, RepLAB implements some methods to store and retrieve information to/from ``.mat`` files.
%
% For now, the following objects implement a pair of ``import``/``export`` methods:
%
% - `+replab.Irreducible` implements `+replab.Irreducible.import` and `+replab.Irreducible.export`. Their usage is documented below.
%
% Beyond those, some of our users report that loading/saving RepLAB objects to ``.mat`` files work in recent versions of MATLAB. Please be careful as the use of ``load`` and ``save`` is currently unsupported. However, if you wish to export computed information about RepLAB objects and those objects are not in the list above, please file an [issue](https://github.com/replab/replab/issues), and we may consider extending our ``import`` / ``export`` mechanism.
%
% ## Importing/exporting information about decompositions
%
% One of the most expensive operations in RepLAB is the decomposition of a representation into irreducibles. Thus, you may wish to save this information in a ``.mat`` file for later use.
%
% RepLAB provides an export mechanism that exposes the information of an irreducible decomposition in plain struct/primitive arrays, so that this information can be saved even when using Octave.
%
% However, the representation being decomposed is not saved by this mechanism, and needs to be recreated using a script when importing the data, as demonstrated below.
%
% ### Save script

addpath([pwd, '/../..']);
replab_init('verbose', 0);
G = replab.S(3);
rep = G.naturalRep; % representation to decompose
dec = rep.decomposition('exact');
decToExport = dec.squeeze % remove empty isotypic components (important when using exact decompositions)
data = decToExport.export;
save decomposition_export_example.mat data

% ### Load script

clear all % pretend we start all over again
replab_init('verbose', 0);
% we construct again the group and its representation
G = replab.S(3);
rep = G.naturalRep; % representation to decompose
data = load('decomposition_export_example.mat');
importedDec = replab.Irreducible.import(rep, data.data)


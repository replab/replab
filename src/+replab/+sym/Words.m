classdef Words
% Sets of 'words' useful for constructing irreducible representations of symmetric groups
%
% The following properties two are described in
% J.D. Wiltshire-Gordon, A. Woo, M. Zajaczkowska, "Specht Polytopes and Specht Matroids" (2017)
% `<https://arxiv.org/abs/1701.05277>`_
%
% There is an easy to construct choice of w1 and w2, as described in Section 1.
% For the examples, we consider partition = [4 2].
%
% Example:
%   >>> partition = [4 2];
%   >>> W = replab.sym.Words(partition);
%   >>> isequal(W.word, [1 1 1 1 2 2])
%       1
%   >>> isequal(W.conjWord, [1 2 3 4 1 2])
%       1
%   >>> isequal(W.dimWord, [4 3 2 1 2 1])
%       1

    properties (SetAccess = protected)
        word % (integer(1,\*)): Word ``w1``
        conjWord % (integer(1,\*)): Word ``w2``
        dimWord % (integer(1,\*)): This is unrelated to the paper, but is used to find the dimension and is constructed similiarly.
    end

    methods

        function self = Words(part,conjPart)
        % Construct the words for a given partition
        %
        % Args:
        %   part (integer(1,\*)): Partition
        %   conjPart (integer(1,\*)): Conjugate partition of part
            if nargin == 1
                conjPart = replab.sym.findPartitions.conjugatePart(part);
            end
            N = numel(part);
            l = sum(part);
            self.word = repelem(1:N,part);
            self.conjWord = zeros(1,l);
            self.dimWord = zeros(1,l);
            index = 1;
            for m = 1:max(conjPart)
                p = sum(conjPart > m - 1);
                self.conjWord(index:index+p-1) = 1:p;
                self.dimWord(index:index+p-1) = p:-1:1;
                index = index + p;
            end
        end

    end

end
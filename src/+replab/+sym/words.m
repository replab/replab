classdef  words
%Sets of 'words' useful for constructing irreducivle
%represnetations of symmetric groups.
    properties (SetAccess = protected)
        %The following properties two are described in 
        %Wiltshire-Gordon, John D.; Woo, Alexander; Zajaczkowska, Magdalena (2017).
        %"Specht Polytopes and Specht Matroids"
        %https://arxiv.org/abs/1701.05277
        %They are an easy to costruct choice of w1 and w2, 
        %as described in Section 1.
        % For the examples, we consider partition = [5 4 2 1]
        word %integer (1,:): Ex: [1 1 1 1 2 2]
        conjWord %integer (1,:): Ex: [1 2 3 4 1 2]
        %This is unrelated to the paper, but is used to find the dimension
        %and is constructed similiarly.
        dimWord %integer (1,:): Ex: [4 3 2 1 2 1]
    end
    
    methods
      function self = words(part,conjPart,type);
      % Construct the words for a given partition
      %
      % Args:
      %   part (integer(1,:)): Partition
      %   conjPart (integer(1,:)): Conjugate partition of part
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
        if nargin == 3
            if type == 'char'
                self.conjWord = char(64+self.conjWord);
                self.word = char(64+self.word);
            end
        end
        end
    end
  
end
classdef findPartitions
    properties
        partCell
        conjPartCell
        partitionHash
        partSet
        nParts
        cycleSizes
        nCycles
    end
    
    methods(Static)
        function conj = conjugatePart(part)
            m = max(part);
            conj = zeros(1,m);
            for j = 1:m
                conj(j) = nnz(j-1<part);
            end
        end
        
        function dim = dimension(part)
            n = sum(part);
            words = replab.sym.words(part,replab.sym.findPartitions.conjugatePart(part));
            columns = zeros(1,n);
            for k = 1:n
                columns(k) = sum(words.conjWord(k+1:n)==words.conjWord(k));
            end    
            dim = round(factorial(n)/prod(words.dimWord+columns));
        end     
        
        function eigVal = eigenvalue(part)
            eigVal = -(sum(part)+(part-2*(1:numel(part)))*part')/2;
        end
        
        function powers = toPowerForm(part)
            n = sum(part);
            powers = zeros(1,n);
            for j = 1:n
                powers(j) = nnz(part == j);
            end
        end
    end
    
    methods
        function self = findPartitions(N)
            powers = generate(N);
            self.partitionHash = replab.perm.Set(N);
            self.partitionHash.insert(powers');  
            powers = fliplr(powers);
            self.nParts = size(powers,1);
            self.cycleSizes = cell(1,self.nParts);
            self.nCycles = cell(1,self.nParts);
            self.nParts = size(powers,1);
            self.partCell = cell(1,self.nParts);
            inds = N:-1:1;
            for k = self.nParts:-1:1
                pow = powers(k,:);
                part =repelem(inds,pow);
                self.nCycles{k} = nonzeros(pow);
                self.cycleSizes{k} = inds(self.nCycles{k});
                self.partCell{k} = part;
                self.conjPartCell{k} = replab.sym.findPartitions.conjugatePart(part);
            end
            function plist = generate(total_sum,candidate_set,max_count,fixed_count)
                % extracts the list of all partitions of a number as integer sums of a list of candidates
                % usage: plist = partitions(total_sum,candidate_set)
                % usage: plist = partitions(total_sum,candidate_set,max_count,fixed_count)
                %
                % PARTITIONS solves the money changing problem. E.g.,
                % how can you make change for one dollar given coins
                % of a given set of denominations. A good reference on
                % the general problem is found here:
                %
                % http://en.wikipedia.org/wiki/Integer_partition
                %
                % PARTITIONS uses a recursive strategy to enumerate all
                % possible partitions of the total_sum. This may be
                % highly intensive for large sums or large sets of
                % candidates.
                %
                % arguments: (input)
                %  total_sum - scalar positive integer (to be partitioned)
                %
                %              BEWARE! a large total_sum can easily cause
                %              stack problems. For example, the number of
                %              partitions of 40 is 37338, a set that took 24
                %              seconds to completely enumerate on my cpu.
                %
                %  candidate_set - (OPTIONAL) vector of (distinct) candidate
                %              positive integers for the partitions.
                %
                %              Efficiency considerations force me to require
                %              that the candidates be sorted in non-decreasing
                %              order. An error is produced otherwise.
                %
                %              DEFAULT: candidate_set = 1:total_sum
                %
                %              BEWARE! large candidate sets can easily cause
                %              stack problems
                %
                %  max_count - (OPTIONAL) the maximum quantity of any
                %              candidate in the final sum.
                %
                %              max_count must be either a vector of the
                %              same length as candidate_set, or a scalar
                %              that applies to all elements in that set.
                %
                %              DEFAULT = floor(total_sum./candidate_set)
                %
                %  fixed_count - (OPTIONAL) Allows you to specify a fixed
                %              number of terms in the partitioned sum.
                %
                %              fixed_count must be a positive integer if
                %              supplied.
                %
                %              DEFAULT = []
                %
                % arguments: (output)
                %  plist - array of partitions of total_sum. This is a list
                %              of the quantity of each element such that
                %              plist*candidate_set(:) yields total_sum
                %
                %
                % Author: John D'Errico
                % e-mail: woodchips@rochester.rr.com
                % Release: 2
                % Release date: 7/15/08
                % default for candidate_set
                if (nargin<2) || isempty(candidate_set)
                candidate_set = 1:total_sum;
                end
                % how many candidates are there
                n = length(candidate_set);
                % error checks
                if any(candidate_set<0)
                error('All members of candidate_set must be >= 0')
                end
                % candidates must be sorted in increasng order
                if any(diff(candidate_set)<0)
                error('Efficiency requires that candidate_set be sorted')
                end
                % check for a max_count. do we supply a default?
                if (nargin<3) || isempty(max_count)
                % how high do we need look?
                max_count = floor(total_sum./candidate_set);
                elseif length(max_count)==1
                % if a scalar was provided, then turn it into a vector
                max_count = repmat(max_count,1,n);
                end
                % check for a fixed_count
                if (nargin<4) || isempty(fixed_count)
                    fixed_count = [];
                elseif (fixed_count<0) || (fixed_count~=round(fixed_count))
                    error('fixed_count must be a positive integer if supplied')
                end
                % check for degenerate cases
                if isempty(fixed_count)
                    if total_sum == 0
                        plist = zeros(1,n);
                    return
                    elseif (n == 0)
                        plist = [];
                    return
                    elseif (n == 1)
                    % only one element in the set. can we form
                    % total_sum from it as an integer multiple?
                        p = total_sum/candidate_set;
                        if (p==fix(p)) && (p<=max_count)
                            plist = p;
                        else
                            plist = [];
                        end
                        return
                    end
                    else
                    % there was a fixed_count supplied
                    if (total_sum == 0) && (fixed_count == 0)
                        plist = zeros(1,n);
                        return
                    elseif (n == 0) || (fixed_count <= 0)
                        plist = [];
                        return
                    elseif (n==1)
                    % there must be a non-zero fixed_count, since
                    % we did not trip the last test. since there
                    % is only one candidate in the set, will it work?
                    if ((fixed_count*candidate_set) == total_sum) && (fixed_count <= max_count)
                        plist = fixed_count;
                    else
                        plist = [];
                    end
                    return
                    end
                end
                % finally, we can do some work. start with the
                % largest element and work backwards
                m = max_count(end);
                % do we need to back off on m?
                c = candidate_set(end);
                m = min([m,floor(total_sum/c),fixed_count]);
                plist = zeros(0,n);
                for i = 0:m
                    temp = generate(total_sum - i*c, ...
                    candidate_set(1:(end-1)), ...
                    max_count(1:(end-1)),fixed_count-i);
                    plist = vertcat(plist,[temp,repmat(i,size(temp,1),1)]);  %#ok
                end
                end
            end
    end
end
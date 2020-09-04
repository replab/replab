function basis = findIsotypicBasis(rep, irrep)
% Returns a basis of an isotypic component, with the multiplicity space identified
% For more background see Serre's textbook Chapter 2 or Faessler and
% Stiefel Chapter 5 on how to generate symmetry adapted basis, in which a
% matrix with the symmetry of the representation becomes block-diagonal.
% 
% Args:
%   irrep (`+replab.Rep`): Irreducible representation
%   rep (`+replab.Rep`): Representation to decompose
%
% Returns:
%   double(\*,\*), may be sparse: Basis of the isotypic component corresponding to ``irrep`` of dimension ``rep.dimension x (multiplicity*irrep.dimension)``
%
% Construct a projector on an irrep and find c_{i} linearly independent
% column vectors of this matrix of dimension ``rep.dimension x
% rep.dimension``: 
Projector = secondProjection(rep,irrep,1,1);
[~,pivot] = rref(Projector{1,1});
vectors = [Projector{1,1}(:,pivot)];
% The choice of the linearly independent columns is arbitrary in general.
% Comment: would it be useful to make it choose say the first 2 linearly
% independent columns that it finds? 
%
% Need to act with the off-diagonal maps to find the other basis vectors
% inside the isotypic component:
id = irrep.dimension; 
v_new = zeros(rep.dimension,1);           
    for j=1:id    
        for k = 1:length(pivot) 
        map = secondProjection(rep,irrep,1,j);
        vec_new = map{1,j}*vectors(:,k);
        v_new = [v_new vec_new];
        end
    end
% This labelling of the basis vectors inside the components is sufficient
% to make the matrix that has the symmetry of the representation block
% diagonal. If you need the images of the representation to be block
% diagonal instead, just interchange the ``for`` loops for j and k. 
% Comment: 
% - we could give the user the chance of selecting which one is of
%   interest for a particular case. 
% - Is there another way to get the images block diagonal without doing
%   this relabeling? 
[~,columns] = size(v_new);
basis = v_new(:,2:columns);
% Returns a basis with the multiplicity space identified in the isotypic
% component, which is isomorphic to a direct sum of c_i copies of
% ``irrep``. 
end
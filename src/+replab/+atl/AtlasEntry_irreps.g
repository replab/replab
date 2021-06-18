LoadPackage( "repsn" );;
LoadPackage( "json" );;

Assert(0, IsBound(G));;

# This script assumes that:
# - The group is given as the variable "G", without making further assumptions on it.
#   It can be a permutation group, a pc-group, etc.
# - The group name is given as the variable "name"

# We use the following letters:
# - G: the original group
# - P: the permutation group isomorphic to G
# - F: the finitely presented group isomorphic to G
# - mu_*: Isomorphisms

# We identify the group

id := IdGroup(G);;

# First, we reduce the number of generators

Ggens := SmallGeneratingSet(G);;

# Then we find a permutation group with a small degree

mu_G_P1 := IsomorphismPermGroup(G);;
mu_P1_P := SmallerDegreePermutationRepresentation(Image(mu_G_P1));;
mu_G_P := CompositionMapping(mu_P1_P, mu_G_P1);;
P := Range(mu_G_P);;

# Then we find the finitely presented group
mu_G_F1 := IsomorphismFpGroupByGenerators(G, Ggens);;
mu_F1_F := IsomorphismSimplifiedFpGroup(Image(mu_G_F1));;
mu_G_F := CompositionMapping(mu_F1_F, mu_G_F1);;
F := Range(mu_G_F);;

# Compute relators of F
F_rels := List(RelatorsOfFpGroup(F), LetterRepAssocWord);;

# Compute generators of P, which match the generators of P
F_gens := GeneratorsOfGroup(F);;
nGens := Size(F_gens);;
G_gens := List(F_gens, f -> PreImageElm(mu_G_F, f));;
P_gens := List(G_gens, f -> ImageElm(mu_G_P, f));;
P_degree := LargestMovedPoint(P_gens);;
permimg := g -> ListPerm(Inverse(g), P_degree);;
P_gen_images := List(P_gens, permimg);;

# Group properties
derivedSeriesAbelianInvariants := List(DerivedSeriesOfGroup(G), AbelianInvariants);;

# Character table
G_ct := CharacterTable(G);;
G_chars := Irr(G_ct);;
G_cc_images := List(ConjugacyClasses(G_ct), c -> permimg(ImageElm(mu_G_P, Representative(c))));;
G_generatorNames := List([ 1 .. nGens ], i -> Concatenation("x", String(i)));;
G_char_values := List(G_chars, c -> List(ValuesOfClassFunction(c), String));;

# G_auto := AutomorphismGroup(prm_grp);;
# innerAuto := InnerAutomorphismsAutomorphismGroup(auto);;
# outerAuto := NaturalHomomorphismByNormalSubgroup(auto, innerAuto);;
# outerImg := IsomorphismPermGroup(Image(outerAuto));;
# outerImg1 := SmallerDegreePermutationRepresentation(Range(outerImg));;

permutationGenerators := P_gen_images;;
classes := G_cc_images;;
characters := G_char_values;;
order := String(Size(G));;
irreps := List(G_chars, c -> IrreducibleAffordingRepresentation(c));;
irrepImages := List(irreps, rep -> List(G_gens, g -> ImageElm(rep, g)));;
irrepStrings := List(irrepImages, lvl1 -> List(lvl1, lvl2 -> List(lvl2, lvl3 -> List(lvl3, String))));;

name := StringFormatted("SmallGroup({},{})", id[1], id[2]);;

group := rec( name := name, generatorNames := G_generatorNames, permutationGenerators := permutationGenerators, relators := F_rels, order := order, classes := classes, derivedSeriesAbelianInvariants := derivedSeriesAbelianInvariants );;
complexCharacterTable := rec( characters := characters, irreps := irrepStrings );;
data := rec( group := group, complexCharacterTable := complexCharacterTable );;
GapToJsonStream(OutputTextUser(), data);;

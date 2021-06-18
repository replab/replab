LoadPackage( "repsn" );;
LoadPackage( "json" );;

Iprm1 := IsomorphismPermGroup(G);;
Iprm2 := SmallerDegreePermutationRepresentation(Image(Iprm1));;
Iprm := CompositionMapping(Iprm2, Iprm1);;

Ifp1 := IsomorphismFpGroup(G);;
Ifp2 := IsomorphismSimplifiedFpGroup(Image(Ifp1));;

rels := List(RelatorsOfFpGroup(Range(Ifp2)), LetterRepAssocWord);;

derivedSeriesAbelianInvariants := List(DerivedSeriesOfGroup(G), g -> AbelianInvariants(g));;

FP_gens := List(GeneratorsOfGroup(Range(Ifp2)), g -> PreImageElm(Ifp2, g));;
G_gens := List(FP_gens, g -> PreImageElm(Ifp1, g));;
prm_gens :=  List(G_gens, g -> ImageElm(Iprm, g));;
prm_grp := Group(prm_gens);;
n := LargestMovedPoint(prm_gens);;
prm_img_gens := List(prm_gens, g -> ListPerm(Inverse(g), n));;

tbl := CharacterTable(G);;
chars := Irr(tbl);;
C := List(ConjugacyClasses(tbl), c -> ListPerm(Inverse(ImageElm(Iprm, Representative(c))), n));;

values := List(chars, c -> List(ValuesOfClassFunction(c), String));;

auto := AutomorphismGroup(prm_grp);;
innerAuto := InnerAutomorphismsAutomorphismGroup(auto);;
outerAuto := NaturalHomomorphismByNormalSubgroup(auto, innerAuto);;
outerImg := IsomorphismPermGroup(Image(outerAuto));;
outerImg1 := SmallerDegreePermutationRepresentation(Range(outerImg));;

generatorNames := List([ 1 .. Size(prm_img_gens) ], i -> Concatenation("x", String(i)));;
permutationGenerators := prm_img_gens;;
classes := C;;
characters := values;;
order := String(Size(G));;
irreps := List(chars, c -> IrreducibleAffordingRepresentation(c));;
irrepImages := List(irreps, rep -> List(prm_gens, g -> ImageElm(rep, PreImageElm(Iprm, g))));;
irrepStrings := List(irrepImages, lvl1 -> List(lvl1, lvl2 -> List(lvl2, lvl3 -> List(lvl3, String))));;


group := rec( generatorNames := generatorNames, permutationGenerators := permutationGenerators, relators := rels, order := order, classes := classes, derivedSeriesAbelianInvariants := derivedSeriesAbelianInvariants );;
complexCharacterTable := rec( characters := characters, irreps := irrepStrings );;
data := rec( group := group, complexCharacterTable := complexCharacterTable );;
GapToJsonStream(OutputTextUser(), data);;

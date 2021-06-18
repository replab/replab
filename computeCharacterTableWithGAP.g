LoadPackage( "repsn" );;
LoadPackage( "json" );;

tbl := CharacterTable(G);;
chars := Irr(tbl);;

irreps := List(chars, c -> IrreducibleAffordingRepresentation(c));;
irrepImages := List(irreps, rep -> List(GeneratorsOfGroup(G), g -> ImageElm(rep, g)));;
irrepStrings := List(irrepImages, lvl1 -> List(lvl1, lvl2 -> List(lvl2, lvl3 -> List(lvl3, String))));;

GapToJsonStream(OutputTextUser(), irrepStrings);;

<TeXmacs|1.99.18>

<style|generic>

<\body>
  <doc-data|<doc-title|Real representations>|<doc-author|<author-data|<author-name|Denis
  Rosset>>>|<doc-date|January 26, 2021>>

  <section|Definitions>

  For a real or complex matrix <math|X>, we write <math|X<rsup|\<top\>>> its
  transpose, <math|<overline|X>> its conjugate and
  <math|X<rsup|\<dag\>>=<wide|X|\<bar\>><rsup|\<top\>>> its conjugate
  transpose.

  We write <math|\<rho\>:G\<rightarrow\><math-ss|GL><around*|(|V|)>> a
  representation with <math|V=<with|font|Bbb|R><rsup|D>>.\ 

  A linear equivariant map <math|f:V\<rightarrow\>V<rprime|'>>, represented
  by a matrix <math|F> is such that <math|f<around*|(|<wide|v|\<vect\>>|)>=F*<wide|v|\<vect\>>>,
  with the equivariant condition <math|\<rho\><rsub|g><rprime|'>\<circ\>f=f\<circ\>*\<rho\><rsub|g>>
  written in matrix form:

  <\equation>
    \<rho\><rsub|g><rprime|'>*F=F*\<rho\><rsub|g><space|2em><text|for all
    >g\<in\>G.
  </equation>

  A specific equivariant space is the commutant space
  <math|<math-ss|C><rsub|\<rho\>>=<around*|{|X:\<rho\><rsub|g>X=X*\<rho\><rsub|g>,\<forall\>g\<in\>G|}>>.

  An antilinear equivariant map <math|f:V\<rightarrow\>V<rprime|'>> is
  represented by a matrix <math|F> such that
  <math|f<around*|(|<wide|v|\<vect\>>|)>=<overline|F*<wide|v|\<vect\>>>>,
  with the equivariant condition <math|\<rho\><rsub|g><rprime|'>\<circ\>f=f\<circ\>*\<rho\><rsub|g>>
  written

  <\equation>
    \<rho\><rsub|g><rprime|'><overline|F<wide|v|\<vect\>>>=<overline|F\<rho\><rsub|g><wide|v|\<vect\>>>,<space|2em><overline|\<rho\><rprime|'><rsub|g>>F=F*\<rho\><rsub|g>.
  </equation>

  We will also consider representations over
  <math|<with|font|Bbb|C><rsup|d>>. For a complex representation
  <math|\<sigma\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|C><rsup|d>|)>>,
  we write <math|<wide|\<sigma\>|\<bar\>>> its conjugate transpose and
  <math|\<sigma\><rsup|+>> with <math|\<sigma\><rsup|+><rsub|g>=<around*|(|\<sigma\><rsub|g<rsup|-1>>|)><rsup|\<top\>>>
  its dual.

  <section|Interpreting the result of EV decomposition>

  We sample a generic element <math|X\<in\><math-ss|C><rsub|\<rho\>>>, and
  write its eigenvalue decomposition either as:

  <\equation>
    U<rsup|-1>*X*U=<matrix|<tformat|<table|<row|<cell|\<lambda\><with|font|Bbb|1><rsub|d><rsub|>>|<cell|>>|<row|<cell|>|<cell|\<ddots\>>>>>><space|2em><text|or><space|2em>U<rsup|-1>*X*U=<matrix|<tformat|<table|<row|<cell|<around*|(|\<alpha\>+\<beta\>
    b|)><with|font|Bbb|1><rsub|d><rsub|>>|<cell|>|<cell|>>|<row|<cell|>|<cell|<around*|(|\<alpha\>-\<beta\>
    b|)><with|font|Bbb|1><rsub|d>>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ddots\>>>>>>,
  </equation>

  where <math|\<lambda\>>, <math|\<alpha\>> and <math|\<beta\>> are real, and
  the <math|<around*|(|\<ddots\>|)>> block denotes other distinct eigenvalues
  of <math|X>, and the corresponding eigenspace describes a subrepresentation
  of <math|\<rho\>> that can be treated separately by induction. We thus
  focus on the other blocks.

  In the first case, <math|\<lambda\>> is a real eigenvalue (multiple, when
  <math|d\<gtr\>1>), and its eigenspace defines a real-type irreducible
  subrepresentation.

  In the second case, <math|\<alpha\>\<pm\>i \<beta\>> is a pair of conjugate
  complex eigenvalues, corresponding to a pair of complex subrepresentations
  <math|\<sigma\>> and <math|\<tau\>>. Those subrepresentations do not
  necessary split also over the reals. We distinguish three cases.

  In the first case, the subrepresentations <math|\<sigma\>> and
  <math|\<tau\>> are inequivalent. Then <math|\<sigma\>\<oplus\>\<tau\>> is a
  complex-type real representation.

  Otherwise, the subrepresentations <math|\<sigma\>:G\<rightarrow\><math-ss|GL><around*|(|W|)>>
  and <math|\<tau\>> are equivalent. We then look for an equivariant
  antilinear map <math|j:W\<rightarrow\>W>, represented by a matrix <math|J>
  such that <math|<overline|\<sigma\><rsub|g>>J=J*\<sigma\><rsub|g>>. We
  distinguish two cases:

  <\itemize-minus>
    <item>We have <math|j<rsup|2>=\<mu\>1>, i.e.
    <math|j<around*|(|j<around*|(|<wide|v|\<vect\>>|)>|)>=\<mu\><wide|v|\<vect\>>>,
    that is <math|j<around*|(|<overline|J<wide|v|\<vect\>>>|)>=<overline|J\<cdot\><around*|(|<overline|J<wide|v|\<vect\>>>|)>>=<wide|J|\<bar\>>J<wide|v|\<vect\>>=\<mu\><wide|v|\<vect\>>>
    and <math|<wide|J|\<bar\>>J=\<mu\><with|font|Bbb|1>> for
    <math|\<mu\>\<gtr\>0>. Then <math|\<sigma\>\<oplus\>\<tau\>> contains two
    real representations of real-type.

    <item>We have <math|j<rsup|2>=-\<mu\>1> for <math|\<mu\>\<gtr\>0>. Then
    <math|\<sigma\>\<oplus\>\<tau\>> contains a real representation of
    quaternionic-type, and <math|J> describes the quaternionic structure.
  </itemize-minus>

  <section|Complex EV, distinguishing the three cases>

  Without loss of generality, we restrict <math|\<rho\>> such that the
  generic <math|X\<in\><math-ss|C><rsub|\<rho\>>> splits with
  <math|\<lambda\>=\<alpha\>+i*\<beta\>>:

  <\equation>
    U<rsup|-1>*X*U=<matrix|<tformat|<table|<row|<cell|\<lambda\><with|font|Bbb|1><rsub|d><rsub|>>|<cell|>>|<row|<cell|>|<cell|<wide|\<lambda\>|\<bar\>><with|font|Bbb|1><rsub|d>>>>>>,
  </equation>

  and we write, with real matrices <math|A,B,Y,Z>:

  <\equation>
    <matrix|<tformat|<table|<row|<cell|Y+i Z>>|<row|<cell|Y-i
    Z>>>>>X<matrix|<tformat|<table|<row|<cell|<around*|(|A+i
    B|)>>|<cell|<around*|(|A-i B|)>>>>>>=<matrix|<tformat|<table|<row|<cell|\<lambda\><with|font|Bbb|1><rsub|d><rsub|>>|<cell|>>|<row|<cell|>|<cell|<wide|\<lambda\>|\<bar\>><with|font|Bbb|1><rsub|d>>>>>>,
  </equation>

  with

  <\equation>
    <matrix|<tformat|<table|<row|<cell|Y>>|<row|<cell|Z>>>>><matrix|<tformat|<table|<row|<cell|A>|<cell|B>>>>>=<frac|1|2><matrix|<tformat|<table|<row|<cell|<with|font|Bbb|1>>|<cell|>>|<row|<cell|>|<cell|-<with|font|Bbb|1>>>>>>.
  </equation>

  We have two subrepresentations

  <\equation>
    \<sigma\><rsub|g>=<around*|(|Y+i*Z|)>\<rho\><rsub|g><around*|(|A+i
    B|)>=<around*|(|Y\<rho\><rsub|g>A-Z\<rho\><rsub|g>B|)>+i<around*|(|Y\<rho\><rsub|g>B+Z*\<rho\><rsub|g>A|)>,
  </equation>

  and

  <\equation>
    \<tau\><rsub|g>=<around*|(|Y-i*Z|)>\<rho\><rsub|g><around*|(|A-i
    B|)>=<around*|(|Y\<rho\><rsub|g>A-Z\<rho\><rsub|g>B|)>-i<around*|(|Y\<rho\><rsub|g>B+Z*\<rho\><rsub|g>A|)>,
  </equation>

  and we see that <math|\<tau\><rsub|g>=<wide|\<sigma\>|\<bar\>><rsub|g>>.
  Now, given a generic <math|X<rprime|'>>, we define <math|F> as:

  <\equation>
    F=<around*|(|Y+i Z|)>X<rprime|'><around*|(|A-i*B|)>,
  </equation>

  and remark that

  <\eqnarray>
    <tformat|<table|<row|<cell|\<sigma\><rsub|g>F>|<cell|=>|<cell|<around*|(|Y+i*Z|)>\<rho\><rsub|g><around*|(|A+i
    B|)><around*|(|Y+i Z|)>X<rprime|'><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i*Z|)><around*|(|A+i
    B|)><around*|(|Y+i Z|)>\<rho\><rsub|g>X<rprime|'><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i
    Z|)>\<rho\><rsub|g>X<rprime|'><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i
    Z|)>X<rprime|'>\<rho\><rsub|g><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i
    Z|)>X<rprime|'>\<rho\><rsub|g><around*|(|A-i*B|)><around*|(|Y-i*Z|)><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i
    Z|)>X<rprime|'><around*|(|A-i*B|)><around*|(|Y-i*Z|)>\<rho\><rsub|g><around*|(|A-i*B|)>>>|<row|<cell|>|<cell|=>|<cell|<around*|(|Y+i
    Z|)>X<rprime|'><around*|(|A-i*B|)>\<tau\><rsub|g>>>|<row|<cell|>|<cell|=>|<cell|F*\<tau\><rsub|g>,>>>>
  </eqnarray>

  as <math|<around*|(|Y\<pm\>i*Z|)><around*|(|A\<pm\>i
  B|)>=<with|font|Bbb|1>>, while <math|\<Pi\><rsub|\<pm\>>\<equiv\><around*|(|A\<pm\>i
  B|)><around*|(|Y\<pm\>i*Z|)>> is a projector on an invariant subspace of
  <math|\<rho\><rsub|g>> and thus commutes with <math|\<rho\>>. Thus <math|F>
  is an equivariant map from <math|\<tau\>> to <math|\<sigma\>>. As we know
  that <math|\<tau\>=<wide|\<sigma\>|\<bar\>>>, <math|F> also describes an
  antilinear map from <math|\<sigma\>> to <math|\<sigma\>>.

  Thus we have <math|<wide|F|\<bar\>>F=\<mu\><with|font|Bbb|1>>, with
  <math|\<mu\>\<in\><with|font|Bbb|R>>. If <math|\<mu\>=0>, then
  <math|\<rho\>> is a real irrep of complex type. If <math|\<mu\>\<gtr\>0>,
  then <math|\<rho\>> contains two copies of a real-type irrep. If
  <math|\<mu\>\<less\>0>, then <math|\<rho\>> is a real irrep of quaternionic
  type.

  Using real arithmetic only:

  <\equation>
    F=<around*|(|Y+i Z|)>X<rprime|'><around*|(|A-i*B|)>=Y X<rprime|'> A+Z
    X<rprime|'>B+i<around*|(|Z X<rprime|'> A-Y X<rprime|'> B|)>,
  </equation>

  and

  <\equation>
    <wide|F|\<bar\>>F=<around*|(|Y X<rprime|'> A+Z
    X<rprime|'>B|)><rsup|2>+<around*|(|Z X<rprime|'> A-Y X<rprime|'>
    B|)><rsup|2>.
  </equation>
</body>

<\initial>
  <\collection>
    <associate|page-height|auto>
    <associate|page-medium|paper>
    <associate|page-type|letter>
    <associate|page-width|auto>
  </collection>
</initial>

<\references>
  <\collection>
    <associate|auto-1|<tuple|1|1|../../.TeXmacs/texts/scratch/no_name_8.tm>>
    <associate|auto-2|<tuple|2|?|../../.TeXmacs/texts/scratch/no_name_8.tm>>
    <associate|auto-3|<tuple|3|?|../../.TeXmacs/texts/scratch/no_name_8.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Definitions>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
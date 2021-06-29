<TeXmacs|1.99.18>

<style|<tuple|generic|old-lengths>>

<\body>
  <doc-data|<doc-title|Error control in RepLAB>|<doc-author|<author-data|<author-name|Denis
  Rosset>>>|<doc-date|January 11, 2021>>

  <section|Definitions>

  We write <math|G> a compact group with measure <math|\<mu\>>. We write
  <math|<big|int><rsub|G>f<around*|(|g|)> d\<mu\><around*|(|g|)>> the Haar
  integral for a measurable function <math|f>.

  We write <math|<with|font|Bbb*|K>=<with|font|Bbb*|R>,<with|font|Bbb*|C>>
  the field on which representations are defined. Our representations are
  finite-dimensional; thus we can identify the <math|d>-dimensional
  representation vector space <math|V> with
  <math|<with|font|Bbb*|K><rsup|d>>, by specifying a basis on <math|V>. In
  the rest of this document, we will always assume
  <math|V=<with|font|Bbb*|K><rsup|d>> for simplicity so that, for example, a
  map <math|F:V\<rightarrow\>V> can be identified with a matrix
  <math|<around*|(|F<rsub|i j>|)><rsub|i j>> with coefficients in
  <math|<with|font|Bbb*|K>>.

  <subsection|Representations>

  However, beware that the Euclidean inner product structure on <math|V> is
  not always directly meaningful.

  We write <math|<math-ss|Aut><around*|(|V|)>> the linear maps
  <math|V\<rightarrow\>V>, which we identify with the space
  <math|<with|font|Bbb*|K><rsup|d*\<times\>d>> of <math|d\<times\>d> matrices
  with coefficients over <math|<with|font|Bbb*|K>>.

  We write <math|<math-ss|GL><around*|(|V|)>> the group of linear
  automorphisms <math|V\<rightarrow\>V>, which we identify with
  <math|<math-ss|GL><rsub|d><around*|(|<with|font|Bbb*|K>|)>> the group of
  <math|d\<times\>d> invertible matrices with coefficients over
  <math|<with|font|Bbb*|K>>.

  Given a map <math|A:V\<rightarrow\>V> or a matrix <math|<around*|(|A<rsub|i
  j>|)><rsub|i j>>, we write <math|A<rsup|\<dag\>>> its Hermitian adjoint
  given by <math|<around*|(|A<rsup|\<dag\>>|)><rsub|i j>=A<rsub|j
  i><rsup|\<ast\>>>, where <math|<rsup|\<ast\>>> denotes the complex
  conjugate. We also write the transpose <math|A<rsup|\<top\>>> with
  <math|<around*|(|A<rsup|\<top\>>|)><rsub|i j>=A<rsub|j i>>.

  We write <math|<math-ss|U><around*|(|V|)>> or
  <math|<math-ss|U><rsub|d><around*|(|<with|font|Bbb*|K>|)>> the group of
  <math|d\<times\>d> unitary matrices with coefficients over <math|\<bbb-K\>>
  (i.e. those with <math|A<rsup|\<dag\>>=A<rsup|-1>>).

  A representation of <math|G> is a map <math|\<rho\>:G\<rightarrow\><math-ss|GL><around*|(|V|)>>,
  <math|g\<mapsto\>\<rho\><rsub|g>>, such that
  <math|\<rho\><rsub|g<rsub|1>g<rsub|2>>=\<rho\><rsub|g<rsub|1>>\<rho\><rsub|g<rsub|2>>>.
  A representation is unitary if <math|\<rho\><rsub|g<rsup|-1>>=<around*|(|\<rho\><rsub|g>|)><rsup|\<dag\>>>.

  <subsection|Numerics>

  A map <math|<wide|\<rho\>|~>:G\<rightarrow\><math-ss|Aut><around*|(|V|)>>
  is an <em|<math|\<varepsilon\>>-representation> if

  <\equation>
    <around*|\<\|\|\>|<wide|\<rho\>|~><rsub|g<rsub|1>><wide|\<rho\>|~><rsub|g<rsub|2>>-<wide|\<rho\>|~><rsub|g<rsub|1>
    g<rsub|2>>|\<\|\|\>><rsub|F>\<leqslant\>\<varepsilon\>
  </equation>

  for all <math|g<rsub|1>,g<rsub|2>\<in\>G>. We denote by
  <math|<around*|\<\|\|\>|A|\<\|\|\>><rsub|F>=<sqrt|tr*A<rsup|\<dag\>>A>> the
  Frobenius norm.

  The <em|condition number> of a matrix <math|A> is the ratio between its
  maximal and minimal singular value:

  <\equation>
    \<kappa\><around*|(|A|)>=\<sigma\><rsub|max><around*|(|A|)>/\<sigma\><rsub|min><around*|(|A|)>.
  </equation>

  By analogy, we define the <em|condition number> of <math|<wide|\<rho\>|~>>

  <\equation>
    \<kappa\><around*|(|<wide|\<rho\>|~>|)>=max<rsub|g\<in\>G>
    \<kappa\><around*|(|<wide|\<rho\>|~><rsub|g>|)>.
  </equation>

  <subsection|Number data type>

  An <em|exact> number, for our purposes, is a member of one of the
  cyclotomic fields <math|<with|font|Bbb*|Q><around*|(|e<rsup|2\<pi\>i/n>|)>>
  for some <math|n\<geqslant\>1>. This field can encode any irreducible
  representation of a finite group over the real or complex numbers.

  An <em|interval> (or <em|intval>) coefficient is given either: over
  <math|<with|font|Bbb*|R>> by an interval <math|<around*|[|inf,sup|]>>, over
  <math|<with|font|Bbb*|C>> by a midpoint and a radius.

  A floating-point <em|double> is represented using the <verbatim|binary64>
  IEEE 754 encoding with 53 significand bits, corresponding to roughly 16
  decimal digits.

  <section|Constructing a representation>

  <subsection|Defining representation>

  If <math|G> itself is <math|<math-ss|U><rsub|d><around*|(|<with|font|Bbb*|K>|)>>,
  that is the unitary group for <math|<with|font|Bbb*|K>=<with|font|Bbb*|R>>
  and the orthogonal group for <math|<with|font|Bbb*|K>=<with|font|Bbb*|C>>,
  then its defining representation is simplfy <math|\<rho\><rsub|g>=g>.

  <paragraph*|How do we sample from <math|<math-ss|U>>? \V>Sample a matrix
  with Gaussian-distributed real or complex entries, use QR decomposition
  (forcing the R diagonal to be positive).

  <paragraph*|Error model. \V>We use an error radius of <math|10<rsup|-15>>
  on individual coefficients (from our experiment on real and complex
  matrices of sizes <math|2\<times\>2>, <math|4\<times\>4>, ...,
  <math|64\<times\>64>). The error seems to be the largest for small
  matrices. The error, expressed in Frobenius norm, is better than

  <\equation>
    <around*|\<\|\|\>|10<rsup|-15><matrix|<tformat|<table|<row|<cell|1>|<cell|\<ldots\>>|<cell|1>>|<row|<cell|\<ldots\>>|<cell|>|<cell|\<ldots\>>>|<row|<cell|1>|<cell|\<ldots\>>|<cell|1>>>>>|\<\|\|\>><rsub|F>=10<rsup|-15>n,
  </equation>

  but we stay with the simpler error model.

  <paragraph*|Interval arithmetic. \V>Not done.

  <paragraph*|Higher precision. \V>If, one day, we allow the use of GEM in
  RepLAB, we'll consider that the floating-point approximate group element
  <math|U> \Prepresents\Q a \Pclose\Q orthogonal/unitary matrix in
  higher-precision, the higher-precision matrix being obtained
  deterministically (for example from Newton iterations). See
  https://en.wikipedia.org/wiki/Orthogonal_matrix#Nearest_orthogonal_matrix

  <subsection|(Signed) permutation representations>

  If <math|G> is a permutation group or a signed permutation group,\ 

  <subsection|>

  <subsection|From images>

  <section|Remarkable spaces>

  <subsection|Remarkable spaces>

  The commutant <math|<math-ss|C><around*|(|\<rho\>|)>\<subseteq\><math-ss|Aut><around*|(|V|)>>
  is the space of all maps/matrices that commute with <math|\<rho\>>:

  <\equation>
    <math-ss|C><around*|(|\<rho\>|)>=<around*|{|A\<in\><with|font|Bbb*|K><rsup|d\<times\>d>:\<rho\><rsub|g>A
    -A\<rho\><rsub|g>=0,\<forall\>g\<in\>G|}>,
  </equation>

  a condition we also write <math|A=\<rho\><rsub|g>A\<rho\><rsub|g><rsup|-1>.>

  The Hermitian invariant <math|<math-ss|H><around*|(|\<rho\>|)>\<subseteq\><math-ss|Aut><around*|(|V|)>>
  is the space of all maps/matrices that obey a slightly different condition

  <\equation>
    <math-ss|H><around*|(|\<rho\>|)>=<around*|{|A\<in\><with|font|Bbb*|K><rsup|d\<times\>d>:\<rho\><rsub|g>*A-A<around*|(|\<rho\><rsub|g><rsup|-1>|)><rsup|\<dag\>>=0,\<forall\>g\<in\>G|}>
  </equation>

  which is also <math|A=\<rho\><rsub|g>A\<rho\><rsub|g><rsup|\<dag\>>>.

  <appendix|Numerical experiment on unitary sampling>

  \;
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
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|4|2>>
    <associate|auto-11|<tuple|2.2|2>>
    <associate|auto-12|<tuple|2.3|2>>
    <associate|auto-13|<tuple|2.4|2>>
    <associate|auto-14|<tuple|3|?>>
    <associate|auto-15|<tuple|3.1|?>>
    <associate|auto-16|<tuple|A|?>>
    <associate|auto-2|<tuple|1.1|1>>
    <associate|auto-3|<tuple|1.2|1>>
    <associate|auto-4|<tuple|1.3|1>>
    <associate|auto-5|<tuple|2|2>>
    <associate|auto-6|<tuple|2.1|2>>
    <associate|auto-7|<tuple|2.1|2>>
    <associate|auto-8|<tuple|2.1|2>>
    <associate|auto-9|<tuple|4|2>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Definitions>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|1tab>|1.1<space|2spc>Representations
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2>>

      <with|par-left|<quote|1tab>|1.2<space|2spc>Numerics
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|1.3<space|2spc>Number data type
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Constructing
      a representation> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Defining representation
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>

      <with|par-left|<quote|4tab>|How do we sample from
      <with|mode|<quote|math>|<rigid|<with|mode|<quote|text>|<with|font-family|<quote|ss>|font-shape|<quote|right>|U>>>>?
      \V <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Error model. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Interval arithmetic. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>(Signed) permutation
      representations <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10>>

      <with|par-left|<quote|1tab>|2.3<space|2spc>From images
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Remarkable
      spaces> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Remarkable spaces
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13>>
    </associate>
  </collection>
</auxiliary>
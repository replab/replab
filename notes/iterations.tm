<TeXmacs|1.99.18>

<style|generic>

<\body>
  <doc-data|<doc-title|Iterative methods for representation
  theory>|<doc-author|<author-data|<author-name|Denis Rosset>>>>

  <section|Notation>

  <paragraph*|Scalars. \V>We write <math|i=<sqrt|-1>> the imaginary unit,
  <math|a+i*b\<in\><with|font|Bbb|C>>; we write a quaternion
  <math|q\<in\><with|font|Bbb|H>> as <math|q=x+y*j> with
  <math|x,y\<in\><with|font|Bbb|C>>, with the relations <math|k\<equiv\>i
  j=-j i>, <math|j k=i>, <math|k i=j>, <math|i<rsup|2>=j<rsup|2>=k<rsup|2>=-1>.

  The conjugate of a complex number <math|x=a+i*b\<in\><with|font|Bbb|C>> is
  <math|<overline|x>=a-i*b>; the conjugate of a quaternion <math|q=x+y j> is
  <math|<overline|q>=<overline|x>-y j>.

  <paragraph*|Linear algebra. \V>For a matrix <math|A> over
  <math|<with|font|Bbb|R>,<with|font|Bbb|C>> or <math|<with|font|Bbb|K>>, we
  write <math|A<rsup|<text|T>>> the transpose, <math|A<rsup|<text|C>>> the
  complex conjugate, <math|A<rsup|<text|H>>=<around*|(|A<rsup|<text|C>>|)><rsup|<text|T>>>
  the conjugate transpose.

  We write <math|A<rsup|-<text|T>>>, <math|A<rsup|-<text|C>>>,
  <math|A<rsup|-<text|H>>> the transpose, conjugate and conjugate transpose
  of <math|A<rsup|-1>>.

  <paragraph*|Representations. \V>We write
  <math|\<rho\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|K><rsup|D>|)>>
  a representation of the compact group <math|G> over the
  <math|d>-dimensional vector space over <math|<with|font|Bbb|K>=<with|font|Bbb|R>,<with|font|Bbb|C>>.
  We write <math|\<sigma\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|K><rsup|d>|)>>
  a subrepresentation of <math|\<rho\>> defined by the injection and
  projection linear maps <math|I> and <math|P>:

  <\equation>
    <label|Eq:SubRep>I:<with|font|Bbb|K><rsup|d>\<rightarrow\><with|font|Bbb|K><rsup|D>,<space|1em>P:<with|font|Bbb|K><rsup|D>\<rightarrow\><with|font|Bbb|K><rsup|d>,<space|2em>\<sigma\><rsub|g>\<equiv\>P*\<rho\><rsub|g>*I.
  </equation>

  The maps <math|I> and <math|P> can of course be written as matrices,
  equipping <math|<with|font|Bbb|K><rsup|D>,<with|font|Bbb|K><rsup|d>> with
  the Euclidean basis. We have <math|P\<cdot\>I=<with|font|Bbb|1>>.

  <paragraph*|Groups. \V>We write <math|<math-ss|GL><around*|(|<with|font|Bbb|K><rsup|D>|)>>
  the general linear group, <math|<math-ss|U><around*|(|<with|font|Bbb|K><rsup|D>|)>>
  the group of <math|d\<times\>d> matrices over <math|<with|font|Bbb|K>> with
  <math|M*M<rsup|<text|H>>=<with|font|Bbb|1>>. We can sample from the Haar
  measure of <math|<math-ss|U>> using various techniques.

  <paragraph*|Types of real irreducible representations. \V>Real irreducible
  representations can be of three different types:

  <\itemize-minus>
    <item>The representation is of real-type, and its Frobenius-Schur
    indicator is <math|1>.

    <item>The representation is of complex-type; its Frobenius-Schur
    indicator is 0 and it can be put in a basis where its images have the
    form:

    <\equation>
      <label|Eq:ComplexEncoding>\<rho\><rsub|g>=<matrix|<tformat|<table|<row|<cell|a<rsub|11>>|<cell|-b<rsub|11>>|<cell|a<rsub|21>>|<cell|-b<rsub|21>>|<cell|\<cdots\>>>|<row|<cell|b<rsub|11>>|<cell|a<rsub|11>>|<cell|b<rsub|21>>|<cell|a<rsub|21>>|<cell|>>|<row|<cell|\<vdots\>>|<cell|>|<cell|>|<cell|>|<cell|>>>>>
    </equation>

    <item>The representation is of quaternionic-type; its Frobenius-Schur
    indicator is <math|-1> and it can be put in a basis where its images have
    the form:

    <\equation>
      <label|Eq:QuaternionEncoding>\<rho\><rsub|g>=<matrix|<tformat|<cwith|1|4|1|4|cell-halign|r>|<table|<row|<cell|a<rsub|11>>|<cell|-b<rsub|11>>|<cell|-c<rsub|11>>|<cell|-d<rsub|11>>|<cell|\<cdots\>>>|<row|<cell|b<rsub|11>>|<cell|a<rsub|11>>|<cell|-d<rsub|11>>|<cell|c<rsub|11>>|<cell|>>|<row|<cell|c<rsub|11>>|<cell|d<rsub|11>>|<cell|a<rsub|11>>|<cell|-b<rsub|11>>|<cell|>>|<row|<cell|d<rsub|11>>|<cell|-c<rsub|11>>|<cell|b<rsub|11>>|<cell|a<rsub|11>>|<cell|>>|<row|<cell|\<vdots\>>|<cell|>|<cell|>|<cell|>|<cell|>>>>>
    </equation>
  </itemize-minus>

  <paragraph*|Encoding real irreducible representations. \V>We generalize
  slightly the definition<nbsp><eqref|Eq:SubRep> to encode the different
  representation types. In particular, we extend the choice of scalars to the
  quaternions, <math|<with|font|Bbb|L>=<with|font|Bbb|R>,<with|font|Bbb|C>,<with|font|Bbb|H>>,
  and write:

  <\equation>
    <label|Eq:SubRep1>I:<with|font|Bbb|L><rsup|d>\<rightarrow\><with|font|Bbb|L><rsup|D>,<space|1em>P:<with|font|Bbb|L><rsup|D>\<rightarrow\><with|font|Bbb|L><rsup|d>,<space|2em>\<sigma\><rsub|g>\<equiv\>P*<around*|(|\<rho\><rsub|g>|)><rsub|<with|font|Bbb|L>>*I.
  </equation>

  where <math|<around*|(|\<rho\><rsub|g>|)><rsub|<with|font|Bbb|L>>> is the
  embedding of <math|\<rho\><rsub|g>\<in\><with|font|Bbb|K><rsup|D\<times\>D>>
  in <math|<with|font|Bbb|L><rsup|D\<times\>D>>.

  Then complex-type real representations of dimension <math|2d> can be
  represented by injection/projection maps of dimension <math|d> over
  <math|<with|font|Bbb|L>=<with|font|Bbb|C>>. Same for quaternion-type real
  representations of dimension <math|4d> over
  <math|<with|font|Bbb|L>=<with|font|Bbb|H>>.

  <section|List of primitives>

  The following operations have an iterative algorithm. In the following, we
  let <math|\<rho\>> be a (pretty good approximation of a) representation of
  <math|G>.\ 

  <paragraph*|Find a projection map. \V>Given an injection map <math|I> and a
  representation <math|\<rho\>>, find a projection map <math|P> such that
  <math|\<sigma\><rsub|g>=P*\<rho\><rsub|g>I> is a subrepresentation.

  See Section<nbsp><reference|Sec:FindProjectionMap>.

  <paragraph*|Refining a subrepresentation. \V>Given an approximate
  subrepresentation <math|<wide|\<sigma\>|~>> defined using a representation
  <math|\<rho\>> and approximate injection and projection maps
  <math|<wide|I|~>>, <math|<wide|P|~>>, find a better approximate
  subrepresentation <math|\<sigma\>> according to approximate maps <math|I>,
  <math|P>.

  See Section<nbsp><reference|Sec:RefineNonUnitary> for the non-unitary
  version and Section<nbsp><reference|Sec:RefineUnitary> for the unitary
  version.

  <strong|Reveal the division algebra present in an irreducible real
  subrepresentation. \V> Given an approximate subrepresentation
  <math|<wide|\<sigma\>|~>> defined using
  <math|<around*|(|\<rho\>,<wide|I|~>,<wide|P|~>|)>>, find an approximate
  subrepresentation <math|\<sigma\>>, defined by
  <math|<around*|(|\<rho\>,I,P|)>> such that <math|\<sigma\><rsub|g>> has the
  proper form<nbsp><eqref|Eq:ComplexEncoding> or
  <eqref|Eq:QuaternionEncoding>.

  See Section<nbsp><reference|Sec:Reveal> for a sketch.

  <paragraph*|Harmonize a subrepresentation. \V>Given an approximate
  subrepresentation <math|<wide|\<sigma\>|~>> (defined by
  <math|<wide|I|~>,<wide|P|~>>) and a (pretty good approximation of a)
  subrepresentation <math|\<sigma\><rsub|0>> (defined by
  <math|I<rsub|0>,P<rsub|0>>), find an approximate subrepresentation
  <math|\<sigma\>> (defined by <math|I,P>) such that
  <math|\<sigma\><around*|(|g|)>\<sim\>\<sigma\><rsub|0><around*|(|g|)>>.

  <strong|Make a subrepresentation unitary. \V> Given an approximate
  subrepresentation <math|<wide|\<sigma\>|~>> (defined by
  <math|<wide|I|~>,<wide|P|~>>), find a subrepresentation <math|\<sigma\>>
  such that <math|\<sigma\><rsub|g<rsup|-1>>=<around*|(|\<sigma\><rsub|g>|)><rsup|<text|H>>>.

  <paragraph*|Split a reducible subrepresentation. \V>Given a reducible
  subrepresentation <math|<wide|\<sigma\>|~>> (defined by <math|<wide|I|~>>,
  <math|<wide|P|~>>), find two subrepresentations
  <math|\<sigma\><rsub|1>,\<sigma\><rsub|2>> defined by
  <math|<around*|(|I<rsub|1>,P<rsub|1>|)>> and
  <math|<around*|(|I<rsub|2>,P<rsub|2>|)>> such that

  <\equation>
    P<rsub|i>I<rsub|j>=\<delta\><rsub|i j><with|font|Bbb|1>.
  </equation>

  <section|Iterative algorithms>

  <subsection|General construction>

  In the algorithms below, we monitor the difference between two successive
  iterates as <math|\<delta\><rsub|k>=<around*|\<\|\|\>|X<rsub|k+1>-X<rsub|k>|\<\|\|\>><rsub|FRO>>,
  which, by analogy with other solvers we call the <em|step size>. We also
  measure an optimality criterion <math|\<omega\><rsub|k>>, which may not be
  completely correlated with proper optimality (this depends on the
  algorithm).

  However, we assume that both the optimality criterion and the step size
  will converge to <math|0> as <math|k\<rightarrow\>\<infty\>>.

  Due to the random nature of the algorithm, we consider the window of the
  last <math|R> samples:

  <\equation>
    <wide|\<omega\>|\<bar\>>=<around*|(|\<omega\><rsub|k-R>,\<ldots\>,\<omega\><rsub|k>|)>,<space|2em><wide|\<delta\>|\<bar\>>=<around*|(|\<delta\><rsub|k-R>,\<ldots\>,\<delta\><rsub|k>|)>.
  </equation>

  We then fit an exponential model of the form
  <math|a\<cdot\>10<rsup|b\<cdot\>k>> on those samples; this corresponds to a
  linear fit on <math|log<rsub|10><wide|\<omega\>|\<bar\>>> and
  <math|log<rsub|10><wide|\<delta\>|\<bar\>>>. However, we observe that
  outliers are present in the sequence. We thus use the Theil-Sen estimator
  on those sequences.

  Basically, the Theil-Sen estimator considers all pairwise slopes
  <math|\<sigma\><rsup|\<omega\>><rsub|<around*|(|i,j|)>>=<frac|log<rsub|10><wide|\<omega\>|\<bar\>><rsub|i>-log<rsub|10><wide|\<omega\>|\<bar\>><rsub|j>|i-j>>
  and returns the median value <math|<wide|\<sigma\>|^><around*|(|<wide|\<omega\>|\<bar\>>|)>>.

  We set a minimum and maximum number of iterations. We also set an
  optimality tolerance and a step size tolerance.

  Obviously, this minimum number is <math|\<geqslant\>R>.

  <\named-algorithm|Test>
    <strong|Input>: Optimality factors <math|<strong|\<omega\>>=<around*|(|\<omega\><rsub|i>|)><rsub|i>>,
    step sizes <math|<strong|\<delta\>>=<around*|(|\<delta\><rsub|i>|)><rsub|i>>,
    current step <math|k>

    <strong|Output>: Exit code, zero if continue, negative if error, positive
    if successful

    <\enumerate-numeric>
      <item><math|exitFlag\<assign\>0>

      <item><strong|if> <math|k\<geqslant\>R> <strong|and>
      <math|k\<geqslant\>minIterations>

      <item><space|1em>Compute the slopes
      <math|<wide|\<sigma\>|^><around*|(|<wide|\<omega\>|\<bar\>>|)>> and
      <math|<wide|\<sigma\>|^><around*|(|<wide|\<delta\>|\<bar\>>|)>>.

      <item><space|1em>Compute <math|max <wide|\<omega\>|\<bar\>>> and
      <math|max <wide|\<delta\>|\<bar\>>>

      <item><space|1em><math|maxSlope\<assign\>-1/maxIterations>

      <item><space|1em><strong|if> <math|max
      <wide|\<omega\>|\<bar\>>\<leqslant\>optimalityTol> <strong|and>
      <math|max <wide|\<delta\>|\<bar\>>\<leqslant\>stepSizeTol> <strong|and>
      <math|<wide|\<sigma\>|^><around*|(|<wide|\<omega\>|\<bar\>>|)>\<geqslant\>maxSlope>
      <strong|and> <math|<wide|\<sigma\>|^><around*|(|<wide|\<delta\>|\<bar\>>|)>\<geqslant\>maxSlope>

      <item><space|2em><math|exitFlag\<assign\>1>

      <item><space|1em><strong|endif>

      <item><strong|endif>

      <item><strong|if> <math|k\<geqslant\>maxIterations> <strong|and>
      <math|exitFlag=0>

      <item><space|1em><math|exitFlag\<assign\>-1>

      <item><strong|endif>
    </enumerate-numeric>
  </named-algorithm>

  <subsection|Building blocks>

  <paragraph*|Matrix inverse. \V>For a matrix <math|A>, the iteration

  <\equation>
    X<rsub|k+1>=2X<rsub|k>-X<rsub|k>A*X<rsub|k>
  </equation>

  converges to <math|A<rsup|-1>> provided the eigenvalues of
  <math|\<Delta\>=<with|font|Bbb|1>-A\<cdot\>X<rsub|0>> have magnitudes less
  than one<nbsp><cite|Schulz1933>. In turn, we have
  <math|<around*|\<\|\|\>|\<Delta\>*<wide|x|\<vect\>>|\<\|\|\>><rsub|2>\<leq\><around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|FRO><around*|\<\|\|\>|x|\<\|\|\>><rsub|2>>,
  so that the spectral radius of <math|\<Delta\>> is bounded by the Frobenius
  norm <math|<around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|FRO>>.

  <paragraph*|Matrix inverse square root. \V>For a matrix <math|A>, the
  following<nbsp><cite|Sherif1991> converges to <math|A<rsup|-1/2>>:

  <\equation>
    S<rsub|k+1>=2S<rsub|k><around*|(|<with|font|Bbb|1>+A*S<rsub|k>*S<rsub|k>|)><rsup|-1>.
  </equation>

  The only condition is that <math|A> does not have real negative
  eigenvalues.

  <paragraph*|Nearest orthogonal matrix. \V>See Wikipedia/Orthogonal matrix.

  <\equation>
    N<rsub|k>=Q<rsub|k><rsup|<text|H>>Q<rsub|k>,<space|2em>P<rsub|k>=<frac|1|2>Q<rsub|k>N<rsub|k>,<space|2em>Q<rsub|k+1>=2Q<rsub|k>+P<rsub|k>N<rsub|k>-3P<rsub|k>.
  </equation>

  Stable when the condition number of <math|Q<rsub|0>> is less than three?

  <paragraph*|Biorthogonalization step. \V>Given a subrepresentation defined
  by <math|<around*|(|\<rho\>,I<rsub|0>,P<rsub|0>|)>>, this updates the
  projection map so that (on average?), <math|<around*|\<\|\|\>|P\<cdot\>I<rsub|0>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>\<leqslant\><around*|\<\|\|\>|P<rsub|0>\<cdot\>I<rsub|0>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>.

  Note that this does not use the group/representation structure, rather
  adjusts the projection map.

  <\named-algorithm|BiorthogonalizationStep (projection version)>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P<rsub|0>\<in\><with|font|Bbb|L><rsup|d\<times\>D>> with
    <math|d\<leqslant\>D>

    <strong|Output>: <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>

    <\enumerate-numeric>
      <item><math|\<Pi\>\<assign\>P<rsub|0>\<cdot\>I<rsub|0>>

      <item><math|\<Delta\>\<assign\><with|font|Bbb|1><rsub|d>-\<Pi\>>

      <item><math|\<varepsilon\>=<around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|FRO>>

      <item><strong|if> <math|\<varepsilon\>\<less\>0.9>

      <item><space|1em><math|P\<assign\>P+\<Delta\>\<cdot\>P>

      <item><strong|else>

      <item><space|1em><math|P\<assign\>\<Pi\><rsup|-1>\<cdot\>P>

      <item><strong|end>
    </enumerate-numeric>
  </named-algorithm>

  <paragraph*|Biorthogonalization, idea behind construction. \V>When the
  convergence is not assured (<math|\<varepsilon\>> is a bound on the
  spectral radius of <math|P*I>), we perform a corrective step suggested
  in<nbsp><cite|Pan1989>.

  <subsection|Find a projection map for an injection
  map><label|Sec:FindProjectionMap>

  In this problem, we are given an injection map <math|I<rsub|0>> and look
  for a corresponding projection map.

  Note: In this algorithm as in the others, we try to avoid the explicit
  computation of images <math|\<rho\><rsub|g>>, but rather assume that these
  two actions can be efficiently computed:
  <math|\<rho\><rsub|g><rsup|row><around*|[|X|]>=\<rho\><rsub|g>X><math|> and
  <math|\<rho\><rsub|g><rsup|col><around*|[|X|]>=X*\<rho\><rsub|g<rsup|-1>>>.

  <\named-algorithm>
    FindProjectionMap
  <|named-algorithm>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>
    with <math|d\<leqslant\>D>

    <strong|Output>: <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>,
    <math|exitFlag>

    <\enumerate-numeric>
      <item>Sample <math|U> from <math|<math-ss|U><around*|(|<with|font|Bbb|L><rsup|d>|)>>

      <item><math|P\<assign\>U\<cdot\><around*|(|I<rsub|0>|)><rsup|<text|H>>>

      <item><math|P\<assign\>><strong|BiorthogonalizationStep><math|<around*|(|I,P|)>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em>Sample a multiset <math|S>,
      <math|<around*|\||S|\|>=N<rsub|samples>>, from the Haar measure of
      <math|G>

      <item><space|1em><math|P<rsub|1>:=<big|sum><rsub|g\<in\>S><around*|[|P\<cdot\><around*|(|\<rho\><rsub|g>\<cdot\>I<rsub|0>|)>|]>\<cdot\><around*|(|P\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>>

      <item><space|1em><math|\<delta\><rsub|k>\<assign\><around*|\<\|\|\>|P<rsub|1>-P|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|k>\<assign\>\<delta\><rsub|k>>

      <item><space|1em><math|P<rsub|1>\<assign\>><strong|BiorthogonalizationStep><math|<around*|(|I,P<rsub|1>|)>>

      <item><space|1em><strong|if> <strong|Test><math|<around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|P\<assign\>P<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  To be discussed later: stopping criterion
  <math|<strong|Test><around*|(|<strong|\<delta\>>,<strong|\<omega\>>|)>>.

  Note that in Step 6, we specify the associativity of the matrix products to
  optimize memory/CPU use.

  Why does this work?

  We start with an approximate <math|<wide|P|~>>, and look for <math|P> such
  that <math|P*I=<with|font|Bbb|1>> and the commutator
  <math|<around*|[|I\<cdot\>P,\<rho\><rsub|g>|]>=0> for all <math|g\<in\>G>.
  Write <math|<wide|P|~>=P+\<Delta\>P>, such ; we then compute
  <math|<wide|X|~>=<big|sum><rsub|g\<in\>S>\<rho\><rsub|g><around*|(|I\<cdot\><wide|P|~>|)>\<rho\><rsub|g<rsup|-1>>>
  which is an approximation of the commutant projection
  <math|X=<big|int>\<mathd\>\<mu\><around*|(|g|)>\<rho\><rsub|g><around*|(|I\<cdot\><wide|P|~>|)>\<rho\><rsub|g<rsup|-1>>=I\<cdot\>P>.
  If we had access to <math|X> then <math|<wide|P|~>*X=<around*|(|P\<cdot\>I+\<Delta\>P\<cdot\>I|)>P=P>
  which is what we look for. We believe we still improve <math|<wide|P|~>> by
  using <math|<wide|X|~>> instead of <math|X>, and take care of maintaing
  biorthogonality/proper scaling during the iterations.

  A formal proof is TBD.

  <subsection|Refine an approximate subrepresentation (non
  unitary)><label|Sec:RefineNonUnitary>

  Here, we are given an approximate representation and aim to improve the
  precision of its injection and projection maps.

  Optionally, we can provide injection/projection maps of previously refined
  subrepresentations in the same multiplicity space, as to preserve
  biorthogonality in the same isotypic component. Biorthogonality with
  inequivalent subrepresentations is not a problem, as this is catered for by
  the group averaging.

  We first present an algorithm based on projection into the commutant.

  <\named-algorithm>
    RefineNonUnitary (medium scale)
  <|named-algorithm>
    <strong|Input>: Approximate subrepresentation
    <math|<around*|(|\<rho\>,I<rsub|0>,P<rsub|0>|)>>,
    <math|\<rho\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|K><rsup|D>|)>>,
    <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P<rsub|0>\<in\><with|font|Bbb|L><rsup|d\<times\>D>>

    <strong|Input> (optional): Known part of the multiplicity space
    <math|<around*|(|\<rho\>,I<rsub|\<perp\>>,P<rsub|\<perp\>>|)>>,
    <math|I<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|D\<times\>e>>,
    <math|P<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|e\<times\>D>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>, <math|exitFlag>

    <\enumerate-numeric>
      <item><math|I:=I<rsub|0>>

      <item><math|P\<assign\>P<rsub|0>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em><math|<wide|F|~>\<assign\>I\<cdot\>P>

      <item><space|1em><math|<wide|F|\<bar\>>\<assign\><big|int>\<mathd\>\<mu\><around*|(|g|)>\<rho\><rsub|g><wide|F|~>\<rho\><rsub|g<rsup|-1>>>

      <item><space|1em><math|I<rsub|1>\<assign\>I>

      <item><space|1em><math|P<rsub|1>\<assign\>P>

      <item><space|1em><strong|for> <math|j=1,\<ldots\>,n<rsub|innerIterations>>

      <item><space|3em><math|I<rsub|1>\<assign\><wide|F|\<bar\>>\<cdot\>I<rsub|1>>

      <item><space|3em><math|P<rsub|1>\<assign\>P<rsub|1>\<cdot\><wide|F|\<bar\>>>

      <item><space|3em><strong|if> <math|e\<gtr\>0>

      <item><space|4em><math|I<rsub|1>\<assign\>I<rsub|1>-I<rsub|\<perp\>>\<cdot\><around*|(|P<rsub|\<perp\>>\<cdot\>I<rsub|1>|)>>

      <item><space|4em><math|P<rsub|1>\<assign\>P<rsub|1>-<around*|(|P<rsub|1>\<cdot\>I<rsub|\<perp\>>|)>\<cdot\>P<rsub|\<perp\>>>

      <item><space|3em><strong|endif>

      <item><space|3em><math|I<rsub|1>\<assign\>I<rsub|1>\<cdot\><around*|(|P<rsub|0>\<cdot\>I<rsub|1>|)><rsup|-1>>

      <item><space|3em><math|P<rsub|1>\<assign\><around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)><rsup|-1>\<cdot\>P<rsub|1>>

      <item><space|1em><strong|endfor>

      <item><space|1em><math|\<delta\><rsub|i>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>/<around*|\<\|\|\>|I<rsub|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|i>\<assign\><around*|\<\|\|\>|<wide|F|\<bar\>>-<wide|F|~>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|exitFlag\<assign\><strong|Test><around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|1em><strong|if> <math|exitFlag\<neq\>0>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><space|1em><math|P\<assign\>P<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  <\named-algorithm>
    RefineNonUnitary (large scale)
  <|named-algorithm>
    <strong|Input>: Approximate subrepresentation
    <math|<around*|(|\<rho\>,I<rsub|0>,P<rsub|0>|)>>,
    <math|\<rho\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|K><rsup|D>|)>>,
    <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P<rsub|0>\<in\><with|font|Bbb|L><rsup|d\<times\>D>>

    <strong|Input> (optional): Known part of the multiplicity space
    <math|<around*|(|\<rho\>,I<rsub|\<perp\>>,P<rsub|\<perp\>>|)>>,
    <math|I<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|D\<times\>e>>,
    <math|P<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|e\<times\>D>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>, <math|exitFlag>

    <\enumerate-numeric>
      <item><math|I:=I<rsub|0>>

      <item><math|P\<assign\>P<rsub|0>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em>Sample a multiset <math|S>,
      <math|<around*|\||S|\|>=N<rsub|samples>>, from the Haar measure of
      <math|G>

      <item><space|1em><math|I<rsub|1>\<assign\><rsub|><big|sum><rsub|g\<in\>S><around*|(|\<rho\><rsub|g>\<cdot\>I|)><around*|[|<around*|(|P\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>I|]>>

      <item><space|1em><math|P<rsub|1>:=<big|sum><rsub|g\<in\>S><around*|[|P\<cdot\><around*|(|\<rho\><rsub|g>\<cdot\>I|)>|]>\<cdot\><around*|(|P\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>>

      <item><space|1em><strong|if> <math|e\<gtr\>0>

      <item><space|2em><math|I<rsub|1>\<assign\>I<rsub|1>-I<rsub|\<perp\>>\<cdot\><around*|(|P<rsub|\<perp\>>\<cdot\>I<rsub|1>|)>>

      <item><space|2em><math|P<rsub|1>\<assign\>P<rsub|1>-<around*|(|P<rsub|1>\<cdot\>I<rsub|\<perp\>>|)>\<cdot\>P<rsub|\<perp\>>>

      <item><space|1em><strong|endif>

      <item><space|1em><math|f\<assign\>trace<around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)>/d>

      <item><space|1em><math|I<rsub|1>\<assign\>I<rsub|1>/<sqrt|<around*|\||f|\|>>>

      <item><space|1em><math|P<rsub|1>\<assign\>P<rsub|1>/<around*|(|<frac|f|<around*|\||f|\|>><sqrt|<around*|\||f|\|>>|)>>

      <item><space|1em><math|\<delta\><rsub|k>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>/<around*|\<\|\|\>|I<rsub|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|k>\<assign\><around*|\<\|\|\>|P<rsub|1>\<cdot\>I<rsub|1>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><strong|if> <math|\<omega\><rsub|k>\<less\>0.9>

      <item><space|2em><math|I<rsub|1>\<assign\>2I<rsub|1>-I<rsub|1>\<cdot\>P<rsub|0>\<cdot\>I<rsub|1>><space|1em>(approximation
      of <math|I<rsub|1>\<assign\>I<rsub|1>\<cdot\><around*|(|P<rsub|0>\<cdot\>I<rsub|1>|)><rsup|-1>>)

      <item><space|2em><math|P<rsub|1>\<assign\>2P<rsub|1>-P<rsub|1>\<cdot\>I<rsub|1>\<cdot\>P<rsub|1>>
      (approximation of <math|P<rsub|1>\<assign\><around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)><rsup|-1>\<cdot\>P<rsub|1>>)

      <item><space|1em><strong|else>

      <item><space|2em><math|I<rsub|1>\<assign\>I<rsub|1>\<cdot\><around*|(|P<rsub|0>\<cdot\>I<rsub|1>|)><rsup|-1>>

      <item><space|2em><math|P<rsub|1>\<assign\><around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)><rsup|-1>\<cdot\>P<rsub|1>>

      <item><space|1em><strong|end>

      <item><space|1em><math|exitFlag\<assign\><strong|Test><around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|1em><strong|if> <math|exitFlag\<neq\>0>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><space|1em><math|P\<assign\>P<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  <strong|TODO>: check that this algorithm works for reducible
  representations as well (in particular, the Newton approximation to the
  inverse could blow up).

  <subsection|Refine an approximate subrepresentation
  (unitary)><label|Sec:RefineUnitary>

  Here, we assume that <math|\<rho\>> is unitary, and that
  <math|P<rsub|0>=<around*|(|I<rsub|0>|)><rsup|<text|H>>>. We return
  <math|I>, understanding that <math|P=I<rsup|<text|H>>>.

  <\named-algorithm>
    RefineUnitary (medium scale)
  <|named-algorithm>
    <strong|Input>: Matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>

    <strong|Input> (optional): Known part of the multiplicity space
    <math|<around*|(|\<rho\>,I<rsub|\<perp\>>,P<rsub|\<perp\>>|)>>,
    <math|I<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|D\<times\>e>>,
    <math|P<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|e\<times\>D>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>, <math|exitFlag>

    <\enumerate-numeric>
      <item><strong|if> <math|I<rsub|0><rsup|<text|H>>\<cdot\>I<rsub|0>\<neq\><with|font|Bbb|1>>

      <item><space|1em><math|I<rsub|0>\<assign\><strong|qr><around*|(|I<rsub|0>|)>>

      <item><strong|endif>

      <item><math|I:=I<rsub|0>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em><math|<wide|F|~>\<assign\>I\<cdot\>I<rsup|<text|H>>>

      <item><space|1em><math|<wide|F|\<bar\>>\<assign\><big|int>\<mathd\>\<mu\><around*|(|g|)>\<rho\><rsub|g><wide|F|~>\<rho\><rsub|g<rsup|-1>>>

      <item><space|1em><math|I<rsub|1>\<assign\>I>

      <item><space|1em><strong|for> <math|j=1,\<ldots\>,n<rsub|innerIterations>>

      <item><space|3em><math|I<rsub|1>\<assign\><wide|F|\<bar\>>\<cdot\>I<rsub|1>>

      <item><space|3em><strong|if> <math|e\<gtr\>0>

      <item><space|4em><math|I<rsub|1>\<assign\>I<rsub|1>-I<rsub|\<perp\>>\<cdot\><around*|(|P<rsub|\<perp\>>\<cdot\>I<rsub|1>|)>>

      <item><space|3em><strong|endif>

      <item><space|3em><math|I<rsub|1>\<assign\><strong|qr><around*|(|I<rsub|1>|)>>

      <item><space|1em><strong|endfor>

      <item><space|1em><math|\<delta\><rsub|i>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>/<around*|\<\|\|\>|I<rsub|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|i>\<assign\><around*|\<\|\|\>|<wide|F|\<bar\>>-<wide|F|~>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|exitFlag\<assign\><strong|Test><around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|1em><strong|if> <math|exitFlag\<neq\>0>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  <\named-algorithm>
    RefineUnitary (large scale)
  <|named-algorithm>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>
    with <math|d\<leqslant\>D>

    <strong|Input> (optional): Known part of the multiplicity space given by
    <math|I<rsub|\<perp\>>>, <math|I<rsub|\<perp\>>\<in\><with|font|Bbb|L><rsup|D\<times\>e>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|exitFlag>

    <\enumerate-numeric>
      <item><math|I:=I<rsub|0>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em>Sample a multiset <math|S>,
      <math|<around*|\||S|\|>=N<rsub|samples>>, from the Haar measure of
      <math|G>

      <item><space|1em><math|I<rsub|1>\<assign\><rsub|><big|sum><rsub|g\<in\>S><around*|(|\<rho\><rsub|g>\<cdot\>I|)><around*|[|<around*|(|I<rsup|<text|H>>\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>I|]>><space|1em>noting
      that <math|<around*|(|I<rsup|<text|H>>\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>=<around*|(|\<rho\><rsub|g>\<cdot\>I|)><rsup|<text|H>>>

      <item><space|1em><strong|if> <math|e\<gtr\>0>

      <item><space|2em><math|I<rsub|1>\<assign\>I<rsub|1>-I<rsub|\<perp\>>\<cdot\><around*|(|<around*|(|I<rsub|\<perp\>>|)><rsup|<text|H>>\<cdot\>I<rsub|1>|)>>

      <item><space|1em><strong|endif>

      <item><space|1em><math|f\<assign\>trace<around*|(|I<rsub|1><rsup|<text|H>>\<cdot\>I<rsub|1>|)>/d>

      <item><space|1em><math|I<rsub|1>\<assign\>I<rsub|1>/<sqrt|f>>

      <item><space|1em><math|\<delta\><rsub|k>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>/<around*|\<\|\|\>|I|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|k>\<assign\><around*|\<\|\|\>|<around*|(|I<rsub|1><rsup|>|)><rsup|<text|H>>\<cdot\>I<rsub|1>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|N\<assign\><around*|(|I<rsub|1>|)><rsup|<text|H>>\<cdot\>I<rsub|1>><space|3em>In
      9-11, iteration to force unitarity

      <item><space|1em><math|R\<assign\>I<rsub|1>\<cdot\>N/2>

      <item><space|1em><math|I<rsub|1>\<assign\>2I<rsub|1>+R\<cdot\>N-3R>

      <item><space|1em><math|exitFlag\<assign\><strong|Test><around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|1em><strong|if> <math|exitFlag\<neq\>0>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  Note: the orthogonalization iteration is only stable for a condition number
  <math|\<less\>3>. We could upper bound the condition number based on the
  value of <math|\<omega\><rsub|k>> and use the QR algorithm to orthogonalize
  <math|I<rsub|1>> \U but we never observed stability problems.

  <subsection|Reveal the division algebra present in a
  representation><label|Sec:Reveal>

  Starting from the subrepresentation defined by <math|<around*|(|I,P|)>>,
  randomly partition the subrepresentation space
  <math|<with|font|Bbb|K><rsup|d>> in <math|<with|font|Bbb|K><rsup|d/2>\<oplus\><with|font|Bbb|K><rsup|d/2>>
  (complex-type) or <math|<with|font|Bbb|K><rsup|d/4>\<oplus\><with|font|Bbb|K><rsup|d/4>\<oplus\><with|font|Bbb|K><rsup|d/4>\<oplus\><with|font|Bbb|K><rsup|d/4>>
  and create a map in <math|<with|font|Bbb|L>=<with|font|Bbb|C>> or
  <math|<with|font|Bbb|L>=<with|font|Bbb|H>>. Use the refinement algorithms
  of Sections<nbsp><reference|Sec:RefineNonUnitary>
  and<nbsp><reference|Sec:RefineUnitary> on these maps.

  (This option has not been exercised much \U we currently use samples from
  the commutant).

  <subsection|Harmonize a subrepresentation (non-unitary)>

  The model subrepresentation <math|\<sigma\><rsub|M>> is defined using the
  representation <math|\<mu\>:G\<rightarrow\><math-ss|GL><around*|(|<with|font|Bbb|K><rsup|D<rsub|M>>|)>>,
  and the injection/projection maps <math|I<rsub|M>:<with|font|Bbb|L><rsup|d>\<rightarrow\><with|font|Bbb|L><rsup|D<rsub|M>>>,
  <math|P<rsub|M>:<with|font|Bbb|L><rsup|D<rsub|M>>\<rightarrow\><with|font|Bbb|L><rsup|d>>.

  We want to harmonize the subrepresentation
  <math|\<sigma\><around*|(|g|)>=P<rsub|0>*\<rho\><rsub|g>*I<rsub|0>> so that
  <math|\<sigma\><around*|(|g|)>\<sim\><around*|(|\<sigma\><rsub|M>|)><around*|(|g|)>=P<rsub|M>\<mu\><rsub|g>I<rsub|M>>.

  <\named-algorithm>
    HarmonizeNonUnitary
  <|named-algorithm>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P<rsub|0>\<in\><with|font|Bbb|L><rsup|d\<times\>D>>,
    <math|I<rsub|M>\<in\><with|font|Bbb|L><rsup|D<rsub|M>\<times\>d>>,
    <math|P<rsub|M>\<in\><with|font|Bbb|L><rsup|d\<times\>D<rsub|M>>> with
    <math|d\<leqslant\>D<rsub|M>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>

    <\enumerate-numeric>
      <item><math|I:=I<rsub|0>>

      <item><math|P\<assign\>P<rsub|0>>

      <item><strong|for> <math|i=1,\<ldots\>,maxIters>

      <item><space|1em>Sample two multisets
      <math|S<rsub|1>>,<math|S<rsub|2>>, <math|<around*|\||S<rsub|1>|\|>=<around*|\||S<rsub|2>|\|>=N<rsub|samples>>,
      from the Haar measure of <math|G>

      <item><space|1em><math|I<rsub|1>\<assign\><rsub|><big|sum><rsub|g\<in\>S<rsub|1>><around*|(|\<rho\><rsub|g>\<cdot\>I|)>\<cdot\><around*|[|<around*|(|P<rsub|M>\<cdot\>\<mu\><rsub|g<rsup|-1>>|)>\<cdot\>I<rsub|M>|]>>

      <item><space|1em><math|P<rsub|1>:=<big|sum><rsub|g\<in\>S<rsub|2>><around*|[|P<rsub|M>\<cdot\><around*|(|\<mu\><rsub|g>\<cdot\>I<rsub|M>|)>|]>\<cdot\><around*|(|P\<cdot\>\<rho\><rsub|g<rsup|-1>>|)>>

      <item><space|1em><math|f\<assign\>trace<around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)>/d>

      <item><space|1em><math|I<rsub|1>\<assign\>I<rsub|1>/<sqrt|<around*|\||f|\|>>>

      <item><space|1em><math|P<rsub|1>\<assign\>P<rsub|1>/<around*|(|<frac|f|<around*|\||f|\|>><sqrt|<around*|\||f|\|>>|)>>

      <item><space|1em><math|\<delta\><rsub|i>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|i>\<assign\><around*|\<\|\|\>|P<rsub|1>\<cdot\>I<rsub|1>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><strong|if> <strong|Test><math|<around*|(|<strong|\<delta\>>,<strong|\<omega\>>|)>>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><space|1em><math|P\<assign\>P<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  Note that we do not force the biorthogonalization
  <math|P\<cdot\>I=<with|font|Bbb|1>>, apart from the scalar factor; thus,
  <math|\<omega\><rsub|i>> has better statistical properties as we evolve
  <math|I> and <math|P> from a different set of samples.

  <subsection|Harmonize a subrepresentation (unitary)>

  Same spirit as above, but now we assume that <math|\<rho\>> is unitary, and
  that <math|<around*|(|I<rsub|0>|)><rsup|<text|H>>I<rsub|0>=<with|font|Bbb|1>>.
  The representation <math|\<mu\>> needs not be unitary, but
  <math|\<sigma\><rsub|M>> is.

  <\named-algorithm>
    HarmonizeUnitary
  <|named-algorithm>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|I<rsub|M>\<in\><with|font|Bbb|L><rsup|D<rsub|M>\<times\>d>>,
    <math|P<rsub|M>\<in\><with|font|Bbb|L><rsup|d\<times\>D<rsub|M>>> with
    <math|d\<leqslant\>D<rsub|M>>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>

    <\enumerate-numeric>
      <item>Sample <math|U> from <math|<math-ss|U><around*|(|<with|font|Bbb|L><rsup|d>|)>>

      <item><math|I:=I<rsub|0>*U>

      <item><strong|for> <math|i=1,\<ldots\>,maxIters>

      <item><space|1em>Sample a multiset <math|S>,
      <math|<around*|\||S|\|>=N<rsub|samples>>, from the Haar measure of
      <math|G>

      <item><space|1em><math|I<rsub|1>\<assign\><rsub|><big|sum><rsub|g\<in\>S<rsub|1>><around*|(|\<rho\><rsub|g>\<cdot\>I|)><around*|[|<around*|(|P<rsub|M>\<cdot\>\<mu\><rsub|g<rsup|-1>>|)>\<cdot\>I<rsub|M>|]>>

      <item><space|1em><math|f\<assign\>trace<around*|(|I<rsub|1><rsup|<text|H>>\<cdot\>I<rsub|1>|)>/d>

      <item><space|1em><math|I<rsub|1>\<assign\>I<rsub|1>/<sqrt|<around*|\||f|\|>>>

      <item><space|1em><math|\<delta\><rsub|i>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>>

      <item><space|1em><math|\<omega\><rsub|i>\<assign\><around*|\<\|\|\>|I<rsub|1><rsup|<text|H>>\<cdot\>I<rsub|1>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><strong|if> <strong|Test><math|<around*|(|<strong|\<delta\>>,<strong|\<omega\>>|)>>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><strong|endfor>

      <item><math|N\<assign\>I<rsup|<text|H>>\<cdot\>I><space|3em>In 15-17,
      iteration to finish off unitarity

      <item><math|R\<assign\>I\<cdot\>N/2>

      <item><math|I\<assign\>2I+R\<cdot\>N-3R>
    </enumerate-numeric>
  </named-algorithm>

  <subsection|Make a subrepresentation unitary>

  Note: we use <math|<around*|(|1+X|)>/2> and <math|<around*|(|3-X|)>/2> as
  approximations of the matrix square root and its inverse.

  <\named-algorithm>
    Unitarize
  <|named-algorithm>
    <strong|Input>: matrix <math|I<rsub|0>\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P<rsub|0>\<in\><with|font|Bbb|L><rsup|d\<times\>D>> with
    <math|d\<leqslant\>D>

    <strong|Output>: <math|I\<in\><with|font|Bbb|L><rsup|D\<times\>d>>,
    <math|P\<in\><with|font|Bbb|L><rsup|d\<times\>D>>

    <\enumerate-numeric>
      <item><math|I:=I<rsub|0>>

      <item><math|P\<assign\>P<rsub|0>>

      <item><strong|for> <math|k=1,\<ldots\>,maxIters>

      <item><space|1em>Sample a multiset <math|S>,
      <math|<around*|\||S|\|>=N<rsub|samples>>, from the Haar measure of
      <math|G>

      <item><space|1em><math|X:=<big|sum><rsub|g\<in\>S><around*|[|P\<cdot\><around*|(|\<rho\><rsub|g>\<cdot\>I|)>|]>\<cdot\><around*|[|P\<cdot\><around*|(|\<rho\><rsub|g>\<cdot\>I|)>|]><rsup|<text|H>>>

      <item><space|1em><math|f\<assign\>trace<around*|(|X|)>/d>

      <item><space|1em><math|X\<assign\>X/f>

      <item><space|1em><math|I<rsub|1>\<assign\>I<rsub|1>\<cdot\><frac|<with|font|Bbb|1>+X|2>>

      <item><space|1em><math|P<rsub|1>\<assign\><frac|3<with|font|Bbb|1>-X|2>\<cdot\>P<rsub|1>>

      <item><space|1em><math|\<omega\><rsub|k>\<assign\><around*|\<\|\|\>|P<rsub|1>\<cdot\>I<rsub|1>-<with|font|Bbb|1>|\<\|\|\>><rsub|FRO>>

      <item><space|1em><strong|if> <math|\<omega\><rsub|k>\<gtr\>0.9>

      <item><space|2em><math|P<rsub|1>\<assign\><around*|(|P<rsub|1>\<cdot\>I<rsub|1>|)><rsup|-1>\<cdot\>P<rsub|1>>

      <item><space|1em><strong|elseif> <math|k> is odd

      <item><space|2em><math|P<rsub|1>\<assign\>2P<rsub|1>-P<rsub|1>\<cdot\>I<rsub|1>\<cdot\>P<rsub|1>>

      <item><space|1em><strong|else> (<math|k> is even)

      <item><space|2em><math|I<rsub|1>\<assign\>2I<rsub|1>-I<rsub|1>\<cdot\>P<rsub|1>\<cdot\>I<rsub|1>>

      <item><space|1em><strong|end>

      <item><space|1em><math|\<delta\><rsub|k>\<assign\><around*|\<\|\|\>|I<rsub|1>-I|\<\|\|\>><rsub|FRO>>

      <item><space|1em><strong|if> <strong|Test><math|<around*|(|<strong|\<delta\>>,<strong|\<omega\>>,k|)>>

      <item><space|2em><strong|break>

      <item><space|1em><strong|endif>

      <item><space|1em><math|I\<assign\>I<rsub|1>>

      <item><space|1em><math|P\<assign\>P<rsub|1>>

      <item><strong|endfor>
    </enumerate-numeric>
  </named-algorithm>

  <\bibliography|bib|tm-plain|../../../Documents/Bibliography/Zotero>
    <\bib-list|3>
      <bibitem*|1><label|bib-Pan1989>V.<nbsp>Pan<localize| and >J.<nbsp>Reif.
      <newblock>Fast and efficient parallel solution of dense linear systems.
      <newblock><with|font-shape|italic|Computers & Mathematics with
      Applications>, 17(11):1481\U1491, jan 1989.<newblock>

      <bibitem*|2><label|bib-Schulz1933>Günther Schulz. <newblock>Iterative
      Berechung der reziproken Matrix. <newblock><with|font-shape|italic|ZAMM
      - Journal of Applied Mathematics and Mechanics / Zeitschrift für
      Angewandte Mathematik und Mechanik>, 13(1):57\U59, 1933.<newblock>

      <bibitem*|3><label|bib-Sherif1991>Nagwa Sherif. <newblock>On the
      computation of a matrix inverse square root.
      <newblock><with|font-shape|italic|Computing>, 46(4):295\U305, dec
      1991.<newblock>
    </bib-list>
  </bibliography>
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
    <associate|Eq:ComplexEncoding|<tuple|2|1>>
    <associate|Eq:QuaternionEncoding|<tuple|3|1>>
    <associate|Eq:SubRep|<tuple|1|1>>
    <associate|Eq:SubRep1|<tuple|4|2>>
    <associate|Sec:FindProjectionMap|<tuple|3.3|4>>
    <associate|Sec:RefineNonUnitary|<tuple|3.4|5>>
    <associate|Sec:RefineUnitary|<tuple|3.5|7>>
    <associate|Sec:Reveal|<tuple|3.6|8>>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-10|<tuple|2|2>>
    <associate|auto-11|<tuple|2|2>>
    <associate|auto-12|<tuple|2|2>>
    <associate|auto-13|<tuple|3|2>>
    <associate|auto-14|<tuple|3.1|2>>
    <associate|auto-15|<tuple|3.2|3>>
    <associate|auto-16|<tuple|3.2|3>>
    <associate|auto-17|<tuple|7|3>>
    <associate|auto-18|<tuple|8|3>>
    <associate|auto-19|<tuple|9|4>>
    <associate|auto-2|<tuple|1|1>>
    <associate|auto-20|<tuple|8|4>>
    <associate|auto-21|<tuple|3.3|4>>
    <associate|auto-22|<tuple|3.4|5>>
    <associate|auto-23|<tuple|3.5|7>>
    <associate|auto-24|<tuple|3.6|8>>
    <associate|auto-25|<tuple|3.7|8>>
    <associate|auto-26|<tuple|3.8|9>>
    <associate|auto-27|<tuple|3.9|10>>
    <associate|auto-28|<tuple|24|10>>
    <associate|auto-3|<tuple|1|1>>
    <associate|auto-4|<tuple|1|1>>
    <associate|auto-5|<tuple|1|1>>
    <associate|auto-6|<tuple|1|1>>
    <associate|auto-7|<tuple|3|2>>
    <associate|auto-8|<tuple|2|2>>
    <associate|auto-9|<tuple|2|2>>
    <associate|bib-Pan1989|<tuple|1|10>>
    <associate|bib-Schulz1933|<tuple|2|10>>
    <associate|bib-Sherif1991|<tuple|3|10>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|bib>
      Schulz1933

      Sherif1991

      Pan1989
    </associate>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Notation>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Scalars. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Linear algebra. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Representations. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Groups. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Types of real irreducible representations.
      \V <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Encoding real irreducible representations.
      \V <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-7><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>List
      of primitives> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-8><vspace|0.5fn>

      <with|par-left|<quote|4tab>|Find a projection map. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-9><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Refining a subrepresentation. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-10><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Harmonize a subrepresentation. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-11><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Split a reducible subrepresentation. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-12><vspace|0.15fn>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Iterative
      algorithms> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-13><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>General construction
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-14>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Building blocks
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-15>>

      <with|par-left|<quote|4tab>|Matrix inverse. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-16><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Matrix inverse square root. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-17><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Nearest orthogonal matrix. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-18><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Biorthogonalization step. \V
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-19><vspace|0.15fn>>

      <with|par-left|<quote|4tab>|Biorthogonalization, idea behind
      construction. \V <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-20><vspace|0.15fn>>

      <with|par-left|<quote|1tab>|3.3<space|2spc>Find a projection map for an
      injection map <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-21>>

      <with|par-left|<quote|1tab>|3.4<space|2spc>Refine an approximate
      subrepresentation (non unitary) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-22>>

      <with|par-left|<quote|1tab>|3.5<space|2spc>Refine an approximate
      subrepresentation (unitary) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-23>>

      <with|par-left|<quote|1tab>|3.6<space|2spc>Reveal the division algebra
      present in a representation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-24>>

      <with|par-left|<quote|1tab>|3.7<space|2spc>Harmonize a
      subrepresentation (non-unitary) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-25>>

      <with|par-left|<quote|1tab>|3.8<space|2spc>Harmonize a
      subrepresentation (unitary) <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-26>>

      <with|par-left|<quote|1tab>|3.9<space|2spc>Make a subrepresentation
      unitary <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-27>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|Bibliography>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-28><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
<TeXmacs|1.99.13>

<style|generic>

<\body>
  <doc-data|<doc-title|Testing for irreducibility>|<doc-author|<author-data|<author-name|Denis
  Rosset>>>>

  Let <math|G> be a compact group. Let <math|\<rho\>:G\<rightarrow\>\<cal-U\><around*|(|<with|font|Bbb*|F<rsup|d>>|)>>
  be a unitary representation over a field
  <math|<with|font|Bbb*|F>\<in\><around*|{|<with|font|Bbb*|R>,<with|font|Bbb*|C>|}>>;
  equipping <math|V=<with|font|Bbb*|F><rsup|d>> with the Euclidean basis, we
  can write the images <math|\<rho\><rsub|g>> using matrices.

  Our problem is to determine whether <math|\<rho\>> is irreducible over
  <math|<with|font|Bbb*|F>>. We use the following test.

  <section|The test>

  Let <math|X> be a matrix sample from the Gaussian Unitary/Orthogonal
  Ensemble.

  In the real case, we have <math|X<rsub|i
  j>\<sim\>\<cal-N\><around*|(|0,1|)>>. In the complex case, we have
  <math|X=R+<sqrt|-1>I> with <math|R<rsub|i
  j>\<sim\>\<cal-N\><around*|(|0,1|)>> and <math|I<rsub|i
  j>\<sim\>\<cal-N\><around*|(|0,1|)>>.

  Then we compute its projection on the space of matrices that commute with
  <math|\<rho\>>:

  <\equation>
    <wide|X|^>=<big|int><rsub|G>\<mathd\>\<mu\><around*|(|g|)>*\<rho\><rsub|g>X\<rho\><rsub|g><rsup|-1>.
  </equation>

  Finally, we set <math|\<mu\>=tr <wide|X|^>/d>, and
  <math|\<Delta\>=<around*|\<\|\|\>|<wide|X|^>-\<mu\><with|font|Bbb*|1><rsub|d>|\<\|\|\>><rsub|F>>.

  <\proposition>
    If <math|\<rho\>> is irreducible, we have <math|\<Delta\>=0>. If
    <math|\<rho\>> is reducible, then the probability of having
    <math|\<Delta\>\<leqslant\>\<delta\>> is upper bounded by
    <math|\<varepsilon\>>:

    <\equation>
      P<around*|(|\<Delta\>\<leqslant\>\<delta\>|)>\<leqslant\>\<varepsilon\>,<space|2em>\<varepsilon\><around*|(|\<delta\>|)>=<choice|<tformat|<table|<row|<cell|erf<around*|(|\<delta\>/<sqrt|2>|)>,>|<cell|<with|font|Bbb*|F>=<with|font|Bbb*|R>,>>|<row|<cell|1-exp<around*|[|-\<delta\><rsup|2>/2|]>,>|<cell|<with|font|Bbb*|F>=<with|font|Bbb*|C>.>>>>>
    </equation>
  </proposition>

  We also have <math|\<delta\>=<sqrt|2>*erf<rsup|-1>\<varepsilon\>> (real)
  and <math|\<delta\>=<sqrt|-2 log 1 - e>>.

  <section|Proof>

  <subsection|Simplified case>

  We first consider the case where <math|\<rho\>> is the sum of two
  inequivalent representations <math|\<sigma\><rsup|1>> and
  <math|\<sigma\><rsup|2>>. Then there exists a matrix <math|Q> such that

  <\equation>
    Q \<rho\><rsub|g> Q<rsup|\<dag\>>=<matrix|<tformat|<table|<row|<cell|\<sigma\><rsub|g><rsup|1>>|<cell|>>|<row|<cell|>|<cell|\<sigma\><rsub|g><rsup|2>>>>>>,
  </equation>

  where <math|\<sigma\><rsup|1>> and <math|\<sigma\><rsup|2>> are two
  representations of <math|G>. By Schur's lemma, <math|<wide|X|^>> has the
  form

  <\equation>
    <wide|X|^>=Q<below|<wide*|<matrix|<tformat|<table|<row|<cell|\<lambda\><rsub|1><with|font|Bbb*|1><rsub|d<rsub|1>>>|<cell|>>|<row|<cell|>|<cell|\<lambda\><rsub|2><with|font|Bbb*|1><rsub|d<rsub|2>>>>>>>|\<wide-underbrace\>>|\<Lambda\>>Q<rsup|\<dag\>>,
  </equation>

  where <math|d<rsub|1>=dim \<sigma\><rsub|1>> and <math|d<rsub|2>=dim
  \<sigma\><rsub|2>>, and <math|d=d<rsub|1>+d<rsub|2>>. Now, the matrix
  <math|<wide|X|^>> obtained during the projection is the one that commutes
  with <math|\<rho\>> and minimizes the Frobenius norm
  <math|<around*|\<\|\|\>|X-<overline|X>|\<\|\|\>><rsub|F>>, and thus

  <\equation>
    \<lambda\><rsub|1>=<frac|1|d<rsub|1>><big|sum><rsub|i=1><rsup|d<rsub|1>><wide|X|^><rsub|i
    i>,<space|2em>\<lambda\><rsub|2>=<frac|1|d<rsub|2>><big|sum><rsub|i=d<rsub|1>+1><rsup|d><wide|X|^><rsub|i
    i>.
  </equation>

  We have furthermore <math|\<mu\>=<frac|tr <wide|X|^>|d>=<frac|tr
  \<Lambda\>|d>=<frac|\<lambda\><rsub|1>d<rsub|1>+\<lambda\><rsub|2>d<rsub|2>|d>>.

  <strong|In the real case>, <math|\<lambda\><rsub|1>> is the average of
  <math|d<rsub|1>> independent variables of distribution
  <math|\<cal-N\><around*|(|0,1|)>>, and <math|\<lambda\><rsub|2>> is the
  average of <math|d<rsub|2>> independent variables of distribution
  <math|\<cal-N\><around*|(|0,1|)>>. Thus
  <math|\<lambda\><rsub|1>\<sim\>\<cal-N\><around*|(|0,1/<sqrt|d<rsub|1>>|)>>,
  <math|\<lambda\><rsub|2>\<sim\>\<cal-N\><around*|(|0,1/<sqrt|d<rsub|2>>|)>>,
  while

  <\equation>
    \<lambda\><rsub|1>-\<lambda\><rsub|2>\<sim\>\<cal-N\><around*|(|0,<sqrt|<frac|d|d<rsub|1>d<rsub|2>>>|)>,
  </equation>

  and

  <\equation>
    P<around*|(|<around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|>\<leqslant\>\<epsilon\>|)>=erf<around*|(|<frac|\<epsilon\>|<sqrt|2>\<sigma\>>|)>,<space|2em>\<sigma\>=<sqrt|d/d<rsub|1>d<rsub|2>>.
  </equation>

  <strong|In the complex case>, <math|\<lambda\><rsub|1,2>=\<lambda\><rsup|R><rsub|1,2>+<sqrt|-1>\<lambda\><rsub|1,2><rsup|I>>
  where <math|\<lambda\><rsub|1><rsup|R,I>\<sim\>\<cal-N\><around*|(|0,1/<sqrt|d<rsub|1>>|)>>
  and <math|\<lambda\><rsub|2><rsup|R,I>\<sim\>\<cal-N\><around*|(|0,1/<sqrt|d<rsub|2>>|)>>.
  Thus

  <\equation>
    Re \<lambda\><rsub|1>-\<lambda\><rsub|2>,Im
    \<lambda\><rsub|1>-\<lambda\><rsub|2>\<sim\>\<cal-N\><around*|(|0,<sqrt|<frac|d|d<rsub|1>d<rsub|2>>>|)>
  </equation>

  and

  <\equation>
    P<around*|(|<around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|>\<leqslant\>\<epsilon\>|)>=P<around*|(|*<around*|(|Re
    \<lambda\><rsub|1>-\<lambda\><rsub|2>|)><rsup|2>+<around*|(|Im
    \<lambda\><rsub|1>-\<lambda\><rsub|2>|)><rsup|2>\<leqslant\>\<epsilon\><rsup|2>|)>=P<around*|(|x<rsup|2>+y<rsup|2>\<leqslant\>\<epsilon\><rsup|2>/\<sigma\><rsup|2>|)>
  </equation>

  where <math|x,y\<sim\>\<cal-N\><around*|(|0,1|)>>. Now,
  <math|x<rsup|2>+y<rsup|2>> is distributed according to the
  <math|\<chi\><rsup|2>> distribution, which has CDF
  <math|f<around*|(|z|)>=1-e<rsup|-z/2>>, and thus

  <\equation>
    P<around*|(|<around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|>\<leqslant\>\<epsilon\>|)>=1-exp<around*|[|-\<epsilon\><rsup|2>/2\<sigma\><rsup|2>|]>.
  </equation>

  Let us come back to the main computation. In the case of irreducible
  <math|\<rho\>>, we have <math|\<lambda\><rsub|1>=\<lambda\><rsub|2>> by
  Schur's lemma; and then <math|\<Delta\>=0> follows. Now, in the case of
  reducible <math|\<rho\>>, we compute

  <\eqnarray>
    <tformat|<table|<row|<cell|\<Delta\>>|<cell|=>|<cell|<around*|\<\|\|\>|<wide|X|^>-\<mu\><with|font|Bbb*|1><rsub|d>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|Q
    \<Lambda\> Q<rsup|\<dag\>>-\<mu\> <with|font|Bbb*|1><rsub|d>*|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|<around*|\<\|\|\>|\<Lambda\>-\<mu\><with|font|Bbb*|1><rsub|d>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|<sqrt|d<rsub|1><around*|\||\<lambda\><rsub|1>-\<mu\>|\|><rsup|2>+d<rsub|2><around*|\||\<lambda\><rsub|2>-\<mu\>|\|><rsup|2>>>>|<row|<cell|>|<cell|=>|<cell|<sqrt|d<rsub|1><around*|\||d\<lambda\><rsub|1>-d<rsub|1>\<lambda\><rsub|1>-d<rsub|2>\<lambda\><rsub|2>|\|><rsup|2>+d<rsub|2><around*|\||d\<lambda\><rsub|2>-d<rsub|1>\<lambda\><rsub|1>-d<rsub|2>\<lambda\><rsub|2>|\|><rsup|2>>/d>>|<row|<cell|>|<cell|=>|<cell|<sqrt|d<rsub|1><around*|\||d<rsub|2>\<lambda\><rsub|1>-d<rsub|2>\<lambda\><rsub|2>|\|><rsup|2>+d<rsub|2><around*|\||d<rsub|1>\<lambda\><rsub|2>-d<rsub|1>\<lambda\><rsub|1><rsub|>|\|><rsup|2>>/d>>|<row|<cell|>|<cell|=>|<cell|<sqrt|d<rsub|1>d<rsub|2><rsup|2><around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|><rsup|2>+d<rsub|1><rsup|2>d<rsub|2><around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|><rsup|2>>/d>>|<row|<cell|>|<cell|=>|<cell|<sqrt|d<rsub|1>d<rsub|2><around*|(|d<rsub|1>+d<rsub|2>|)><around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|><rsup|2>>/d>>|<row|<cell|>|<cell|=>|<cell|<sqrt|<frac|d<rsub|1>d<rsub|2>|d>><around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|>=<frac|<around*|\||\<lambda\><rsub|1>-\<lambda\><rsub|2>|\|>|\<sigma\>>>>>>
  </eqnarray>

  Then <math|P<around*|(|\<Delta\>\<leqslant\>\<delta\>|)>=erf<around*|(|\<delta\>/<sqrt|2>|)>>
  in the real case, and <math|P<around*|(|\<Delta\>\<leqslant\>\<delta\>|)>=1-exp<around*|[|-\<delta\><rsup|2>|]>>
  in the complex case.

  <subsection|General case>

  In the case where <math|\<sigma\><rsup|1>> and <math|\<sigma\><rsup|2>> are
  equivalent, the off-diagonal blocks can now be non-zero, which cannot
  decrease <math|\<Delta\>>. Our bound thus still holds.

  In the case where <math|\<rho\>> is not unitary, we have
  <math|\<Delta\>=<around*|\<\|\|\>|<wide|X|^>-\<mu\><with|font|Bbb*|1><rsub|d>|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|A<around*|(|<wide|X|^><rsub|U>-\<mu\><with|font|Bbb*|1><rsub|d>|)>A<rsup|-1>|\<\|\|\>><rsub|F>>
  where <math|A> is a change of basis matrix such that
  <math|A<rsup|-1>\<rho\><rsub|g>A> is unitary; and <math|<wide|X|^><rsub|U>>
  is the projection on the commutant of that unitary representation. Let
  <math|\<Delta\><rsub|U>> be the random variable arising out of the test for
  unitary representations. Then <math|\<Delta\>\<leqslant\>cond<around*|(|A|)>\<Delta\><rsub|U>>.

  We have not considered formally the case where <math|\<rho\>> splits into
  more than two subrepresentations; we believe our bound still holds, as
  <math|\<delta\>> and <math|d> are all small.
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
    <associate|auto-2|<tuple|2|1>>
    <associate|auto-3|<tuple|2.1|1>>
    <associate|auto-4|<tuple|2.2|3>>
    <associate|auto-5|<tuple|3|3>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>The
      test> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Proof>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Simplified case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <with|par-left|<quote|1tab>|2.2<space|2spc>General case
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Handling
      errors> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
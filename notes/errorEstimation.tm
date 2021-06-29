<TeXmacs|1.99.13>

<style|generic>

<\body>
  <section|Unitary representations>

  Let <math|\<rho\>:G\<rightarrow\>\<cal-U\><around*|(|V|)>> be the unitary
  parent representation, which we assume to be exact, and
  <math|\<sigma\>:G\<rightarrow\>\<cal-U\><around*|(|W|)>> with
  <math|\<sigma\><rsub|g>=E<rsup|\<dag\>>\<rho\><rsub|g>E> be the
  subrepresentation in the RepLAB sense. We write <math|E:W\<rightarrow\>V>
  the injection map which is also an isometry
  (<math|E<rsup|\<dag\>>E=<with|font|Bbb*|1>>). We define <math|D=dim V> and
  <math|d=dim W>.

  We write <math|\<Sigma\><rsub|\<rho\>><around*|[|X|]>=<big|int><rsub|G>\<mathd\>\<mu\><around*|(|g|)>
  \<rho\><rsub|g>X\<rho\><rsub|g><rsup|-1>> the operator that projects X in
  the commutant.

  Let <math|A:W\<rightarrow\>V> be the approximate injection map; we assume
  that <math|A> is an isometry as well, <math|A<rsup|\<dag\>>A=<with|font|Bbb*|1>>.
  We define <math|<wide|\<sigma\>|~><rsub|g>=A<rsup|\<dag\>>\<rho\><rsub|g>A>
  which is an approximate representation.

  We write the projectors <math|\<Pi\><rsub|E>=E E<rsup|\<dag\>>> and
  <math|\<Pi\><rsub|A>=A A<rsup|\<dag\>>>; we write
  <math|<around*|\<\|\|\>|\<cdot\>|\<\|\|\>><rsub|2>> the operator norm
  (largest singular value) and <math|<around*|\<\|\|\>|\<cdot\>|\<\|\|\>><rsub|F>>
  the Frobenius norm (<math|<around*|\<\|\|\>|X|\<\|\|\>><rsub|F>=<sqrt|tr
  X<rsup|\<dag\>>X>=<sqrt|tr X X<rsup|\<dag\>>>>).

  We know that, with <math|\<Pi\><rsub|A>=\<Pi\><rsub|E>+\<Delta\>>,

  <\equation>
    \<delta\>=<around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|\<Pi\><rsub|A>-\<Pi\><rsub|E>|\<\|\|\>><rsub|F>\<leqslant\><sqrt|<frac|\<pi\><rsup|2>
    d|2>><around*|(|<frac|\<varepsilon\><rsub|\<Sigma\>>|1-\<varepsilon\><rsub|\<Sigma\>>>|)>,<space|2em>\<varepsilon\><rsub|\<Sigma\>>=<around*|\<\|\|\>|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>-\<Pi\><rsub|A>|\<\|\|\>><rsub|F>.
  </equation>

  Then, <math|<wide|\<sigma\>|~><rsub|g>> is a
  <math|\<varepsilon\>>-approximate representation, with

  <\equation>
    \<varepsilon\>=min<rsub|U> max<rsub|g><around*|\<\|\|\>|<wide|\<sigma\>|~><rsub|g>-U\<sigma\><rsub|g>U<rsup|\<dag\>>|\<\|\|\>><rsub|F>=min<rsub|U>
    max<rsub|g><around*|\<\|\|\>|A<rsup|\<dag\>>\<rho\><rsub|g>A-U*E<rsup|\<dag\>>\<rho\><rsub|g>E*U<rsup|\<dag\>>|\<\|\|\>><rsub|F>=min<rsub|V>
    max<rsub|g><around*|\<\|\|\>|\<Pi\><rsub|A>\<rho\><rsub|g>\<Pi\><rsub|A>-V
    \<Pi\><rsub|E>\<rho\><rsub|g>\<Pi\><rsub|E>V<rsup|\<dag\>>|\<\|\|\>><rsub|F>
  </equation>

  where <math|U> and <math|V> are unitary matrices. That means that

  <\eqnarray>
    <tformat|<table|<row|<cell|\<varepsilon\>>|<cell|\<leqslant\>>|<cell|max<rsub|g><around*|\<\|\|\>|\<Pi\><rsub|A>\<rho\><rsub|g>\<Pi\><rsub|A>-\<Pi\><rsub|E>\<rho\><rsub|g>\<Pi\><rsub|E>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|max<rsub|g><around*|\<\|\|\>|<around*|(|\<Pi\><rsub|E>+\<Delta\>|)>\<rho\><rsub|g><around*|(|\<Pi\><rsub|E>+\<Delta\>|)>-\<Pi\><rsub|E>\<rho\><rsub|g>\<Pi\><rsub|E>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|max<rsub|g><around*|\<\|\|\>|\<Delta\>\<rho\><rsub|g>\<Pi\><rsub|E>+\<Pi\><rsub|E>\<rho\><rsub|g>\<Delta\>+\<Delta\>\<rho\><rsub|g>\<Delta\>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|max<rsub|g><around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|F>*\<cdot\><around*|\<\|\|\>|\<rho\><rsub|g>\<Pi\><rsub|E>|\<\|\|\>><rsub|2>+<around*|\<\|\|\>|\<Pi\><rsub|E>\<rho\><rsub|g>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<rho\><rsub|g>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|\<Delta\>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|2\<delta\>+\<delta\><rsup|2>.>>>>
  </eqnarray>

  <section|Non-unitary representations>

  Let <math|\<rho\>:G\<rightarrow\>GL<around*|(|V|)>> a possibly non-unitary
  representation. Then there exists an isomorphism <math|T:U\<rightarrow\>V>
  and a unitary <math|\<rho\><rsub|U>:G\<rightarrow\>\<cal-U\><around*|(|U|)>>
  such that <math|\<rho\>=T*\<rho\><rsub|U>T<rsup|-1>>. Let
  <math|<around*|(|I,P|)>> be an exact injection-projection pair, and
  <math|<around*|(|<wide|I|~>,<wide|P|~>|)>> an approximate one.

  Then

  <\equation>
    \<sigma\>=P \<rho\>I=P T \<rho\><rsub|U>T<rsup|-1>I,<space|2em><wide|\<sigma\>|~>=<wide|P|~>\<rho\><wide|I|~>=<wide|P|~>T
    \<rho\><rsub|U>T<rsup|-1><wide|I|~>.
  </equation>

  Let <math|R> be a matrix such that <math|Q<rsup|\<dag\>>=R P T> and
  <math|Q=T<rsup|-1>I R<rsup|-1>> (<with|color|red|prove that such matrix
  exists>).

  Similarly, let <math|<wide|R|~>> be a matrix such that
  <math|<wide|Q|~><rsup|\<dag\>>=<wide|R|~><wide|P|~>T> and
  <math|<wide|Q|~>=T<rsup|-1><wide|I|~><wide|R|~><rsup|-1>>
  (<with|color|red|prove that such matrix exists>).

  <\equation>
    <around*|\<\|\|\>|<wide|\<sigma\>|~>-\<sigma\>|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|R<rsup|-1>Q<rsup|\<dag\>>\<rho\><rsub|U>Q
    R-<wide|R|~><rsup|-1><wide|Q|~><rsup|\<dag\>>\<rho\><rsub|U><wide|Q|~><wide|R|~>|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|<wide|R|~><rsup|-1><around*|[|Q<rsup|\<dag\>>\<rho\><rsub|U>Q-<wide|Q|~><rsup|\<dag\>>\<rho\><rsub|U><wide|Q|~>|]><wide|R|~>|\<\|\|\>><rsub|F>,
  </equation>

  the inequality being due to the fact that <math|R> and <math|<wide|R|~>>
  provide the optimal bound, while using <math|<wide|R|~>> for both may be
  suboptimal while using the exact representation <math|\<sigma\>> closest to
  <math|<wide|\<sigma\>|~>>.

  Now

  <\equation>
    <around*|\<\|\|\>|<wide|\<sigma\>|~>-\<sigma\>|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|<wide|R|~><rsup|-1><around*|[|Q<rsup|\<dag\>>\<rho\><rsub|U>Q-<wide|Q|~><rsup|\<dag\>>\<rho\><rsub|U><wide|Q|~>|]><wide|R|~>|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|<wide|R|~><rsup|-1>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|<wide|R|~>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|Q<rsup|\<dag\>>\<rho\><rsub|U>Q-<wide|Q|~><rsup|\<dag\>>\<rho\><rsub|U><wide|Q|~>|\<\|\|\>><rsub|F>,
  </equation>

  where the last term can be estimated using the unitary representation
  bound:

  <\equation>
    <around*|\<\|\|\>|Q<rsup|\<dag\>>\<rho\><rsub|U>Q-<wide|Q|~><rsup|\<dag\>>\<rho\><rsub|U><wide|Q|~>|\<\|\|\>><rsub|F>\<leqslant\><around*|[|2\<delta\><rsub|U>+\<delta\><rsup|2><rsub|U>|]>
  </equation>

  where

  <\equation>
    \<delta\><rsub|U>=<around*|\<\|\|\>|\<Pi\><rsub|A><rsup|<around*|(|U|)>>-\<Pi\><rsub|E><rsup|<around*|(|U|)>>|\<\|\|\>><rsub|F>\<leqslant\><sqrt|<frac|\<pi\><rsup|2>
    d|2>><around*|(|<frac|\<varepsilon\><rsub|U>|1-\<varepsilon\><rsub|U>>|)>
  </equation>

  with

  <\equation>
    \<varepsilon\><rsub|U>=<around*|\<\|\|\>|\<Sigma\><rsub|\<rho\><rsub|U>><around*|[|\<Pi\><rsub|A><rsup|<around*|(|U|)>>|]>-\<Pi\><rsub|A><rsup|<around*|(|U|)>>|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|T<rsup|-1><around*|(|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>-\<Pi\><rsub|A>|)>T|\<\|\|\>><rsub|F>.
  </equation>

  If only the condition number of <math|T> is known, then

  <\equation>
    \<varepsilon\><rsub|U>\<leqslant\><around*|\<\|\|\>|T|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|T<rsup|-1>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>-\<Pi\><rsub|A>|\<\|\|\>><rsub|F>=cond<around*|(|T|)><around*|\<\|\|\>|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>-\<Pi\><rsub|A>|\<\|\|\>><rsub|F>,
  </equation>

  with <math|cond<around*|(|X|)>=<around*|\<\|\|\>|X|\<\|\|\>><rsub|2><around*|\<\|\|\>|X<rsup|-1>|\<\|\|\>><rsub|2>>
  the condition number of <math|X>. We can now deduce the bound on
  <math|\<delta\><rsub|U>>, and then

  <\equation>
    <around*|\<\|\|\>|<wide|\<sigma\>|~>-\<sigma\>|\<\|\|\>><rsub|F>\<leqslant\>cond<around*|(|<wide|R|~>|)><around*|[|2\<delta\><rsub|U>+\<delta\><rsup|2><rsub|U>|]>.
  </equation>

  If the matrix <math|T> is known, then <math|<wide|R|~>> can be computed by
  the QR decomposition. If it isn't known, then a bound is given by
  <math|cond<around*|(|<wide|R|~>|)>\<leqslant\>cond<around*|(|<wide|P|~>|)>cond<around*|(|T|)>>
  (prove it).

  <section|Inexact parent representation>

  We assume now that <math|<wide|\<sigma\>|~><rsub|g>=<wide|P|~><wide|\<rho\>|~><rsub|g><wide|I|~>>
  where <math|<wide|\<rho\>|~><rsub|g>> is approximate too, with
  <math|<wide|\<rho\>|~><rsub|g>=\<rho\><rsub|g>+<around*|(|\<Delta\>\<rho\>|)><rsub|g>>.
  Then we expand:

  <\equation>
    <around*|\<\|\|\>|<wide|\<sigma\>|~>-\<sigma\>|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|<wide|P|~><around*|(|\<rho\>+\<Delta\>\<rho\>|)><wide|I|~>-P\<sigma\>I|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|<wide|P|~>\<Delta\>\<rho\><wide|I|~>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|<wide|P|~>\<rho\><wide|I|~>-P\<sigma\>I|\<\|\|\>><rsub|F>\<leqslant\>cond<around*|(|<wide|P|~>|)><around*|\<\|\|\>|<wide|\<rho\>|~>-\<rho\>|\<\|\|\>><rsub|F>+cond<around*|(|<wide|R|~>|)><around*|[|2\<delta\><rsub|U>+\<delta\><rsup|2><rsub|U>|]>
  </equation>

  where the second term is computed using the machinery above that assumes an
  exact <math|\<rho\>>.

  However, the projection <math|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>>
  will not be available. Instead, we obtain
  <math|<wide|\<Sigma\>|~>=\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>+\<Delta\><rsub|\<Sigma\>>>
  and an estimation <math|<around*|\<\|\|\>|\<Delta\><rsub|\<Sigma\>>|\<\|\|\>><rsub|F>>
  of the error in the projection. Then, we compute
  <math|\<varepsilon\><rsub|U>> as

  <\equation>
    \<varepsilon\><rsub|U>\<leqslant\>cond<around*|(|T|)><around*|\<\|\|\>|\<Sigma\><rsub|\<rho\>><around*|[|\<Pi\><rsub|A>|]>-\<Pi\><rsub|A>|\<\|\|\>><rsub|F>\<leqslant\>cond<around*|(|T|)>\<cdot\><around*|[|<around*|\<\|\|\>|<wide|\<Sigma\>|~>-\<Pi\><rsub|A>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<Delta\><rsub|\<Sigma\>>|\<\|\|\>><rsub|F>|]>.
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
    <associate|auto-2|<tuple|2|1|../../.TeXmacs/texts/scratch/no_name_8.tm>>
    <associate|auto-3|<tuple|3|?|../../.TeXmacs/texts/scratch/no_name_8.tm>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Unitary
      representations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Non-unitary
      representations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>
    </associate>
  </collection>
</auxiliary>
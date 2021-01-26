<TeXmacs|1.99.13>

<style|generic>

<\body>
  <doc-data|<doc-title|Estimating error on
  projections>|<doc-author|<author-data|<author-name|Denis
  Rosset>|<\author-affiliation>
    November 6, 2020
  </author-affiliation>>>>

  Let <math|<wide|\<rho\>|~>> and <math|<wide|\<sigma\>|~>> be two
  approximate representations of a finite group <math|G>, corresponding to
  the exact representations <math|\<rho\>> and <math|\<sigma\>>. We have
  <math|<around*|\<\|\|\>|<wide|\<rho\>|~><rsub|g>-\<rho\><rsub|g>|\<\|\|\>><rsub|FRO>\<leq\>\<varepsilon\><rsub|\<rho\>>>
  and <math|<around*|\<\|\|\>|<wide|\<sigma\>|~><rsub|g>-\<sigma\><rsub|g>|\<\|\|\>><rsub|FRO>\<leq\>\<varepsilon\><rsub|\<sigma\>>>,
  valid for all <math|g\<in\>G>.

  We also define the <em|condition number> <math|\<gamma\><rsub|\<rho\>>>,
  such that <math|<around*|\<\|\|\>|<wide|\<rho\>|~><rsub|g>|\<\|\|\>><rsub|2>\<cong\><around*|\<\|\|\>|\<rho\><rsub|g>|\<\|\|\>><rsub|2>\<leqslant\>\<gamma\><rsub|\<rho\>>>
  for all <math|g\<in\>G>. When <math|\<rho\><rsub|g>> is unitary,
  <math|\<gamma\><rsub|\<rho\>>=1>; otherwise, it is given by the condition
  number of the matrix <math|A> such that <math|A \<rho\><rsub|g>A<rsup|-1>>
  is a unitary representation (TODO: prove that
  <math|\<gamma\><rsub|\<rho\>>> does not depend on <math|A>).

  The group <math|G> has a decomposition as a cartesian product of sets
  <math|G=T<rsub|n>\<times\>\<cdots\>\<times\>T<rsub|1>>, with
  <math|<around*|\||G|\|>=<big|prod><rsub|i><around*|\||T<rsub|i>|\|>>. Given
  a matrix <math|X>, we compute the approximation projection
  <math|<wide|X|~>=<wide|X|~><rsub|n>> using the recursion:

  <\equation>
    <wide|X|~><rsub|i>=<big|sum><rsub|t\<in\>T<rsub|i>><wide|\<rho\>|~><rsub|t><wide|X|~><rsub|i-1><wide|\<sigma\>|~><rsub|t<rsup|-1>>,<space|2em><wide|X|~><rsub|0>=X.
  </equation>

  We now want to compute the error associated with <math|<wide|X|~><rsub|i>>,
  compared to the exact projection <math|X<rsub|i>=<big|sum><rsub|t>\<rho\><rsub|t>X<rsub|i-1>\<sigma\><rsub|t>>.
  We write <math|\<xi\><rsub|i>=<around*|\<\|\|\>|X<rsub|i>-<wide|X|~><rsub|i>|\<\|\|\>><rsub|FRO>>,
  and set <math|\<xi\><rsub|0>=0>. Now, we bound:

  <\eqnarray>
    <tformat|<table|<row|<cell|\<xi\><rsub|i>=<around*|\<\|\|\>|<wide|X|~><rsub|i>-X<rsub|i>|\<\|\|\>><rsub|F>>|<cell|\<cong\>>|<cell|<big|sum><rsub|t><around*|\<\|\|\>|<around*|(|\<rho\><rsub|t>-<wide|\<rho\>|~><rsub|t>|)>X<rsub|i-1>\<sigma\><rsub|t<rsup|-1>>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<rho\><rsub|t><around*|(|<wide|X|~><rsub|i-1>-X<rsub|i-1>|)>\<sigma\><rsub|t<rsup|-1>>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<rho\><rsub|t>X<rsub|i-1><around*|(|\<sigma\><rsub|t>-\<sigma\><rsub|t<rsup|-1>>|)>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|t><around*|\<\|\|\>|<around*|(|\<rho\><rsub|t>-<wide|\<rho\>|~><rsub|t>|)>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|X<rsub|i-1>\<sigma\><rsub|t<rsup|-1>>|\<\|\|\>><rsub|2>+<around*|\<\|\|\>|\<rho\><rsub|t>|\<\|\|\>><rsub|2>*\<cdot\><around*|\<\|\|\>|<wide|X|~><rsub|i-1>-X<rsub|i-1>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<sigma\><rsub|t<rsup|-1>>|\<\|\|\>><rsub|2>+>>|<row|<cell|>|<cell|>|<cell|+<around*|\<\|\|\>|\<rho\><rsub|t>X<rsub|i-1>|\<\|\|\>><rsub|2><around*|\<\|\|\>|\<sigma\><rsub|t>-\<sigma\><rsub|t<rsup|-1>>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|\<leqslant\>>|<cell|<big|sum><rsub|t>\<varepsilon\><rsub|\<rho\>>\<cdot\><around*|\<\|\|\>|X<rsub|i-1>|\<\|\|\>><rsub|2>\<cdot\>\<gamma\><rsub|\<sigma\>>+\<gamma\><rsub|\<rho\>>\<xi\><rsub|i-1>\<gamma\><rsub|\<sigma\>>+\<gamma\><rsub|\<rho\>>\<cdot\><around*|\<\|\|\>|X<rsub|i-1>|\<\|\|\>>\<cdot\>\<varepsilon\><rsub|\<sigma\>>.>>>>
  </eqnarray>

  Now, we assume <math|<around*|\<\|\|\>|X<rsub|i-1>|\<\|\|\>><rsub|2>\<leqslant\>\<mu\>=<around*|\<\|\|\>|X|\<\|\|\>><rsub|2>>
  (TODO: prove). Then

  <\equation>
    \<xi\><rsub|i>\<leqslant\><around*|\||T<rsub|i>|\|><around*|(|\<varepsilon\><rsub|\<rho\>>\<mu\>\<gamma\><rsub|\<sigma\>>+\<gamma\><rsub|\<rho\>>\<xi\><rsub|i-1>\<gamma\><rsub|\<sigma\>>+\<gamma\><rsub|\<rho\>>\<mu\>\<varepsilon\><rsub|\<sigma\>>|)>.
  </equation>

  This gives terrible error bounds, of the order of
  <math|\<varepsilon\><rsub|n>=10<rsup|-1>> for <math|G=S<rsub|6>> and both
  <math|\<rho\>> and <math|\<sigma\>> being the dimension 16 irreducible
  representation corresponding to the partition <math|6=3+2+1> in the
  standard Specht basis. While its coefficients are in
  <math|<around*|{|-1,0,+1|}>>, we assumed an error on the image coefficients
  of the order of the machine precision.

  The true error <math|<around*|\<\|\|\>|<wide|X|~>-X|\<\|\|\>><rsub|F>> is
  of the order of <math|10<rsup|-15>>. We had estimated
  <math|\<varepsilon\><rsub|\<rho\>>=\<varepsilon\><rsub|\<sigma\>>\<cong\>10<rsup|-14>>,
  which is pretty pessimistic due to the internal architecture; the true
  value is around <math|10<rsup|-16>>. The nonunitary of the representation
  is <math|\<gamma\><rsub|\<rho\>>=\<gamma\><rsub|\<sigma\>>\<cong\>9>.

  In the error computation, most of the error comes from the nonunitary in
  the <math|\<gamma\><rsub|\<rho\>>\<xi\><rsub|i-1>\<gamma\><rsub|\<sigma\>>>
  term. Replacing <math|\<gamma\><rsub|\<rho\>>\<xi\><rsub|i-1>\<gamma\><rsub|\<sigma\>>\<rightarrow\>\<xi\><rsub|i-1>>
  gives a bound <math|\<xi\><rsub|i>\<simeq\>10<rsup|-8>.>

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
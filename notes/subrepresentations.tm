<TeXmacs|1.99.13>

<style|generic>

<\body>
  <doc-data|<doc-title|Subrepresentations>|<doc-author|<author-data|<author-name|Denis
  Rosset>|<\author-affiliation>
    November 24, 2020
  </author-affiliation>>>>

  <section|Definitions>

  Let <math|V=<with|font|Bbb*|F><rsup|D>> be a vector space where
  <math|<with|font|Bbb*|F>=<with|font|Bbb*|R>,<with|font|Bbb*|C>>, and
  <math|\<rho\>> a representation of a group <math|G> such that

  <\equation>
    \<rho\>:G\<rightarrow\>GL<around*|(|V|)>.
  </equation>

  While <math|V> is the Euclidean space, we consider its inner product
  structure as somehow arbitrary. We do not assume that <math|\<rho\>> is
  unitary. Equivalently, we do not assume that the Euclidean inner product
  <math|<around*|(|\<cdot\>,\<cdot\>|)><rsub|E>> is <math|\<rho\>>-invariant:
  i.e. <math|<around*|(|\<rho\><rsub|g><wide|x|\<vect\>>,\<rho\><rsub|g><wide|y|\<vect\>>|)><rsub|E>\<neq\><around*|(|<wide|x|\<vect\>>,<wide|y|\<vect\>>|)><rsub|E>>
  in general. At any place we use the inner product, or the \Pconjugate
  transpose\Q, any other inner product/transpose map can be substituted.

  <\definition>
    For a map <math|X:V\<rightarrow\>V>, not necessarily equivariant, we
    define <math|\<Sigma\><rsub|\<rho\>><around*|[|X|]>=<big|int><rsub|G>\<mathd\>\<mu\><rsub|g>
    \<rho\><rsub|g>X\<rho\><rsub|g><rsup|-1>>. Consequently,
    <math|\<Sigma\><rsub|\<rho\>><around*|[|X|]>\<rho\><rsub|g>=\<rho\><rsub|g>\<Sigma\><rsub|\<rho\>><around*|[|X|]>>
    is equivariant.
  </definition>

  We come to our main construction.

  <\proposition>
    <label|Prop:IP>Let <math|W=<with|font|Bbb*|F><rsup|d>>. Let
    <math|<around*|(|I,P|)>> be a pair composed of the (injective)
    <em|injection map> <math|I:W\<rightarrow\>V> and <em|projection map>
    <math|P:V\<rightarrow\>W>, with <math|\<Pi\>\<equiv\>I P>. If

    <\equation>
      P I=<with|font|Bbb*|1>,<space|2em><text|and><space|2em>\<Pi\>\<rho\><rsub|g>=\<rho\><rsub|g>\<Pi\>,<space|1em>\<forall\>g\<in\>G,
    </equation>

    then <math|\<sigma\><rsub|g>\<equiv\>P \<rho\><rsub|g>I> defines a
    representation <math|\<sigma\>:G\<rightarrow\>GL<around*|(|W|)>> of
    <math|G>.
  </proposition>

  <\proof>
    We have <math|\<sigma\><rsub|g>\<sigma\><rsub|h>=P\<rho\><rsub|g>I
    P\<rho\><rsub|h>I=P\<rho\><rsub|g>\<rho\><rsub|h>I P I=P\<rho\><rsub|g
    h>I=\<sigma\><rsub|g h>>.
  </proof>

  Note that our requirements imply that the <math|\<Pi\><rsup|2>=I P I P=I
  P=\<Pi\>> is a projector. Let <math|V<rsub|I>=range<around*|(|I|)>>;
  observe that <math|\<Pi\>> is a projector on <math|V<rsub|I>>, and that
  <math|V<rsub|I>> is a <math|\<rho\>>-invariant subspace of <math|V>. Thus,
  <math|<around*|(|I,P|)>> also defines a subrepresentation of <math|\<rho\>>
  in <math|V> on the invariant subspace <math|V<rsub|I>.>

  <\remark>
    If the representation <math|\<rho\>> is unitary, we can assume that
    <math|I> is an isometry (<math|I<rsup|\<dag\>>I=<with|font|Bbb*|1>>) and
    take <math|P=I<rsup|\<dag\>>>; but we do not assume this below. Our
    derivation is better understood keeping <math|I> and <math|P> distinct.
  </remark>

  <\note>
    In the \Pfolklore\Q, representations are usually defined using a basis
    <math|<around*|{|<wide|v|\<vect\>><rsub|1>,\<ldots\>,<wide|v|\<vect\>><rsub|d>|}>>.
    By writing <math|<wide|v|\<vect\>><rsub|i>=I <wide|e|\<vect\>><rsub|i>>
    for the Euclidean basis vectors <math|<around*|{|<wide|e|\<vect\>><rsub|i>|}>>
    of <math|W>, the basis is given by the columns of <math|I> written in
    matrix form. However, the corresponding projection <math|P> is not
    necessarily unique, nor the projector <math|\<Pi\>=I P>. How come?

    It turns out that we did not equip our space <math|V> with an inner
    product; in the textbook treatment, that is one of the first steps:
    construct a <math|G>-invariant inner product. However, that product is
    not necessarily be unique. Indeed, consider the case of a trivial group
    <math|G=<around*|{|e|}>> and a trivial representation <math|\<rho\>> with
    <math|D=dim V\<gtr\>1>. Then any inner product on <math|V> is
    <math|G>-invariant. Now, consider a <math|d=dim V<rsub|I>=1> invariant
    subspace spanned by <math|<wide|v|\<vect\>><rsub|1>> with <math|I> fully
    specified by <math|I*<wide|e|\<vect\>><rsub|1>=<wide|v|\<vect\>><rsub|1>.>
    Then there are infinitely many maps <math|P:V\<rightarrow\>W> such that
    <math|P <wide|v|\<vect\>><rsub|1>=<wide|e|\<vect\>><rsub|1>>.
  </note>

  <section|Incomplete information>

  There are three ways to define a subrepresentation of
  <math|\<rho\>:G\<rightarrow\>GL<around*|(|V|)>>:

  <\enumerate-numeric>
    <item>Provide a subspace <math|V<rsub|I>\<subseteq\>V>.

    <item>Consider a projector <math|\<Pi\>:V\<rightarrow\>V> commuting with
    <math|\<rho\>>.

    <item>Consider a pair of maps <math|<around*|(|I,P|)>> satisfying the
    requirements of Proposition<nbsp><reference|Prop:IP>.
  </enumerate-numeric>

  Note that <math|V<rsub|I>=range \<Pi\>>, and <math|\<Pi\>=I P>, and thus
  the quantity of information increases from 1. to 3.

  Conversely, in <math|2.> the subrepresentation <math|\<sigma\><rsub|g>> is
  defined only up to a similarity transformation
  <math|\<sigma\><rsub|g>\<rightarrow\>T\<sigma\><rsub|g>T<rsup|-1>>. In 1.,
  the projector <math|\<Pi\>> is not uniquely defined in case of nontrivial
  multiplicities.

  <subsection|Reconstructing missing information>

  To recover <math|\<Pi\>> from <math|V<rsub|I>>, we have the following
  proposition.

  <\proposition>
    <label|Prop:Proj>Let <math|V<rsub|I>\<subseteq\>V> be a
    <math|\<rho\>>-invariant subspace. Then the map
    <math|\<Pi\>:V\<rightarrow\>V>

    <\equation>
      \<Pi\>=\<Sigma\><rsub|\<rho\>><around*|[|I<around*|(|I<rsup|\<dag\>>I|)><rsup|-1>I<rsup|\<dag\>>|]>
    </equation>

    is a projector with range <math|V<rsub|I>> which commutes with
    <math|\<rho\>>. We use an arbitrary inner product to compute a conjugate
    transpose <math|I<rsup|\<dag\>>:V\<rightarrow\>W>.
  </proposition>

  <\proof>
    We verify that <math|I<rsup|\<dag\>>> is surjective and that
    <math|I<rsup|\<dag\>>I> is not singular. Then
    <math|<wide|\<Pi\>|^>=I<around*|(|I<rsup|\<dag\>>I|)><rsup|-1>I<rsup|\<dag\>>>
    has range <math|V<rsub|I>> by construction and is a projector:

    <\equation>
      <wide|\<Pi\>|^><rsup|2>=I<around*|(|I<rsup|\<dag\>>I|)><rsup|-1>I<rsup|\<dag\>>I<around*|(|I<rsup|\<dag\>>I|)><rsup|-1>I<rsup|\<dag\>>=<wide|\<Pi\>|^>.
    </equation>

    Now <math|<wide|\<Pi\>|^>> does not commute with <math|\<rho\>>, but
    <math|\<Pi\>=\<Sigma\><rsub|\<rho\>><around*|[|<wide|\<Pi\>|^>|]>> does,
    and has for range <math|V<rsub|I>> (<with|color|red|prove this last
    thing>!).
  </proof>

  To reconstruct a pair <math|<around*|(|I,P|)>> from a projector
  <math|\<Pi\>>, we have the following proposition.

  <\proposition>
    Let <math|\<Pi\>:V\<rightarrow\>V> be a <math|G>-equivariant projector.
    Then there exists a pair <math|<around*|(|I,P|)>> that obeys the
    \ requirements of Proposition<nbsp><reference|Prop:IP> such that
    <math|\<Pi\>=I P>.
  </proposition>

  <\proof>
    Let <math|W=<with|font|Bbb*|F><rsup|d>> with <math|d> the rank of
    <math|\<Pi\>>. Construct the injection map <math|I:W\<rightarrow\>V>
    using either:

    <\itemize-minus>
      <item>a set of <math|d> linearly independent columns of the matrix
      defining <math|\<Pi\>:<with|font|Bbb*|F><rsup|D>\<rightarrow\><with|font|Bbb*|F><rsup|D>>,

      <item>sampling a random map <math|<wide|I|^>:W\<rightarrow\>V> such
      that <math|I=\<Pi\> <wide|I|^>> has range <math|V<rsub|I>>.
    </itemize-minus>

    Then there exists a unique <math|P> such that <math|\<Pi\>=I P>.
  </proof>

  For fun, consider Maschke's theorem.

  <\proposition>
    Let <math|V<rsub|I>> be a <math|G>-invariant subspace of <math|V>. Then
    there exists a <math|G>-invariant subspace <math|V<rsub|I<rprime|'>>>
    such that <math|V=V<rsub|I>\<oplus\>V<rsub|I<rprime|'>>>.
  </proposition>

  <\proof>
    Construct <math|\<Pi\>> with range <math|V<rsub|I>> using
    Proposition<nbsp><reference|Prop:Proj>. Let
    <math|\<Pi\><rprime|'>=<with|font|Bbb*|1>-\<Pi\>>. Then
    <math|\<Pi\><rprime|'>\<Pi\>=\<Pi\>\<Pi\><rprime|'>=0> and thus
    <math|V<rsub|I<rprime|'>>=range \<Pi\><rprime|'>> is a complement of
    <math|V<rsub|I>>.
  </proof>

  <section|Approximate subrepresentations>

  We will need the following.

  <\definition>
    The <em|condition number> <math|C<rsub|\<rho\>>> of a representation
    <math|\<rho\>> is given by

    <\equation>
      C<rsub|\<rho\>>=min<rsub|T> cond<around*|(|T|)>,<space|2em><text|s.t.
      >T\<rho\>T<rsup|-1><text| is unitary>.
    </equation>

    It is also an upper bound on the operator norm
    <math|<around*|\<\|\|\>|\<rho\><rsub|g>|\<\|\|\>><rsub|2>>.
  </definition>

  <subsection|Computing the error affecting a subrepresentation>

  In general, we do not have access to an exact pair
  <math|<around*|(|I,P|)>>, rather to an approximation
  <math|<around*|(|<wide|I|~>,<wide|P|~>|)>>. We shall neglect the floating
  point errors in the relation <math|<wide|P\<cdot\>|~><wide|I|~>=<with|font|Bbb*|1>>.
  We assume that the main source of error is that the relation
  <math|<wide|\<Pi\>|~>\<rho\><rsub|g>=\<rho\><rsub|g><wide|\<Pi\>|~>>, with
  <math|<wide|\<Pi\>|~>=<wide|I|~> <wide|P|~>> is satisfied only
  approximately.

  To construct an error model for <math|<around*|(|<wide|I|~>,<wide|P|~>|)>>,
  we write <math|<wide|I|~>=I+\<Delta\><rsub|I>> and
  <math|<wide|P|~>=P+\<Delta\><rsub|P>> where <math|<around*|(|I,P|)>>
  satisfies the conditions of Proposition<nbsp><reference|Prop:IP>. Thus,
  <math|I> and <math|<wide|I|~>> have the same dimensions, and our
  construction works only if an invariant subspace <math|V<rsub|I>> of
  dimension <math|d=rank<around*|(|<wide|I|~>|)>> exists.

  We assume <math|\<Delta\><rsub|P>\<cdot\>I=0> (or
  <math|<wide|P|~>\<cdot\>I=<with|font|Bbb*|1>>); otherwise, we write
  <math|A=<wide|P|~>\<cdot\>I>, and set <math|I\<rightarrow\>I\<cdot\>A<rsup|-1>>
  and <math|P\<rightarrow\>A\<cdot\>P>.

  We write the exact representation <math|\<sigma\><rsub|g>=P \<rho\><rsub|g>
  I> and the approximate representation <math|<wide|\<sigma\>|~><rsub|g>=<wide|P|~>\<rho\><rsub|g><wide|I|~>>.

  Then, we measure the error by the Frobenius norms of
  <math|<around*|\<\|\|\>|P\<Delta\><rsub|I>|\<\|\|\>><rsub|F>>,
  <math|<around*|\<\|\|\>|\<Delta\><rsub|P>|\<\|\|\>><rsub|F>> and
  <math|<around*|\<\|\|\>|\<Delta\><rsub|I>|\<\|\|\>><rsub|F>>. We compute:

  <\equation>
    \<varepsilon\>=<around*|\<\|\|\>|<wide|\<sigma\>|~><rsub|g>-\<sigma\><rsub|g>|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|<wide|P|~>\<rho\><rsub|g><wide|I|~>-P\<rho\><rsub|g>I|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|<around*|(|P+\<Delta\><rsub|P>|)>\<rho\><rsub|g><around*|(|I+\<Delta\><rsub|I>|)>-P\<rho\><rsub|g>I|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|\<Delta\><rsub|P>
    \<rho\><rsub|g> I+P \<rho\><rsub|g> \<Delta\><rsub|I>+\<Delta\><rsub|P>
    \<rho\><rsub|g> \<Delta\><rsub|I>|\<\|\|\>><rsub|F>.
  </equation>

  Now, <math|\<Delta\><rsub|P>\<rho\><rsub|g>I=\<Delta\><rsub|P>I\<sigma\><rsub|g>=0>.
  Then

  <\equation>
    \<varepsilon\>=<around*|\<\|\|\>|P \<rho\><rsub|g>\<Delta\><rsub|I>+\<Delta\><rsub|P>\<rho\><rsub|g>\<Delta\><rsub|I>|\<\|\|\>><rsub|F>\<leqslant\><around*|\<\|\|\>|\<sigma\><rsub|g>P\<Delta\><rsub|I>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<Delta\><rsub|P>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<rho\><rsub|g>\<Delta\><rsub|I>|\<\|\|\>><rsub|2>\<leqslant\><around*|\<\|\|\>|\<sigma\><rsub|g>|\<\|\|\>><rsub|2>\<cdot\><around*|\<\|\|\>|P\<Delta\><rsub|I>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<Delta\><rsub|P>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<Delta\><rsub|I>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<rho\><rsub|g>|\<\|\|\>><rsub|2>,
  </equation>

  and <math|\<varepsilon\>\<leqslant\>C<rsub|\<sigma\>><around*|\<\|\|\>|P*\<Delta\><rsub|I>|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|\<Delta\><rsub|P>|\<\|\|\>><rsub|F>\<cdot\><around*|\<\|\|\>|\<Delta\><rsub|I>|\<\|\|\>><rsub|F>\<cdot\>C<rsub|\<rho\>>>.

  Finally, <math|C<rsub|\<sigma\>>\<leqslant\>cond<around*|(|P|)>\<cdot\>C<rsub|\<rho\>>>,
  where <math|cond<around*|(|P|)>=cond<around*|(|I|)>\<cong\>cond<around*|(|<wide|P|~>|)>\<cong\>cond<around*|(|<wide|I|~>|)>>
  is the condition number of the matrices defining <math|<around*|(|I,P|)>>
  and <math|<around*|(|<wide|I|~>,<wide|P|~>|)>>.

  <subsection|Finding a pair <math|<around*|(|I,P|)>> close to an approximate
  pair <math|<around*|(|<wide|I|~>,<wide|P|~>|)>>>

  Given <math|<around*|(|<wide|I|~>,<wide|P|~>|)>>, we set
  <math|<wide|P|~><rsub|1>=<wide|P|~>> and
  <math|<wide|I|~><rsub|1>=<wide|I|~>>. We then perform the iterations below,
  starting with <math|k=1>. In the notation below,
  <math|V=V<rsub|I>\<oplus\>V<rsub|I<rprime|'>>\<oplus\>V<rsub|0>> is
  decomposed into invariant subspaces <math|V<rsub|I>>,
  <math|V<rsub|I<rprime|'>>> and <math|V<rsub|0>>.

  We assume that <math|V<rsub|I>> is the invariant subspace we are looking
  for; that <math|V<rsub|0>> is an invariant subspace containing additional
  multiplicity spaces for the irreducible representations present in
  <math|V<rsub|I>> (we can always choose <math|V<rsub|I>> and
  <math|V<rsub|0>> such that <math|<wide|P|~>> does not have support in
  <math|V<rsub|0>>); that <math|V<rsub|I<rprime|'>>> is an invariant subspace
  whose irreducible representations do not overlap with <math|V<rsub|I>>.

  <\enumerate-numeric>
    <item>We decompose <math|<wide|I|~><rsub|k>=I<rsub|k>+\<Delta\><rsub|I<rsub|k>>>,
    for unknown <math|I<rsub|k>> and <math|\<Delta\><rsub|I<rsub|k>>> with
    <math|range I<rsub|k>=V<rsub|I>> and <math|range
    \<Delta\>I<rsub|k>=V<rsub|I<rprime|'>>>.

    <item>We decompose <math|<wide|P|~><rsub|k>=P<rsub|k>+\<Delta\><rsub|P<rsub|k>>>,
    for unknown <math|P<rsub|k>> and <math|\<Delta\><rsub|P<rsub|k>>> such
    that <math|P<rsub|k>\<cdot\>\<Delta\>I<rsub|k>=0> and
    <math|\<Delta\>P<rsub|k>\<cdot\>I<rsub|k>=0.>

    (We do not necessarily have <math|P<rsub|k>I<rsub|k>=<with|font|Bbb*|1>>.)

    <item>We compute <math|<wide|F|~><rsub|k>=<wide|I|~><rsub|k><wide|P|~><rsub|k>=I<rsub|k>P<rsub|k>+I<rsub|k>\<Delta\><rsub|P<rsub|k>>+\<Delta\><rsub|I<rsub|k>>P<rsub|k>+\<Delta\><rsub|I<rsub|k>>\<Delta\><rsub|P<rsub|k>>>.

    <item>We compute <math|<overline|F<rsub|k>>=\<Sigma\><rsub|\<rho\>><around*|(|<wide|F|~><rsub|k>|)>=I<rsub|k>P<rsub|k>+\<Delta\><rsub|I<rsub|k>>\<Delta\><rsub|P<rsub|k>>>
    (by Schur's lemma).

    <item>We monitor the norm <math|<around*|\<\|\|\>|<wide|F|~><rsub|k>-<overline|F<rsub|k>>|\<\|\|\>><rsub|F>>,
    which should decrease until machine precision is reached.

    <item>We compute, by solving a linear least-squares problem:

    <\eqnarray>
      <tformat|<table|<row|<cell|I<rprime|'><rsub|k>>|<cell|=>|<cell|<below|<text|argmin>|I<rsub|k><rprime|'>>
      <around*|\<\|\|\>|<overline|F<rsub|k>>-I<rsub|k><rprime|'><wide|P|~><rsub|k>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|<below|<text|argmin>|I<rsub|k><rprime|'>><around*|\<\|\|\>|I<rsub|k>P<rsub|k>+\<Delta\><rsub|I<rsub|k>>\<Delta\><rsub|P<rsub|k>>-I<rsub|k><rprime|'><around*|(|P<rsub|k>+\<Delta\><rsub|P<rsub|k>>|)>|\<\|\|\>><rsub|F>>>|<row|<cell|>|<cell|=>|<cell|<below|<text|argmin>|I<rsub|k><rprime|'>><around*|\<\|\|\>|<around*|(|I<rsub|k>-I<rsub|k><rprime|'>|)>P<rsub|k>+<around*|(|\<Delta\><rsub|I<rsub|k>>-I<rsub|k><rprime|'>|)>\<Delta\><rsub|P<rsub|k>>|\<\|\|\>><rsub|F>>>>>
    </eqnarray>

    and if <math|P<rsub|k>\<gg\>\<Delta\>P<rsub|k>>, then
    <math|I<rsub|k><rprime|'>\<approx\>I<rsub|k>>.

    <item>Similarly, we compute <math|><math|P<rsub|k><rprime|'>=<below|<text|argmin>|P<rsub|k><rprime|'>>
    <around*|\<\|\|\>|<overline|F<rsub|k>>-<wide|I|~><rsub|k>P<rprime|'><rsub|k>|\<\|\|\>><rsub|F>>,
    and <math|P<rsub|k><rprime|'>\<approx\>P<rsub|k>>.

    <item>We set <math|<wide|I|~><rsub|k+1>=I<rsub|k><rprime|'><around*|(|<wide|P|~>*\<cdot\>I<rprime|'><rsub|k>|)><rsup|-1>>
    to have <math|<wide|P|~>\<cdot\><wide|I|~><rsub|k+1>=<with|font|Bbb*|1>>,
    where <math|<wide|P|~>> is the original approximate projection.

    <item>We set <math|<wide|P|~><rsub|k+1>=<around*|(|P<rprime|'><rsub|k>\<cdot\><wide|I|~><rsub|k+1>|)>P<rsub|k><rprime|'>>
    to have <math|<wide|P|~><rsub|k+1><wide|I|~><rsub|k+1>=<with|font|Bbb*|1>>.

    <item>We monitor <math|<around*|\<\|\|\>|<wide|I|~><rsub|k+1>-<wide|I|~><rsub|k>|\<\|\|\>><rsub|F>>
    and <math|<around*|\<\|\|\>|<wide|P|~><rsub|k+1>-<wide|P|~><rsub|k>|\<\|\|\>><rsub|F>>
    to stop the iterations. (Stopping criterion to be fixed)
  </enumerate-numeric>

  Question: can we repeat steps 6-7 without recomputing
  <math|<overline|F<rsub|k>>>? Use iterative refinement?

  <subsection|A faster algorithm, unitary version>

  Let <math|<wide|Q|~><rsub|1>:W\<rightarrow\>V> be an isometry whose range
  approximates a subspace of <math|V> invariant under <math|\<rho\>>, and
  <math|\<rho\>> be unitary. Let <math|K\<geqslant\>0> be an integer. One
  step of the refinement process works as follows.

  <\enumerate-numeric>
    <item>For the step <math|k>:

    <item>If <math|K=0>, compute <math|<overline|F<rsub|k>>=\<Sigma\><rsub|\<rho\>><around*|(|<wide|Q|~><rsub|k>
    <wide|Q|~><rsup|\<dag\>><rsub|k>|)>>. If <math|K\<geqslant\>1>, compute
    <math|<overline|F<rsub|k>>=<around*|(|<big|sum><rsub|i=1><rsup|K>\<rho\><rsub|g<rsub|i>><wide|Q|~><rsub|k><wide|Q|~><rsub|k><rsup|\<dag\>>\<rho\><rsub|g<rsub|i>><rsup|-1>|)>/K>,
    where <math|<around*|{|g<rsub|i>|}>> are sampled uniformly randomly from
    <math|G>.

    <item>Let <math|X> minimize <math|<around*|\<\|\|\>|
    <wide|Q|~><rsub|k>Y<rsub|k><rsup|\<dag\>>-<overline|F<rsub|k>>|\<\|\|\>><rsub|2>>.
    The solution is <math|Y<rsub|k><rsup|\<dag\>>=<wide|Q|~><rsub|k><rsup|\<dag\>><overline|F<rsub|k>>>.
    Set <math|<wide|Q|~><rsub|k+1>> to an isometry that has the same range as
    <math|Y<rsub|k>> (using one of Gram-Schmidt, QR decomposition, SVD).

    <item>Either repeat steps 1. and 2.; or repeat only step 2 a few times by
    incrementing <math|k> but reusing <math|<overline|F<rsub|k+1>>=<overline|F<rsub|k>>.>
  </enumerate-numeric>

  Finally, remark that the product <math|Y<rsub|k>=<overline|F<rsub|k>><wide|Q|~><rsub|k>>
  can be computed faster, for <math|k\<geqslant\>1>, as:

  <\equation>
    Y<rsub|k>=<big|sum><rsub|i>\<rho\><rsub|g<rsub|i>><wide|Q|~><rsub|k><wide|Q|~><rsub|k><rsup|\<dag\>>\<rho\><rsub|g<rsub|i>><rsup|-1><wide|Q|~><rsub|k>=<big|sum><rsub|i><below|<wide*|<around*|(|\<rho\><rsub|g<rsub|i>><wide|Q|~><rsub|k>|)>|\<wide-underbrace\>>|n\<times\>m><below|<wide*|<around*|(|<wide|Q|~><rsub|k><rsup|\<dag\>>\<rho\><rsub|g<rsub|i>><rsup|-1><wide|Q|~><rsub|k>|)>|\<wide-underbrace\>>|m\<times\>m>,
  </equation>

  assuming that we have a fast method to compute products such as
  <math|A\<rho\><rsub|g>B> and <math|\<rho\><rsub|g>B>. This is the case for
  permutation representations, tensor products, etc. In that case, we perform
  step 2. only once (there is no point in reusing the same samples).

  \;

  <\note>
    Now, we prove that the solution of <math|min<around*|\<\|\|\>|A
    X-B|\<\|\|\>><rsub|F>> is <math|X=A<rsup|\<dag\>>B > if <math|A> is an
    isometry <math|<around*|(|A<rsup|\<dag\>>A=<with|font|Bbb*|1>|)>>. Let
    <math|A<rsup|\<perp\>>> be the matrix orthogonal to <math|A> such that
    the concatenation <math|<around*|(|A,A<rsup|\<perp\>>|)>> is a unitary
    matrix.

    Then

    <\equation>
      <around*|\<\|\|\>|A X- B|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|<around*|(|A,A<rsup|\<perp\>>|)><matrix|<tformat|<table|<row|<cell|X>>|<row|<cell|0>>>>>-B|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|<around*|(|A,A<rsup|\<perp\>>|)><matrix|<tformat|<table|<row|<cell|X>>|<row|<cell|0>>>>>-<around*|(|A,A<rsup|\<perp\>>|)><matrix|<tformat|<table|<row|<cell|A<rsup|\<dag\>>>>|<row|<cell|<around*|(|A<rsup|\<perp\>>|)><rsup|\<dag\>>>>>>>B|\<\|\|\>><rsub|F>
    </equation>

    Now, as the Frobenius norm is invariant under unitary transformations:

    <\equation>
      min <around*|\<\|\|\>|A X- B|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|<matrix|<tformat|<table|<row|<cell|X>>|<row|<cell|0>>>>>-<matrix|<tformat|<table|<row|<cell|A<rsup|\<dag\>>>>|<row|<cell|<around*|(|A<rsup|\<perp\>>|)><rsup|\<dag\>>>>>>>B|\<\|\|\>><rsub|F>=<around*|\<\|\|\>|X-A<rsup|\<dag\>>B|\<\|\|\>><rsub|F>+<around*|\<\|\|\>|<around*|(|A<rsup|\<perp\>>|)><rsup|\<dag\>>B|\<\|\|\>><rsub|F>.
    </equation>

    Thus <math|X=A<rsup|\<dag\>>B>.
  </note>

  <subsection|Faster algorithm, non unitary version>

  Let <math|<around*|(|<wide|I|~>,<wide|P|~>|)>> be an approximate pair.

  <\enumerate-numeric>
    <item>Set <math|I<rsub|1>=<wide|I|~>>,<math|P<rsub|1>=<wide|P|~>>. Set
    <math|k\<leftarrow\>1>.

    <item>Compute temporaries <math|J<rsub|k>=<big|sum><rsub|i=1><rsup|K><around*|(|\<rho\><rsub|g<rsub|i>>I<rsub|k>|)><around*|(|P<rsub|k>\<rho\><rsub|g<rsub|i>><rsup|-1>I<rsub|k>|)>>
    and <math|Q<rsub|k>=<big|sum><rsub|i=1><rsup|K><around*|(|P<rsub|k>\<rho\><rsub|g<rsub|i>>I<rsub|k>|)><around*|(|\<rho\><rsub|g<rsub|i>><rsup|-1>P<rsub|k>|)>>.

    <item>Set <math|I<rsub|k+1>=J<rsub|k><around*|(|<wide|P|~>J<rsub|k>|)><rsup|-1>>,
    to recover <math|\<Delta\>P\<cdot\>I<rsub|k+1>=0>.

    <item>Set <math|P<rsub|k+1>=<around*|(|Q<rsub|k>I<rsub|k+1>|)><rsup|-1>Q<rsub|k>>
    to recover <math|P<rsub|k+1>I<rsub|k+1>=<with|font|Bbb*|1>>.

    <item><math|k\<leftarrow\>k+1> and loop.
  </enumerate-numeric>

  <subsection|Error estimation>

  We cannot check <math|<around*|\<\|\|\>|<wide|F|~><rsub|k>-<overline|F<rsub|k>>|\<\|\|\>><rsub|F>>
  anymore as we only approximate it. However, we can check how the subspaces
  spanned by <math|<wide|Q|~><rsub|k>> and <math|<wide|Q|~><rsub|k+1>> are
  aligned. We construct the matrix <math|\<Lambda\>>:

  <\equation>
    \<Lambda\>=<wide|Q|~><rsub|k><rsup|\<dag\>><wide|Q|~><rsub|k+1>.
  </equation>

  Then, if <math|<wide|Q|~><rsub|k>> and <math|<wide|Q|~><rsub|k+1>> span the
  same space, they will be related by a unitary transformation, and the
  singular values of <math|\<Lambda\>> will all be one. We thus take

  <\equation>
    \<delta\>=2<around*|\<\|\|\>|<wide|Q|~><rsub|k><rsup|\<dag\>><wide|Q|~><rsub|k+1><wide|Q|~><rsub|k+1><rsup|\<dag\>><wide|Q|~><rsub|k>-<with|font|Bbb*|1>|\<\|\|\>><rsub|F>
  </equation>

  to compute <math|<sqrt|<big|sum><rsub|i><around*|(|\<lambda\><rsub|i>-1|)><rsup|2>>>,
  where the <math|\<lambda\><rsub|1>> are the singular values of
  <math|\<Lambda\>> (TODO: prove this, and why is there a factor 2).
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
    <associate|Prop:IP|<tuple|2|1>>
    <associate|Prop:Proj|<tuple|5|2>>
    <associate|auto-1|<tuple|1|1>>
    <associate|auto-2|<tuple|2|2>>
    <associate|auto-3|<tuple|2.1|2>>
    <associate|auto-4|<tuple|3|3>>
    <associate|auto-5|<tuple|3.1|3>>
    <associate|auto-6|<tuple|3.2|3>>
    <associate|auto-7|<tuple|3.3|?>>
    <associate|auto-8|<tuple|3.4|?>>
    <associate|auto-9|<tuple|3.5|?>>
  </collection>
</references>

<\auxiliary>
  <\collection>
    <\associate|toc>
      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|1<space|2spc>Definitions>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-1><vspace|0.5fn>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|2<space|2spc>Incomplete
      information> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-2><vspace|0.5fn>

      <with|par-left|<quote|1tab>|2.1<space|2spc>Reconstructing missing
      information <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-3>>

      <vspace*|1fn><with|font-series|<quote|bold>|math-font-series|<quote|bold>|3<space|2spc>Approximate
      subrepresentations> <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-4><vspace|0.5fn>

      <with|par-left|<quote|1tab>|3.1<space|2spc>Computing the error
      affecting a subrepresentation <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-5>>

      <with|par-left|<quote|1tab>|3.2<space|2spc>Finding a pair
      <with|mode|<quote|math>|<around*|(|I,P|)>> close to an approximate pair
      <with|mode|<quote|math>|<around*|(|<wide|I|~>,<wide|P|~>|)>>
      <datoms|<macro|x|<repeat|<arg|x>|<with|font-series|medium|<with|font-size|1|<space|0.2fn>.<space|0.2fn>>>>>|<htab|5mm>>
      <no-break><pageref|auto-6>>
    </associate>
  </collection>
</auxiliary>
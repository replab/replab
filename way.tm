<TeXmacs|1.99.5>

<style|generic>

<\body>
  <doc-data|<doc-title|Hybrid techniques for representation decomposition>>

  Let <math|G> be a finite group, and <math|\<rho\>:G\<rightarrow\>GL<rsub|N><around*|(|<with|math-font|Bbb*|Q>|)>>
  a representation of <math|G>. Our aim is to obtain a decomposition of
  <math|\<rho\>> into irreducible representations either over
  <math|<with|math-font|Bbb*|Q>> or over <math|<with|math-font|Bbb*|R>>, with
  enough information such that the projection

  <\equation>
    \<Pi\><rsub|\<rho\>><around*|[|M|]>=<frac|1|<around*|\||G|\|>><big|sum><rsub|g\<in\>G>\<rho\><rsub|g>M\<rho\><rsub|g><rsup|\<dag\>>
  </equation>

  can be cheaply computed.

  We proceed as follows.

  Step 1.

  We compute the conjugacy classes of <math|G>.

  We compute the approximate change of basis matrix <math|<wide|U|~>> that
  reveals the isotypic components of <math|\<rho\>> over
  <math|<with|math-font|Bbb*|R>>, such that

  <\equation>
    <wide|U|~><rsup|\<dag\>>\<rho\><around*|(|g|)><wide|U|~>=<matrix|<tformat|<table|<row|<cell|<with|math-font|Bbb*|1>\<otimes\>\<tau\><rsub|1><around*|(|g|)>>|<cell|>|<cell|>>|<row|<cell|>|<cell|<with|math-font|Bbb*|1>\<otimes\>\<tau\><rsub|2><around*|(|g|)>>|<cell|>>|<row|<cell|>|<cell|>|<cell|\<ldots\>>>>>>
  </equation>

  For this change of basis and the conjugacy classes, we compute the
  approximate character table of the representations appearing in the
  decomposition of <math|\<rho\>>.

  We detect representations with rational characters by: testing whether\ 
</body>

<initial|<\collection>
</collection>>
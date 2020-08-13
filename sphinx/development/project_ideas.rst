Project ideas
=============

This lists ideas of self-contained projects to improve RepLAB; we put those here instead of clogging the GitHub issue list.


Faster computation of conjugacy classes
---------------------------------------

Currently, the search for conjugacy classes is pretty slow for groups of large order (we use random search).

This could be improved dramatically.

The most accessible method seems now the Butler1994 reference below. It needs the computation of Sylow subgroups, which could be implemented following Holt2005, chapter 4.

See references:

- G. Butler, “An Inductive Schema for Computing Conjugacy Classes in Permutation Groups,” Mathematics of Computation, vol. 62, no. 205, pp. 363–383, 1994, doi: 10.2307/2153415.

- A. Hulpke, “Computing conjugacy classes of elements in matrix groups,” Journal of Algebra, vol. 387, pp. 268–286, Aug. 2013, doi: 10.1016/j.jalgebra.2013.02.043.

- J. Cannon and B. Souvignier, “On the computation of conjugacy classes in permutation groups,” in Proceedings of the 1997 international symposium on Symbolic and algebraic computation, Kihei, Maui, Hawaii, USA, Jul. 1997, pp. 392–399, doi: 10.1145/258726.258855.

- D. F. Holt, B. Eick, and E. A. O’Brien, Handbook of Computational Group Theory. CRC Press, 2005.

Graph automorphisms
-------------------

RepLAB has a decent implementation of (non-partition-based) backtrack search. Thus, naive algorithms could be implemented, with a possible nauty/saucy/bliss backend.

- `<https://users.cecs.anu.edu.au/~bdm/nauty/>`_
- `<http://www.tcs.hut.fi/Software/bliss/>`_
- `<http://vlsicad.eecs.umich.edu/BK/SAUCY/>`_

Applications of (computational) representation theory
-----------------------------------------------------

- A. Zingoni, “Group-theoretic exploitations of symmetry in computational solid and structural mechanics,” International Journal for Numerical Methods in Engineering, vol. 79, no. 3, pp. 253–289, 2009, doi: 10.1002/nme.2576.

- A. Zingoni, “Group-theoretic applications in solid and structural mechanics: a review,” in Computational structures technology, GBR: Civil-Comp press, 2002, pp. 283–317.

Applied to differential equations: `<https://www.math.ubc.ca/~bluman/>`_

Support for efficient rational arithmetic
-----------------------------------------

Some computations always involve rational representations/bases, for example when decomposing permutation representations of (various products of) symmetric groups.

- We can first use an adhoc Matlab implementation that delegates to double-encoded integers or ``java.math.BigInteger``. See the ``replab.Rational`` class in this branch `<https://github.com/replab/replab/pull/356>`_.

- SuiteSparse has a rational linear equation system solver with Matlab support, `<https://github.com/DrTimothyAldenDavis/SuiteSparse/blob/master/SLIP_LU/MATLAB/SLIP_demo.m>`_

Implement matrices with elements in the cyclotomic field
--------------------------------------------------------

This would enable RepLAB to perform exact computations involving representations/characters.

One of our collaborators has an implementation lifted from GAP System, see `<https://github.com/denisrosset/cyclo>`_, see `<https://www.gap-system.org/Manuals/doc/ref/chap18.html>`_.

SAGE had an implementation, see `<https://trac.sagemath.org/ticket/8327>`_ and the top-level note on `<https://doc.sagemath.org/html/en/reference/number_fields/sage/rings/universal_cyclotomic_field.html>`_.

- T. Breuer, "Integral Bases for Subfields of Cyclotomic Fields", AAECC, vol. 8, no. 4, pp. 279–289, Apr. 1997, doi: 10.1007/s002000050065, `<https://link.springer.com/article/10.1007/s002000050065>`_.

Integrate support for manipulations of tableaux
-----------------------------------------------

Right now, we perform those computations in an adhoc way.

- Discuss with the author of `<https://www.mathworks.com/matlabcentral/fileexchange/62142-matrep-a-matlab-representation-theory-toolbox-symmetric-groups?focused=7410697&tab=function>`_ to integrate the functions in RepLAB?

Recognize direct products
-------------------------

- N. Kayal and T. Nezhmetdinov, "Factoring Groups Efficiently", in Automata, Languages and Programming, Berlin, Heidelberg, 2009, pp. 585–596, doi: 10.1007/978-3-642-02927-1_49.

Recognition of symmetric/alternating groups for large orders
------------------------------------------------------------

For larger orders, implement one of the randomized algorithms below:

- S. Bratus and I. Pak, “Fast Constructive Recognition of a Black Box Group Isomorphic to Sn or An using Goldbach’s Conjecture,” Journal of Symbolic Computation, vol. 29, no. 1, pp. 33–57, Jan. 2000, doi: 10.1006/jsco.1999.0295.

- R. Beals, C. R. Leedham-Green, A. C. Niemeyer, C. E. Praeger, and A. K. Seress, “A BLACK-BOX GROUP ALGORITHM FOR RECOGNIZING FINITE SYMMETRIC AND ALTERNATING GROUPS, I,” p. 17.

- S. Jambor, M. Leuner, A. C. Niemeyer, and W. Plesken, “Fast recognition of alternating groups of unknown degree,” Journal of Algebra, vol. 392, pp. 315–335, Oct. 2013, doi: 10.1016/j.jalgebra.2013.06.005.

Faster construction of irreducible representations of the symmetric group
-------------------------------------------------------------------------

See `<https://sage.math.leidenuniv.nl/src/combinat/yang_baxter_graph.py>`_

based on

- A. Lascoux, "Youngs's representations of the symmetric group," in Symmetry and Structural Properties of Condensed Matter, 0 vols., WORLD SCIENTIFIC, 2001, pp. 94–104.

Implement EPIMORPHISMS from Holt
--------------------------------

This enables the computation of isomorphisms between groups, and the automorphism group of a group.

Real representations
--------------------

Have better support/documentation for real representations

See `<https://www.maths.manchester.ac.uk/~jm/wiki/Representations/Representations>`

#import "@preview/rubber-article:0.5.2": *
#import "@preview/lovelace:0.3.1": *

// Layout and styling
#show: article.with(
  cols: none, // Tip: use #colbreak() instead of #pagebreak() to seamlessly toggle columns
  eq-chapterwise: true,
  eq-numbering: "(1.1)",
  header-display: true,
  header-title: "FEM Homework",
  lang: "en",
  page-margins: 1.75in,
  page-paper: "us-letter",
)

// Frontmatter
#maketitle(
  title: "FEM Homework",
  authors: ("Jules Herrmann",),
  date: datetime.today().display("[day]. [month repr:long] [year]"),
)

In this report, I explain my work for the homework of FEM class.

I chose to work on Project 4 : Elimination tree and direct sparse Cholesky solver.

I also implemented another sphere mesh, based on the geodesic subdivision of an icosahedron.



= Sparse Cholesky Decomposition

This work follows the method proposed in #cite(<scott2023algorithms>,form:"prose").

In the following I use $A^a$ for the "anti-transpose" of the matrix $A$. That is, the matrix obtained by mirroring along the anti-diagonal, or equivalently the matrice $tilde(I) A^T tilde(I)$ where $tilde(I)$ is a matrix with only ones on the antidiagonal.

== Symbolic phase

Function `void elimination_tree(const CSRMatrix &A,TArray<uint32_t> &parent)` is my
implementation of algorithm 4.2 of @scott2023algorithms.

The pseudo code uses an indexing scheme starting at $1$ and uses $0$ to denote the absence of a parent. Since I transcribe the algorithm to use $0$-based indexing, I chose to denote the absence of parent with the maximum integer `UINT32_MAX`.

#v(20pt)

Function `cholesky_sparsity_pattern` is my implementation of algorithme 4.3 of @scott2023algorithms.

If $A$ is the input sparse matrix, and $L L^T$ its Cholesky factorization, then this function returns both the sparsity pattern of $L$ and $L^a$. It procedes by creating the graph of $L$ in the adjacency list data structure (resp. `rowL` and `colL`). It fills it in according the algorithm, then build the sparsity patterns.


== Sparse triangular solver <sparsetrisolve>

Let's first descrive the data structure used to store a sparse vector.

Let $v$ be a vector with $"nnz"$ non-zero coefficients. This vector is decribed with a buffer `v` containg contiguously all the non-zero coefficients. A second buffer `v_ind` of type `uint32_t *` contains the index of the non-zero coefficients in such a way that if $v_i = alpha$, then there exist a $k$ such that $"v[k]" = alpha$ and $"v_ind[k]" = i$. Finally `v_size` stores $"nnz"$.

#v(20pt)

Two functions `sparse_tri_solve` have been implemented. One deals with dense vectors and a sparse matrix, while the other deals with sparse vectors.

The algorithm is a simple forward substitution algorithm, where, at step $"row"$, we compute

$
x_"row" = (b_"row" - sum_("col" = 1)^("row"-1) L_("row","col") x_"col")/L_("row","row")
$

The only difficulty is the computation of the sum in the case of sparse vectors. We want to iterate over $x$ and a row of $L$, which are both sparse vectors, so some care is required for the iteration to stay synchronized.

== Numerical phase

The numerical phase of the Cholesky factorization uses the _up-looking_ cholesky algorithm.

Indeed this algorithm can be implemented in such a way that the only iterations over consecutive values of the matrix happen along a row. Indeed, I use the Compressed Sparse Row (CSR) scheme to store my matrices and the scheme doesn't allow for fast iterations over a column.

#align(center,
pseudocode-list(booktabs: true,title:[Up-looking Cholesky factorization])[
  + $L_(11) <- (A_(11))^(1/2)$
  + *for* $i$ from $2$ to $N$
    + Solve $L_(1:i-1,1:i-1) z = A_(i,1:i-1)$
    + $L_(i,1:i-1) <- z$
    + $L_(i,i) <- (A_(i i) - z z^T)^(1/2)$
  + *end*
])

This method uses the triangular solver that I presented in @sparsetrisolve.

Since the rows of $A$, as well as their indices, are stored contiguously in memory, no data movement is required to call the sparse triangular solver at line 3. I just needed to compute the size of the sparse vector, that is, the number of non-zero coefficients of index inferior to $i$, which I do in linear time.

Function `sparse_cholesky` is my implementation of the algorithm.

== Sparse Cholesky solver

I created a `CholeskySolver` struct that contains the lower subdiagonal matrix $L$, as well as $L^a$.

$L^a$ is also lower subdiagonal, this way, it can be used with the forward substitution algithm I implemented.

$L^a$ can be used to solve linear problem with $L^T$ in the following way :

$
L^T x = b <=> tilde(I) L^T tilde(I) tilde(I) x = tilde(I) b <=> L^a x^r = b^r
$

where $v^r = tilde(I) v$ is the vector in reversed order.

Therefore the following algorithm is used to solve $A x = b$

#align(center,
pseudocode-list(booktabs: true,title:[Cholesky solver])[
  + Solve $L y = b$
  + Solve $L z = y^r$
  + return $z^r$
])

Method `CholeskySolver::solve` implements this algorithm.

A method `CholeskySolver::update_same_pattern` has also been implemented which allows to update the coefficients of the solver without recomputing the sparsity patterns. This is required because the interface of the demonstration program allows the user to change the value of $d t$ and $nu$.

= Better Sphere mesh

= Code Structure

= Build Instructions

#bibliography("refs.bib")

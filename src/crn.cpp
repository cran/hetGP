// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // Computes (Ut x Ux) diag(Lambda) (Ut x Ux) (for inverse computation, Utx kron(Ut, Ux)) 
// // [[Rcpp::export]]
// NumericMatrix KronInv(NumericMatrix Utx, NumericVector Lambda){
//   int n = Lambda.length();
//   NumericMatrix res(n, n);
//   
//   for(int i = 0; i < n; i++){
//     for(int j = i; j < n; j++){
//       for(int k = 0; k < n; k++)
//       res(i,j) = res(j, i) += Utx(i, k) * Lambda(k) * Utx(j, k);
//     }
//   }
//   
//   return res;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// /*** R
// library(hetGP)
// nx <- 30
// nt <- 20 
// x <- matrix(runif(nx), ncol = 1)
// t <- matrix(runif(nt), ncol = 1)
// Cx <- cov_gen(x, theta = 0.03)
// Ct <- cov_gen(t, theta = 0.01)
// g <- 1e-6
// C <- kronecker(Ct, Cx)
// SCx <- svd(Cx)
// SCt <- svd(Ct)
// 
// Ki <- tcrossprod(kronecker(SCt$u, SCx$u) * rep(1/(kronecker(SCt$d, SCx$d) + g), each = nx*nt), kronecker(SCt$u, SCx$u))
// Ki2 <- chol2inv(chol(C + diag(g, nx * nt)))
// Ki3 <- KronInv(kronecker(SCt$u, SCx$u), 1/(kronecker(SCt$d, SCx$d) + g))
// 
// range(Ki - Ki2)
// range(Ki - Ki3)
// 
// microbenchmark::microbenchmark(Ki <- tcrossprod(kronecker(SCt$u, SCx$u) * rep(1/(kronecker(SCt$d, SCx$d) + g), each = nx*nt), kronecker(SCt$u, SCx$u)),
//                                Ki2 <- chol2inv(chol(C + diag(g, nx * nt))),
//                                Ki3 <- KronInv(kronecker(SCt$u, SCx$u), 1/(kronecker(SCt$d, SCx$d) + g))
// )
// */

#include "cholesky.h"
#include <cmath>
#include <array.h>
#include <matrix.h>
#include <sparse_matrix.h>
#include <iostream>


int main(){
    CSRMatrix A;
    A.symmetric = true;
    A.cols = 2;
    A.rows = 2;
    A.nnz = 4;
    A.row_start = new uint32_t[3]{0, 2, 4};
    A.col = new uint32_t[4]{0, 1, 0, 1};
    A.data.resize(4);
    A.data[0] = 4.0;
    A.data[1] = 1.0;
    A.data[2] = 1.0;
    A.data[3] = 3.0;

    std::cout << "Matrix A: " << std::endl;
    for(size_t i = 0; i < A.rows; ++i){
        size_t start = A.row_start[i];
        size_t stop = A.row_start[i + 1];
        for(size_t k = start; k < stop; ++k){
            std::cout << "A(" << i << ", " << A.col[k] << ") = " << A.data[k] << std::endl;
        }
    }

    std::cout << "========================================" << std::endl;
    std::cout << "Test de Symbolic Cholesky " << std::endl;
    std::cout << "========================================" << std::endl;

    etree T;
    symbolic_cholesky(A, T);
    std::cout << "Elimination tree: " << std::endl;
    for(size_t i = 0; i < T.size; ++i){
        std::cout << "T[" << i << "] = " << T[i] << std::endl;
    }

    std::cout << "===============================" << std::endl;
    std::cout << "Test de L_Pattern: " << std::endl;
    std::cout << "===============================" << std::endl;
    //Resultat attendu : T[0] = 1, T[1] = 0
    CSRPattern L_Pattern;
    L_pattern(A, T, L_Pattern);
    std::cout << "Sparsity pattern of L: " << std::endl;
    for(size_t i = 0; i < L_Pattern.rows; ++i){
        size_t start = L_Pattern.row_start[i];
        size_t stop = L_Pattern.row_start[i + 1];
        for(size_t k = start; k < stop; ++k){
            std::cout << "L(" << i << ", " << L_Pattern.col[k] << ") is non zero" << std::endl;
        }
    }
    //Resultat attendu : L(0, 0) is non zero, L(1, 0) is non zero, L(1, 1) is non zero

    std::cout << "===============================" << std::endl;
    std::cout << "Testing Cholesky factorization" << std::endl;
    std::cout << "===============================" << std::endl;

    CSRMatrix L;

    cholesky_fact(A, L, L_Pattern);

    std::cout << "Résultat de la factorisation:" << std::endl;
    for(size_t i = 0; i < L.rows; ++i){
        size_t start = L.row_start[i];
        size_t stop = L.row_start[i + 1];
        for(size_t k = start; k < stop; ++k){
            std::cout << "L(" << i << ", " << L.col[k] << ") = " << L.data[k] << std::endl;
        }
    }

    std::cout<<"Matrice L attendue : " << std::endl;
    std::cout << "L(0, 0) = 2.0" << std::endl;
    std::cout << "L(1, 0) = 0.5" << std::endl;
    std::cout << "L(1, 1) = " <<sqrt(11)/2.0 << std::endl;
    
    std::cout << "========================================" << std::endl;
    std::cout << "Test des fonctions de résolution LL^Tx=b" << std::endl;
    std::cout << "========================================" << std::endl;

    TArray<double> b(2,1.0);
    TArray<double> x(2, 0.0);

    cholesky_solve(L,b,x);

    std::cout << "La solution est égale à:" << std::endl;
    for(size_t i = 0; i < x.size; ++i){
        std::cout << "x["<< i << "] = " << x[i] << std::endl;
    }

    std::cout << "Le résultat attendu est x[0] = " << (2.0/11.0) << " et x[1] = " << (3.0/11.0) << std::endl;
}
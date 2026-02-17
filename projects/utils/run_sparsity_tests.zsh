#!/usr/bin/env zsh

EXEC_SYMBOLIC=./build/release/test_symbolic
EXEC_CHOLESKY=./build/release/test_cholesky

OUT_SYMBOLIC=sparsity_results.csv
OUT_CHOLESKY=cholesky_timings.csv

START=2
END=20

echo "subdiv,n,P_nnz,L_nnz,L_D_nnz" > $OUT_SYMBOLIC
echo "subdiv,n,chol_time,chol_nested_time" > $OUT_CHOLESKY

for s in {$START..$END}; do
    echo "Running subdiv=$s"

    ############################################################
    # SYMBOLIC TEST
    ############################################################

    symbolic_output=$($EXEC_SYMBOLIC -s $s)

    subdiv=$(echo "$symbolic_output" | awk '/subdiv:/ {print $2}')
    n=$(echo "$symbolic_output" | awk '/n:/ {print $2}')
    P_nnz=$(echo "$symbolic_output" | awk '/P.nnz/ {print $3}')
    L_nnz=$(echo "$symbolic_output" | awk '/L.nnz/ {print $3}')
    L_D_nnz=$(echo "$symbolic_output" | awk '/L_D.nnz/ {print $3}')

    echo "$subdiv,$n,$P_nnz,$L_nnz,$L_D_nnz" >> $OUT_SYMBOLIC


    ############################################################
    # CHOLESKY TIMING TEST
    ############################################################

    chol_output=$($EXEC_CHOLESKY --subdiv $s)

subdiv_chol=$(echo "$chol_output" | awk '/subdiv:/ {print $2}')
n_chol=$(echo "$chol_output" | awk '/n:/ {print $2}')
chol_time=$(echo "$chol_output" | awk '/chol_time:/ {print $2}')
chol_nested_time=$(echo "$chol_output" | awk '/chol_nested_time:/ {print $2}')


    echo "$subdiv_chol,$n_chol,$chol_time,$chol_nested_time" >> $OUT_CHOLESKY

done

echo "Done."
echo "Symbolic results -> $OUT_SYMBOLIC"
echo "Cholesky timings -> $OUT_CHOLESKY"

python plots.py

import matplotlib.pyplot as plt
import pandas as pd

sparse = pd.read_csv('sparsity_results.csv')
times = pd.read_csv('cholesky_timings.csv')

fig, ax = plt.subplots()

plt.plot(sparse['n'], sparse[['P_nnz', 'L_nnz', 'L_D_nnz']], label = ['M', 'L', 'L - nested disscetion'])
plt.legend(title="Number of non-zeros")
plt.title('Number of non-zero entries as a function of matrix size')
plt.xlabel('Matrix size')
plt.ylabel('Non-zero entries')
plt.yscale('log')
plt.yscale('log')
plt.savefig('nonzero')

plt.figure()

plt.plot(times['n'], times[['chol_time', 'chol_nested_time']], label = ['Default ordering', 'Nested dissection'])
plt.legend()
plt.title('Computation time of solution using Cholesky decomposition')
plt.xlabel('Matrix size')
plt.ylabel('Time (s)')
plt.yscale('log')
plt.savefig('time')
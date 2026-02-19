import matplotlib.pyplot as plt
import pandas as pd

sparse = pd.read_csv('sparsity_results.csv')
times = pd.read_csv('cholesky_timings.csv')
times_ns = pd.read_csv('NS_timings.csv')

sparse['fill in L'] = sparse['L_nnz']/sparse['P_nnz']
sparse['fill in L_D'] = sparse['L_D_nnz']/sparse['P_nnz']
times_ns['n'] = times['n']


plt.plot(sparse['n'], sparse[['P_nnz', 'L_nnz', 'L_D_nnz']], label = ['M', 'L', 'L - nested disscetion'])
plt.legend(title="Number of non-zeros")
plt.title('Number of non-zero entries as a function number of degrees of freedom (Ndof)')
plt.xlabel('Ndof')
plt.ylabel('Non-zero entries')
plt.yscale('log')
plt.xscale('log')
plt.savefig('nonzero')

plt.figure()

plt.plot(sparse['n'], sparse[['fill in L', 'fill in L_D']], label = ['Default ordering', 'Nested dissection'])
plt.title('Fill-in as a function of a function number of degrees of freedom (Ndof)')
plt.legend()
plt.xlabel('Ndof')
plt.ylabel('Fill in')
plt.xscale('log')
plt.savefig('fill_in')

plt.figure()

plt.plot(times['n'], times[['chol_time_avg', 'chol_nested_time_avg']], label = ['Default ordering', 'Nested dissection'])
plt.legend()
plt.title('Computation time of Cholesky factor as a function of a function number of degrees of freedom (Ndof)')
plt.xlabel('Ndof')
plt.ylabel('Time (s)')
plt.yscale('log')
plt.xscale('log')
plt.savefig('time')

plt.figure()


plt.plot(times_ns['subdiv'], times_ns[['ns_time_avg', 'cholesky_time_avg']], label = ['CG', 'Cholesky'])
plt.title('Computation time of one time step in Navier Stokes Solver')
plt.xlabel('Ndof')
plt.ylabel('Time (s)')
plt.yscale('log')
plt.legend(title="Solution method")

plt.savefig('ns_time')
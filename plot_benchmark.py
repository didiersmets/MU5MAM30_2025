import matplotlib.pyplot as plt

# DoF values
dof = [2402, 5402, 9602, 15002, 21602]

# CG
cg = [
    1.083242495e-03,   
    3.673716966667e-03, 
    7.988859856667e-03,  
    1.591887660870e-02, 
    2.823359677258e-02 
]

# Cholesky on cube
chol_cube = [
    2.03600610e-03,      
    9.336062873333e-03,  
    2.716612513000e-02,  
    7.860209744000e-02,  
    6.140079360201e-02   
]

# Cholesky on cube2
chol_cube2 = [
    1.06682903e-03,      
    3.932715110000e-03, 
    9.560879926667e-03, 
    1.891312615667e-02,  
    3.743330821333e-02   
]

# CG on cube2
cg_cube2 = [
    1.016766866221e-03,  
    3.464933672241e-03,  
    9.269718566667e-03,  
    1.558798111000e-02,  
    2.921652596667e-02  
]

# Ignore missing values
def valid_xy(x, y):
    xx, yy = [], []
    for xi, yi in zip(x, y):
        if yi is not None:
            xx.append(xi)
            yy.append(yi)
    return xx, yy

# Plot
plt.figure(figsize=(8, 5))

x, y = valid_xy(dof, cg)
plt.plot(x, y, marker='o', label='CG (cube)')

x, y = valid_xy(dof, chol_cube)
plt.plot(x, y, marker='o', label='Cholesky (cube)')

x, y = valid_xy(dof, chol_cube2)
plt.plot(x, y, marker='o', label='Cholesky (cube2)')

x, y = valid_xy(dof, cg_cube2)
plt.plot(x, y, marker='o', linestyle='--', label='CG (cube2)')

plt.xlabel('DoF')
plt.ylabel('Average step time (s)')
plt.title('Navier-Stokes benchmark: CG vs Cholesky')
plt.grid(True)
plt.legend()
plt.tight_layout()

plt.yscale('log')
plt.xscale('log')


plt.savefig('benchmark_runtime_vs_dof.png', dpi=200)
plt.show()
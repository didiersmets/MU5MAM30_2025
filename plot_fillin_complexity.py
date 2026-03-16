import math
import matplotlib.pyplot as plt

N = [2402, 5402, 9602, 15002, 21602]

fill_cube = [
    244509,
    855169,
    2063429,
    4073289,
    7088749
]

fill_cube2 = [
    69318,
    174442,
    349628,
    585768,
    861292
]

# Keep only non-missing values
def valid_xy(x, y):
    xx, yy = [], []
    for xi, yi in zip(x, y):
        if yi is not None:
            xx.append(xi)
            yy.append(yi)
    return xx, yy

# Reference curve C * N * log(N)
x2, y2 = valid_xy(N, fill_cube2)

if len(x2) == 0:
    raise ValueError("No valid cube2 fill-in data.")

C = y2[0] / (x2[0] * math.log(x2[0]))
ref_nlogn = [C * x * math.log(x) for x in x2]

# N^(4/3)
C43 = y2[0] / (x2[0] ** (4.0 / 3.0))
ref_n43 = [C43 * (x ** (4.0 / 3.0)) for x in x2]

# Plot
plt.figure(figsize=(8, 5))

x1, y1 = valid_xy(N, fill_cube)
if len(x1) > 0:
    plt.plot(x1, y1, marker='o', label='Measured fill-in (cube)')

plt.plot(x2, y2, marker='o', label='Measured fill-in (cube2)')
plt.plot(x2, ref_nlogn, '--', label=r'Reference $C\,N\log N$')
plt.plot(x2, ref_n43, ':', label=r'Reference $C\,N^{4/3}$')

plt.xscale('log')
plt.yscale('log')

plt.xlabel('DoF $N$')
plt.ylabel(r'Fill-in = $nnz(L)-nnz(A)$')
plt.title('Measured fill-in vs theoretical growth')
plt.grid(True, which='both', linestyle='--', linewidth=0.5)
plt.legend()
plt.tight_layout()

plt.savefig('fillin_complexity.png', dpi=200)
plt.show()

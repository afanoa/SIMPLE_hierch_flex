import sys, numpy as np
# compare_latent.py A B NODE: per-dimension |correlation| and relative RMS difference of two arms' focused-flex latent at NODE (rows matched by rows.txt)
A, B, node = sys.argv[1:4]
def load(run):
    rows = np.loadtxt(f"{run}/1_hetero_tree/{node}/rows.txt", dtype=int)
    x = np.loadtxt(f"{run}/1_hetero_tree/{node}/fflex/flex_pca_coordinates.txt")
    if x.ndim == 1: x = x[:, None]
    if x.shape[0] != rows.size and x.shape[0] > rows.size: x = x[rows - 1] if x.shape[0] >= rows.max() else x
    return rows, x
ra, xa = load(A); rb, xb = load(B)
common, ia, ib = np.intersect1d(ra, rb, return_indices=True)
xa, xb = xa[ia], xb[ib]
print(f"{node}: {common.size} common particles of {ra.size} / {rb.size}; dims {xa.shape[1]} / {xb.shape[1]}")
for d in range(min(xa.shape[1], xb.shape[1])):
    a, b = xa[:, d], xb[:, d]; c = np.corrcoef(a, b)[0, 1]
    s = np.sign(c) if c != 0 else 1.
    rel = np.sqrt(np.mean((a - s * b) ** 2)) / np.sqrt(np.mean(a ** 2))
    print(f"  dim {d+1}: |corr| {abs(c):.4f}  rel rms diff {rel:.4f}  sd {a.std():.4g} / {b.std():.4g}")

import matplotlib.pyplot as plt
from tqdm import tqdm

outfile = "atoms.png"

def plot_grid(filename='release/grid.dat'):
  with open(filename) as grid_file:
    fig, ax = plt.subplots()
    fig.set_size_inches(10, 4)
    ax.set_xlim((-1, 500))
    ax.set_xticks([tick for tick in range(0, 501, 100)])
    max_j = -1
    for line in tqdm(grid_file):
      if line.startswith('#'):
        continue
      i, j, atom_type, *_ = line.split()
      i, j, atom_type = int(i), int(j), int(atom_type)
      if j > max_j and atom_type in (1, 2):
        max_j = j
      # ax.scatter(i, j, color='blue' if atom_type == 1 else 'green', s=1)
      ax.add_patch(plt.Circle((i, j), 0.5, color='blue' if atom_type == 1 else 'green'))
    ax.set_ylim((31, max_j+1))
    ax.set_yticks([tick for tick in range(31, max_j+1)])
    print("Saving...")
    plt.savefig(outfile)
    print(f"Plot saved at {outfile}")

plot_grid()

import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import argparse
import logging
import os
from pathlib import Path

time_format = '%H:%M:%S'
log_format = '[%(levelname)s] %(asctime)s.%(msecs)03d [%(funcName)s]: %(message)s'
logging.basicConfig(format=log_format, datefmt=time_format)
logger = logging.getLogger("main")
logger.setLevel(logging.INFO)

def plot_grid(args):
  infile_name = "grid.dat"
  outfile_name = "atoms.png"

  infile_path = args.input_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep'), sep=r'\s+')
  logger.info("Reading done")

  logger.info(f"Generating grid plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.empty((img_height, img_width))
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])
    grid[y_idx, x_idx] = row['type']
  
  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  cmap = mpl.colors.ListedColormap(('white', 'blue', 'green'))
  norm = mpl.colors.BoundaryNorm((0, 1, 2, 3), cmap.N)
  ax.set_xticks([1] + list(range(100, x_max+2, 100)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  ax.imshow(grid, origin='lower', cmap=cmap, norm=norm, interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  logger.info(f"Generation done")

  outfile_path = args.output_dir / outfile_name
  logger.info(f"Saving grid plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info(f"Plot saved")
  
  if not args.quiet:
    plt.show()
  plt.close()

def plot_energies(args):
  infile_name = "grid.dat"
  outfile_name = "energies.png"
  
  infile_path = args.input_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep'), sep=r'\s+')
  logger.info("Reading done")
  
  logger.info(f"Generating energy plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.empty((img_height, img_width))
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])
    grid[y_idx, x_idx] = row['ep']
  
  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  ax.set_xticks([1] + list(range(100, x_max+2, 100)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  im = ax.imshow(grid, origin='lower', cmap="Reds", interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  cax = fig.add_axes((ax.get_position().x1+0.01, ax.get_position().y0,0.02, ax.get_position().height))
  plt.colorbar(im, cax=cax, ticks=np.floor(np.linspace(df['ep'].min(), df['ep'].max(), 3)*1000)/1000, label='Energy [$eV$]')
  logger.info(f"Generation done")

  outfile_path = args.output_dir / outfile_name
  logger.info(f"Saving energy plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info("Plot saved")

  if not args.quiet:
    plt.show()
  plt.close()

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("-i", "--input_dir", type=Path, required=False, default="release")
  parser.add_argument("-o", "--output_dir", type=Path, required=False, default=".")
  parser.add_argument("-q", "--quiet", action="store_true", help="Pass this flag to prevent program from showing interactive plots")
  return parser.parse_args()

if __name__ == "__main__":
  args = parse_args()
  if not (os.path.exists(args.output_dir) and os.path.isdir(args.output_dir)):
    os.mkdir(args.output_dir)
  plot_grid(args)
  plot_energies(args)

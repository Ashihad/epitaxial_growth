import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

import argparse
import logging
import os
import sys
from pathlib import Path

time_format = '%H:%M:%S'
log_format = '[%(levelname)s] %(asctime)s.%(msecs)03d [%(funcName)s]: %(message)s'
logging.basicConfig(format=log_format, datefmt=time_format)
logger = logging.getLogger("main")
logger.setLevel(logging.INFO)

def plot_grid(args):
  infile_name = "grid.dat"
  outfile_name = "atoms.png"

  infile_path = args.work_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep'), sep=r'\s+')
  logger.info("Reading done")

  logger.info(f"Generating grid plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.zeros((img_height, img_width))
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])
    if y_idx != 0 and y_idx != y_max:
      grid[y_idx-1, x_idx-1] = row['type']
  
  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  cmap = mpl.colors.ListedColormap(('white', 'navy', 'goldenrod'))
  norm = mpl.colors.BoundaryNorm((0, 1, 2, 3), cmap.N)
  ax.set_xticks([1] + list(range(25, x_max+2, 25)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  ax.imshow(grid, origin='lower', cmap=cmap, norm=norm, interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  logger.info(f"Generation done")

  outfile_path = args.work_dir / outfile_name
  logger.info(f"Saving grid plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info(f"Plot saved")
  
  if not args.quiet:
    plt.show()
  plt.close()

def plot_energies(args):
  infile_name = "grid.dat"
  outfile_name = "energies.png"
  
  infile_path = args.work_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep'), sep=r'\s+')
  logger.info("Reading done")
  
  logger.info(f"Generating energy plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.empty((img_height, img_width))
  grid.fill(np.nan)
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])
    if y_idx != 0 and y_idx != y_max:
      grid[y_idx-1, x_idx-1] = row['ep']
  

  cmap = plt.get_cmap('gnuplot').copy()
  cmap.set_bad(color='white')
  grid_masked = np.ma.masked_where(grid == np.nan, grid)

  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  ax.set_xticks([1] + list(range(25, x_max+2, 25)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  im = ax.imshow(grid_masked, origin='lower', cmap=cmap, interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  cax = fig.add_axes((ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height))
  plt.colorbar(im, cax=cax, ticks=np.linspace(df['ep'].min(), df['ep'].max(), 5), label='Energy [$eV$]')
  logger.info(f"Generation done")

  outfile_path = args.work_dir / outfile_name
  logger.info(f"Saving energy plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info("Plot saved")

  if not args.quiet:
    plt.show()
  plt.close()

def plot_energies_contour(args):
  infile_name = "grid.dat"
  outfile_name = "energies_contour.png"
  
  infile_path = args.work_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep'), sep=r'\s+')
  logger.info("Reading done")
  
  logger.info(f"Generating energy contour plot...")

  # filter out lowest row
  df = df[df['y'] > 0]

  x_unique = np.sort(df['x'].unique())
  y_unique = np.sort(df['y'].unique())
  X, Y = np.meshgrid(x_unique, y_unique)
  Z = df.pivot(index='y', columns='x', values='ep').values
  x_max = df['x'].max()
  y_max = df['y'].max()

  # fill the missing row twice
  row_of_nans = np.empty_like(Z[-1])
  row_of_nans.fill(np.nan)
  Z = np.append(Z, [row_of_nans], axis=0)
  Z = np.append(Z, [row_of_nans], axis=0)

  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  ax.set_xticks([1] + list(range(25, x_max+2, 25)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  levels = 10

  im = ax.contourf(Z, levels=levels, origin='lower', cmap="gnuplot")#, extent=(0, x_max+0, 0, y_max+0))
  # ax.contour(Z, levels=levels, colors='black', linewidths=0.5, extent=(0, x_max+0, 0, y_max+0))
  cax = fig.add_axes((ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height))
  plt.colorbar(im, cax=cax, label='Energy [$eV$]')
  logger.info(f"Generation done")

  outfile_path = args.work_dir / outfile_name
  logger.info(f"Saving energy plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info("Plot saved")

  if not args.quiet:
    plt.show()
  plt.close()

def plot_dE_duv(args):
  infile_name = "energies.dat"
  outfile_name = "dEduv.png"
  
  infile_path = args.work_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'ep', 'dEduv'), sep=r'\s+')
  # filter out vacuum
  df = df[df['type'] != 0]
  logger.info("Reading done")
  
  logger.info(f"Generating energy derivative plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.empty((img_height, img_width))
  grid.fill(np.nan)
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])

    if y_idx != 0 and y_idx != y_max:
      grid[y_idx-1, x_idx-1] = row['dEduv']

  cmap = plt.get_cmap('plasma').copy()
  cmap.set_bad(color='white')
  grid_masked = np.ma.masked_where(grid == np.nan, grid)

  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  ax.set_xticks([1] + list(range(25, x_max+2, 25)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  im = ax.imshow(grid_masked, origin='lower', cmap=cmap, interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  cax = fig.add_axes((ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height))
  cb = fig.colorbar(im, cax=cax)
  cb.set_label(r'$\frac{dE}{duv}$', rotation=0, labelpad=15, fontsize=15)
  logger.info(f"Generation done")

  outfile_path = args.work_dir / outfile_name
  logger.info(f"Saving energy plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info("Plot saved")

  if not args.quiet:
    plt.show()
  plt.close()

def plot_displacements(args):
  infile_name = "tmp.dat"
  outfile_name = "displacements.png"
  
  infile_path = args.work_dir / infile_name
  logger.info(f"Reading from {infile_path}...")
  df = pd.read_csv(infile_path, header=1, names=('x', 'y', 'type', 'u', 'v', 'boundary1', 'boundary2', 'gradx', 'grady', 'ep'), sep=r'\s+')
  logger.info("Reading done")

  # filter out vacuum
  df = df[df['type'] != 0]
  
  logger.info(f"Generating displacements plot...")
  x_min, x_max = df['x'].min(), df['x'].max()
  y_min, y_max = df['y'].min(), df['y'].max()
  img_height = y_max - y_min + 1
  img_width = x_max - x_min + 1

  grid = np.empty((img_height, img_width))
  grid.fill(np.nan)
  for _, row in df.iterrows():
    x_idx = int(row['x'])
    y_idx = int(row['y'])
    if y_idx != 0 and y_idx != y_max:
      grid[y_idx-1, x_idx-1] = np.sqrt(row['u']**2 + row['v']**2)
  

  cmap = plt.get_cmap('cividis').copy()
  cmap.set_bad(color='white')
  grid_masked = np.ma.masked_where(grid == np.nan, grid)

  fig, ax = plt.subplots()
  fig.set_size_inches(15, 6)
  ax.set_xticks([1] + list(range(25, x_max+2, 25)))
  ax.set_yticks([1] + list(range(10, y_max+2, 10)))
  im = ax.imshow(grid_masked, origin='lower', cmap=cmap, interpolation=None, extent=(1, x_max+1, 1, y_max+1))
  cax = fig.add_axes((ax.get_position().x1+0.01, ax.get_position().y0, 0.02, ax.get_position().height))
  plt.colorbar(im, cax=cax, label=r'$\vec{U}$')
  logger.info(f"Generation done")

  outfile_path = args.work_dir / outfile_name
  logger.info(f"Saving energy plot at {outfile_path}...")
  plt.savefig(outfile_path)
  logger.info("Plot saved")

  if not args.quiet:
    plt.show()
  plt.close()

def parse_args():
  parser = argparse.ArgumentParser()
  parser.add_argument("-w", "--work_dir", type=Path, required=False, default="release")
  parser.add_argument("-q", "--quiet", action="store_true", help="Pass this flag to prevent program from showing interactive plots")
  return parser.parse_args()

if __name__ == "__main__":
  args = parse_args()
  if not (os.path.exists(args.work_dir) and os.path.isdir(args.work_dir)):
    logger.critical(f"Invalid \"-w\" argument (argument is not a directory, got \"{args.work_dir}\")")
    sys.exit(1)
  # plot_grid(args)
  # plot_energies(args)
  # plot_energies_contour(args)
  plot_dE_duv(args)
  # plot_displacements(args)

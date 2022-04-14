#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
for {frame number} vs {average intensity} data.
"""


import os

import matplotlib.pyplot as plt
import numpy as np


DATA_DIRECTORY = 'data'


def load_data():
  data_from_file_name = {}
  for path, _, file_names in os.walk(DATA_DIRECTORY):
    for file_name in sorted(file_names):
      full_file_name = os.path.join(path, file_name)
      data_from_file_name[file_name] = np.loadtxt(full_file_name, skiprows=1)
  return data_from_file_name


def make_plot(file_name, data):
  
  frame_data = data[:, 0]
  intensity_data = data[:, 1]
  
  figure, axes = plt.subplots()
  axes.plot(frame_data, intensity_data)
  axes.set(
    title=file_name,
    xlabel='Frame number',
    ylabel='Average intensity',
  )
  plt.show()


def main():
  data_from_file_name = load_data()
  for file_name, data in data_from_file_name.items():
    make_plot(file_name, data)


if __name__ == '__main__':
  main()

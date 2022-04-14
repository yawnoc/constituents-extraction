#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
for {frame number} vs {average intensity} data.
"""


import os
import numpy as np


DATA_DIRECTORY = 'data'


def load_data():
  data_from_file_name = {}
  for path, _, file_names in os.walk(DATA_DIRECTORY):
    for file_name in file_names:
      full_file_name = os.path.join(path, file_name)
      data_from_file_name[file_name] = np.loadtxt(full_file_name, skiprows=1)
  return data_from_file_name


def main():
  data_from_file_name = load_data()
  print(data_from_file_name)


if __name__ == '__main__':
  main()

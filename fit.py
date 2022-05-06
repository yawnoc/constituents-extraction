#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
from frame vs intensity data in `data/`.
"""


import os

import numpy as np


DATA_DIRECTORY = 'data'
DATA_COMMENT_LINE_COUNT = 1


def load_data_points():
  """
  Load (frame, intensity) data points from files.
  """
  
  data_points_from_file_name = {}
  
  for path, _, file_names in os.walk(DATA_DIRECTORY):
    for file_name in sorted(file_names):
      full_file_name = os.path.join(path, file_name)
      data_points_from_file_name[file_name] = \
              np.loadtxt(full_file_name, skiprows=DATA_COMMENT_LINE_COUNT)
  
  return data_points_from_file_name


def normalise(data):
  """
  Normalise data (for a variable) to between 0 and 1.
  """
  
  data_min = data.min()
  data_max = data.max()
  normalised_data = (data - data_min) / (data_max - data_min)
  
  return data_min, data_max, normalised_data


def main():
  
  data_points_from_file_name = load_data_points()
  for file_name, data_points in data_points_from_file_name.items():
    
    frame_data = data_points[:, 0]
    intensity_data = data_points[:, 1]
    
    frame_min, frame_max, x_data = normalise(frame_data)
    intensity_min, intensity_max, y_data = normalise(intensity_data)


if __name__ == '__main__':
  main()

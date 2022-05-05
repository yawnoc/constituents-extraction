#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
for {frame number} vs {average intensity} data.
"""


import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as st


DATA_DIRECTORY = 'data'
OUTPUT_DIRECTORY = 'output'


def load_data():
  data_from_file_name = {}
  for path, _, file_names in os.walk(DATA_DIRECTORY):
    for file_name in sorted(file_names):
      full_file_name = os.path.join(path, file_name)
      data_from_file_name[file_name] = np.loadtxt(full_file_name, skiprows=1)
  return data_from_file_name


def background_model(x, m, b):
  """
  Straight line for normalised background intensity.
  """
  
  return m * x + b


def peak_model(x, k, mu, sigma):
  """
  Exponentially modified Gaussian for each peak of normalised intensity.
  
  Given by
    f(x, K) = 1/(2K) exp(1/(2K^2) - x/K) erfc(-(x - 1/K) / sqrt(2)),
  transforms to the Wikipedia version according to K = 1/(sigma lambda).
  See <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html>
  and <https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution>.
  """
  return st.exponnorm.pdf(x, k, mu, sigma)


def normalise(frame_data, intensity_data):
  """
  Normalise frame and intensity data to between 0 and 1.
  
  For brevity, the normalised variables are simply called x and y.
  """
  
  frame_min = frame_data.min()
  frame_max = frame_data.max()
  intensity_min = intensity_data.min()
  intensity_max = intensity_data.max()
  
  x_data = (frame_data - frame_min) / (frame_max - frame_min)
  y_data = (intensity_data - intensity_min) / (intensity_max - intensity_min)
  
  return frame_min, frame_max, intensity_min, intensity_max, x_data, y_data


def make_plot(file_name, data):
  
  frame_data = data[:, 0]
  intensity_data = data[:, 1]
  
  frame_min, frame_max, intensity_min, intensity_max, x_data, y_data = \
          normalise(frame_data, intensity_data)
  
  figure, axes = plt.subplots()
  axes.plot(frame_data, intensity_data)
  axes.set(
    title=file_name,
    xlabel='Frame number',
    ylabel='Average intensity',
  )
  plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'{file_name}.pdf'))
  plt.show()


def main():
  data_from_file_name = load_data()
  for file_name, data in data_from_file_name.items():
    make_plot(file_name, data)


if __name__ == '__main__':
  main()

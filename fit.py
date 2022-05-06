#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
from frame vs intensity data in `data/`.
"""


import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.special as sp


DATA_DIRECTORY = 'data'
DATA_COMMENT_LINE_COUNT = 1
OUTPUT_DIRECTORY = 'output'


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


def background_function(x, m, b):
  """
  Straight line for the background intensity.
  """
  return m * x + b


def peak_function(x, h, mu, sigma, tau):
  """
  Exponentially modified Gaussian for each peak in the intensity.
  
  We use the chromatography version, see Kalambet et al. (2011),
  "Reconstruction of chromatographic peaks using the exponentially modified
  Gaussian function", Journal of Chemometrics, 25(7), 352-356,
  <https://doi.org/10.1002/cem.1343>:
          f(x) =
                  h sigma/tau sqrt(pi)/2
                  . exp[1/2 (sigma/tau)^2 - (x - mu)/tau]
                  . erfc[1/sqrt(2) (sigma/tau - (x - mu)/sigma)]
  where
          h is Gaussian amplitude,
          mu is Gaussian mean,
          sigma is Gaussian standard deviation,
          tau is exponential relaxation time.
  """
  return (
    h * sigma/tau * np.sqrt(np.pi)/2
      * np.exp(1/2 * (sigma/tau)**2 - (x - mu)/tau)
      * sp.erfc(1/np.sqrt(2) * (sigma/tau - (x - mu)/sigma))
  )


def heuristic_background_parameter_guesses(x_data, y_data):
  """
  Heuristically obtain guesses for the background function parameters.
  """
  
  window_size = 4
  x1 = np.mean(x_data[:window_size])
  y1 = np.mean(y_data[:window_size])
  x2 = np.mean(x_data[-window_size:])
  y2 = np.mean(y_data[-window_size:])
  
  m_guess = (y2 - y1) / (x2 - x1)
  b_guess = y1 - m_guess * x1
  
  return m_guess, b_guess


def main():
  
  data_points_from_file_name = load_data_points()
  for file_name, data_points in data_points_from_file_name.items():
    
    frame_data = data_points[:, 0]
    intensity_data = data_points[:, 1]
    
    frame_min, frame_max, x_data = normalise(frame_data)
    intensity_min, intensity_max, y_data = normalise(intensity_data)
    
    m_guess, b_guess = heuristic_background_parameter_guesses(x_data, y_data)
    
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(
      x_data,
      background_function(x_data, m_guess, b_guess),
      label='background guess',
      linestyle='dotted',
    )
    axes.set(
      title=file_name,
      xlabel=f'Normalised frame number [{int(frame_min)}, {int(frame_max)}]',
      ylabel=f'Normalised intensity [{intensity_min:.3}, {intensity_max:.3}]',
    )
    axes.legend()
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'{file_name}.pdf'))


if __name__ == '__main__':
  main()

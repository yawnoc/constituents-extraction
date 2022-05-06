#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
from frame vs intensity data in `data/`.
"""


import collections
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
  Heuristically make guesses for the background function parameters.
  
  1. Determines the centre of mass for a cluster at the start of the data
     and a cluster at the end of the data.
  2. Fits a straight line between these two points.
  
  Dumb, but works.
  """
  
  cluster_size = 10
  x1 = np.mean(x_data[:cluster_size])
  y1 = np.mean(y_data[:cluster_size])
  x2 = np.mean(x_data[-cluster_size:])
  y2 = np.mean(y_data[-cluster_size:])
  
  m_guess = (y2 - y1) / (x2 - x1)
  b_guess = y1 - m_guess * x1
  
  return m_guess, b_guess


def heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess):
  """
  Heuristically find peak and valley locations.
  
  1-peak case: find (foot1, peak, foot2).
  2-peak case: find (foot1, peak1, valley, peak2, foot2).
  
  1. Finds the two feet.
  2. Climbs upward from the two feet,
     checking if the secant line hits a valley.
  
  WARNING: assumes the data is such that there is either 1 peak or 2 peaks.
  No sanity checking is implemented.
  """
  
  y_data_foreground = y_data - background_function(x_data, m_guess, b_guess)
  
  x_climb = collections.deque(x_data)
  y_climb_foreground = collections.deque(y_data_foreground)
  
  foot_y_threshold = 0.05
  while y_climb_foreground[0] < foot_y_threshold:
    x_climb.popleft()
    y_climb_foreground.popleft()
    x_foot1 = x_climb[0]
    y_foot1_foreground = y_climb_foreground[0]
  while y_climb_foreground[-1] < foot_y_threshold:
    x_climb.pop()
    y_climb_foreground.pop()
    x_foot2 = x_climb[-1]
    y_foot2_foreground = y_climb_foreground[-1]
  y_foot1 = y_foot1_foreground + background_function(x_foot1, m_guess, b_guess)
  y_foot2 = y_foot2_foreground + background_function(x_foot2, m_guess, b_guess)
  
  distinct_x_threshold = 0.08
  x_valley = y_valley_foreground = None
  while len(x_climb) > 2:
    x1 = x_climb[0]
    y1_foreground = y_climb_foreground[0]
    x2 = x_climb[-1]
    y2_foreground = y_climb_foreground[-1]
    valley_candidates = \
            [
              (x, y_foreground)
                for x, y_foreground in zip(x_climb, y_climb_foreground)
                if x1 + distinct_x_threshold < x < x2 - distinct_x_threshold
                if (
                  y1_foreground
                    +
                  (y2_foreground - y1_foreground)/(x2 - x1) * (x - x1)
                          >
                  y_foreground
                )
            ]
    try:
      x_valley, y_valley_foreground = \
              min(valley_candidates, key=lambda pair: pair[1])
      break
    except ValueError:
      x_climb.popleft()
      y_climb_foreground.popleft()
      x_climb.pop()
      y_climb_foreground.pop()
  
  try:
    y_valley = \
            (
              y_valley_foreground
              + background_function(x_valley, m_guess, b_guess)
            )
  except TypeError:
    y_valley = None
  return x_foot1, y_foot1, x_valley, y_valley, x_foot2, y_foot2


def main():
  
  data_points_from_file_name = load_data_points()
  for file_name, data_points in data_points_from_file_name.items():
    
    frame_data = data_points[:, 0]
    intensity_data = data_points[:, 1]
    
    frame_min, frame_max, x_data = normalise(frame_data)
    intensity_min, intensity_max, y_data = normalise(intensity_data)
    
    m_guess, b_guess = heuristic_background_parameter_guesses(x_data, y_data)
    x_foot1, y_foot1, x_valley, y_valley, x_foot2, y_foot2 = \
            heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess)
    
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(
      x_data,
      background_function(x_data, m_guess, b_guess),
      label='background guess',
      linestyle='dotted',
    )
    axes.plot(
      [x_foot1, x_foot2, x_valley],
      [y_foot1, y_foot2, y_valley],
      'rx',
    )
    axes.set(
      title=file_name,
      xlabel=f'Normalised frame number [{int(frame_min)}, {int(frame_max)}]',
      ylabel=f'Normalised intensity [{intensity_min:.3}, {intensity_max:.3}]',
    )
    axes.legend()
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'guesses-{file_name}.pdf'))


if __name__ == '__main__':
  main()

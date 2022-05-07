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
  Heuristically make guesses for background function parameters.
  
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


def extract_y_coordinate(xy):
  return xy[1]


def heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess):
  """
  Heuristically find peak and valley locations.
  
  - 1 peak locations: (foot1, peak, foot2).
  - 2 peak locations: (foot1, peak1, valley, peak2, foot2).
  
  1. Finds the two feet.
  2. Climbs upward from the two feet,
     checking if the secant line hits a valley.
  
  WARNING:
  - Not robust; got lucky with `foot_y_threshold` and `distinct_x_threshold`.
  - No sanity checking is implemented;
    assumes the data is such that there be exactly 1 or 2 peaks.
  """
  
  y_data_foreground = y_data - background_function(x_data, m_guess, b_guess)
  
  x_climb = collections.deque(x_data)
  y_climb_foreground = collections.deque(y_data_foreground)
  
  foot_y_threshold = 0.05
  x_foot1 = y_foot1_foreground = None
  x_foot2 = y_foot2_foreground = None
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
              min(valley_candidates, key=extract_y_coordinate)
      break
    except ValueError:
      x_climb.popleft()
      y_climb_foreground.popleft()
      x_climb.pop()
      y_climb_foreground.pop()
  
  try: # 2 peaks
    
    y_valley = \
            (
              y_valley_foreground
              + background_function(x_valley, m_guess, b_guess)
            )
    x_peak1, y_peak1 = \
            max(
              [
                (x, y)
                  for x, y in zip(x_data, y_data)
                  if x_foot1 < x < x_valley
              ],
              key=extract_y_coordinate,
            )
    x_peak2, y_peak2 = \
            max(
              [
                (x, y)
                  for x, y in zip(x_data, y_data)
                  if x_valley < x < x_foot2
              ],
              key=extract_y_coordinate,
            )
    
    return (
      x_foot1, y_foot1,
      x_peak1, y_peak1,
      x_valley, y_valley,
      x_peak2, y_peak2,
      x_foot2, y_foot2,
    )
  
  except TypeError: # 1 peak
    x_peak, y_peak = max(zip(x_data, y_data), key=extract_y_coordinate)
    return x_foot1, y_foot1, x_peak, y_peak, x_foot2, y_foot2


def heuristic_peak_parameter_guesses(x_foot1, x_peak, y_peak, x_foot2):
  """
  Heuristically make guesses for peak function parameters.
  
  No mathematical theory underlies these guesses,
  which are based purely on intuition.
  """
  
  h_guess = 1.5 * y_peak * np.sqrt((x_foot2 - x_peak) / (x_peak - x_foot1))
  mu_guess = x_peak - 0.2 * (x_peak - x_foot1)
  sigma_guess = 0.4 * (x_peak - x_foot1)
  tau_guess = 0.33 * (x_foot2 - x_peak)
  
  return h_guess, mu_guess, sigma_guess, tau_guess


def main():
  
  data_points_from_file_name = load_data_points()
  for file_name, data_points in data_points_from_file_name.items():
    
    frame_data = data_points[:, 0]
    intensity_data = data_points[:, 1]
    
    frame_min, frame_max, x_data = normalise(frame_data)
    intensity_min, intensity_max, y_data = normalise(intensity_data)
    
    ################################
    # Background guess
    ################################
    m_guess, b_guess = heuristic_background_parameter_guesses(x_data, y_data)
    y_background_fit_guess = background_function(x_data, m_guess, b_guess)
    
    try: # 2 peaks
      
      (
        x_foot1, y_foot1,
        x_peak1, y_peak1,
        x_valley, y_valley,
        x_peak2, y_peak2,
        x_foot2, y_foot2,
      ) = \
              heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess)
      peak_valley_x_locations = \
              [x_foot1, x_peak1, x_valley, x_peak2, x_foot2]
      peak_valley_y_locations = \
              [y_foot1, y_peak1, y_valley, y_peak2, y_foot2]
      
      h_guess = mu_guess = sigma_guess = tau_guess = None
      y_peak_fit_guess_with_background = None
      
      ################################
      # Peak 1 guess
      ################################
      x_peak1_foot1 = x_foot1
      x_peak1_foot2 = \
              x_peak1 + (x_valley - x_peak1) / (y_peak1 - y_valley) * y_peak1
      h1_guess, mu1_guess, sigma1_guess, tau1_guess = \
              heuristic_peak_parameter_guesses(
                x_peak1_foot1,
                x_peak1, y_peak1,
                x_peak1_foot2
              )
      y_peak1_fit_guess = \
              peak_function(
                x_data,
                h1_guess, mu1_guess, sigma1_guess, tau1_guess,
              )
      y_peak1_fit_guess_with_background = \
              y_peak1_fit_guess + y_background_fit_guess
      
      ################################
      # Peak 2 guess
      ################################
      x_peak2_foot1 = \
              x_peak2 - (x_peak2 - x_valley) / (y_peak2 - y_valley) * y_peak2
      x_peak2_foot2 = x_foot2
      h2_guess, mu2_guess, sigma2_guess, tau2_guess = \
              heuristic_peak_parameter_guesses(
                x_peak2_foot1,
                x_peak2, y_peak2,
                x_peak2_foot2
              )
      y_peak2_fit_guess = \
              peak_function(
                x_data,
                h2_guess, mu2_guess, sigma2_guess, tau2_guess,
              )
      y_peak2_fit_guess_with_background = \
              y_peak2_fit_guess + y_background_fit_guess
    
    except ValueError: # 1 peak
      
      x_foot1, y_foot1, x_peak, y_peak, x_foot2, y_foot2 = \
              heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess)
      peak_valley_x_locations = [x_foot1, x_peak, x_foot2]
      peak_valley_y_locations = [y_foot1, y_peak, y_foot2]
      
      ################################
      # Peak guess
      ################################
      h_guess, mu_guess, sigma_guess, tau_guess = \
              heuristic_peak_parameter_guesses(
                x_foot1,
                x_peak, y_peak,
                x_foot2,
              )
      y_peak_fit_guess = \
              peak_function(x_data, h_guess, mu_guess, sigma_guess, tau_guess)
      y_peak_fit_guess_with_background = \
              y_peak_fit_guess + y_background_fit_guess
      
      h1_guess = mu1_guess = sigma1_guess = tau1_guess = None
      y_peak1_fit_guess_with_background = None
      
      h2_guess = mu2_guess = sigma2_guess = tau2_guess = None
      y_peak2_fit_guess_with_background = None
    
    ################################
    # Guess plots
    ################################
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(
      x_data,
      y_background_fit_guess,
      label='\n'.join(
        [
          'background guess',
          f'  m={m_guess:.4}',
          f'  b={b_guess:.4}',
        ]
      ),
      linestyle='dotted',
    )
    axes.plot(peak_valley_x_locations, peak_valley_y_locations, 'rx')
    if y_peak_fit_guess_with_background is not None:
      axes.plot(
        x_data,
        y_peak_fit_guess_with_background,
        label='\n'.join(
          [
            'peak guess',
            f'  h={h_guess:.4}',
            f'  μ={mu_guess:.4}',
            f'  σ={sigma_guess:.4}',
            f'  τ={tau_guess:.4}',
          ]
        ),
        linestyle='dotted',
      )
    if y_peak1_fit_guess_with_background is not None:
      axes.plot(
        x_data,
        y_peak1_fit_guess_with_background,
        label='\n'.join(
          [
            'peak guess',
            f'  h1={h1_guess:.4}',
            f'  μ1={mu1_guess:.4}',
            f'  σ1={sigma1_guess:.4}',
            f'  τ1={tau1_guess:.4}',
          ]
        ),
        linestyle='dotted',
      )
    if y_peak2_fit_guess_with_background is not None:
      axes.plot(
        x_data,
        y_peak2_fit_guess_with_background,
        label='\n'.join(
          [
            'peak guess',
            f'  h2={h2_guess:.4}',
            f'  μ2={mu2_guess:.4}',
            f'  σ2={sigma2_guess:.4}',
            f'  τ2={tau2_guess:.4}',
          ]
        ),
        linestyle='dotted',
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

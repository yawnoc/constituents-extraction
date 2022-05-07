#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
from frame vs intensity data in `data/`.
"""


import csv
import collections
import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
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


def one_peak_model(x, m, b, h, mu, sigma, tau):
  return background_function(x, m, b) + peak_function(x, h, mu, sigma, tau)


def two_peak_model(x, m, b, h1, mu1, sigma1, tau1, h2, mu2, sigma2, tau2):
  return (
    background_function(x, m, b)
      + peak_function(x, h1, mu1, sigma1, tau1)
      + peak_function(x, h2, mu2, sigma2, tau2)
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


def heuristic_valley_side_foot_location(x_peak, y_peak, x_valley, y_valley):
  """
  Heuristically guess the location of the valley-side foot of a peak.
  
  This is the foot location that would be visible
  if the other peak were removed.
  
  1. Straight-line extrapolates to y = 0.
  2. Makes sure that the foot isn't too far away from the peak.
  """
  
  extrapolated_foot_displacement = \
          (x_valley - x_peak) / (y_peak - y_valley) * y_peak
  
  extrapolated_foot_direction = np.sign(extrapolated_foot_displacement)
  
  extrapolated_foot_distance = np.abs(extrapolated_foot_displacement)
  extrapolated_foot_distance = \
          min(
            extrapolated_foot_distance,
            1.5 * np.abs(x_valley - x_peak),
          )
  
  x_foot = x_peak + extrapolated_foot_direction * extrapolated_foot_distance
  
  return x_foot


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
              heuristic_valley_side_foot_location(
                x_peak1, y_peak1,
                x_valley, y_valley,
              )
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
              heuristic_valley_side_foot_location(
                x_peak2, y_peak2,
                x_valley, y_valley,
              )
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
      
      ################################
      # 2-peak fit
      ################################
      (
        m_fit, b_fit,
        h1_fit, mu1_fit, sigma1_fit, tau1_fit,
        h2_fit, mu2_fit, sigma2_fit, tau2_fit,
      ), _ = \
              opt.curve_fit(
                two_peak_model,
                x_data, y_data,
                (
                  m_guess, b_guess,
                  h1_guess, mu1_guess, sigma1_guess, tau1_guess,
                  h2_guess, mu2_guess, sigma2_guess, tau2_guess,
                )
              )
      y_fit = \
              two_peak_model(
                x_data,
                m_fit, b_fit,
                h1_fit, mu1_fit, sigma1_fit, tau1_fit,
                h2_fit, mu2_fit, sigma2_fit, tau2_fit,
              )
      
      y_background_fit = background_function(x_data, m_fit, b_fit)
      y_peak1_fit = \
              peak_function(x_data, h1_fit, mu1_fit, sigma1_fit, tau1_fit)
      y_peak1_fit_with_background = y_peak1_fit + y_background_fit
      y_peak2_fit = \
              peak_function(x_data, h2_fit, mu2_fit, sigma2_fit, tau2_fit)
      y_peak2_fit_with_background = y_peak2_fit + y_background_fit
      
      h_fit = mu_fit = sigma_fit = tau_fit = None
      y_peak_fit_with_background = None
    
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
      # 1-peak fit
      ################################
      (m_fit, b_fit, h_fit, mu_fit, sigma_fit, tau_fit), _ = \
              opt.curve_fit(
                one_peak_model,
                x_data, y_data,
              )
      y_fit = \
              one_peak_model(
                x_data,
                m_fit, b_fit,
                h_fit, mu_fit, sigma_fit, tau_fit,
              )
      
      y_background_fit = background_function(x_data, m_fit, b_fit)
      y_peak_fit = peak_function(x_data, h_fit, mu_fit, sigma_fit, tau_fit)
      y_peak_fit_with_background = y_peak_fit + y_background_fit
      
      h1_fit = mu1_fit = sigma1_fit = tau1_fit = None
      h2_fit = mu2_fit = sigma2_fit = tau2_fit = None
      y_peak1_fit_with_background = y_peak2_fit_with_background = None
    
    ################################
    # Guess plots
    ################################
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(
      x_data, y_background_fit_guess,
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
        x_data, y_peak_fit_guess_with_background,
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
        x_data, y_peak1_fit_guess_with_background,
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
        x_data, y_peak2_fit_guess_with_background,
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
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'guess-{file_name}.pdf'))
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'guess-{file_name}.png'))
    
    ################################
    # Fit plots
    ################################
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(x_data, y_fit, label='fit')
    axes.plot(
      x_data, y_background_fit,
      label='\n'.join(
        [
          'background fit',
          f'  m={m_fit:.4}',
          f'  b={b_fit:.4}',
        ]
      ),
      linestyle='dotted',
    )
    if y_peak_fit_with_background is not None:
      axes.plot(
        x_data, y_peak_fit_with_background,
        label='\n'.join(
          [
            'peak fit',
            f'  h={h_fit:.4}',
            f'  μ={mu_fit:.4}',
            f'  σ={sigma_fit:.4}',
            f'  τ={tau_fit:.4}',
          ]
        ),
        linestyle='dotted',
      )
    if y_peak1_fit_with_background is not None:
      axes.plot(
        x_data, y_peak1_fit_with_background,
        label='\n'.join(
          [
            'peak1 fit',
            f'  h1={h1_fit:.4}',
            f'  μ1={mu1_fit:.4}',
            f'  σ1={sigma1_fit:.4}',
            f'  τ1={tau1_fit:.4}',
          ]
        ),
        linestyle='dotted',
      )
    if y_peak2_fit_with_background is not None:
      axes.plot(
        x_data, y_peak2_fit_with_background,
        label='\n'.join(
          [
            'peak2 fit',
            f'  h2={h2_fit:.4}',
            f'  μ2={mu2_fit:.4}',
            f'  σ2={sigma2_fit:.4}',
            f'  τ2={tau2_fit:.4}',
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
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.pdf'))
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.png'))
    
    csv_file_name = os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.csv')
    with open(csv_file_name, 'w', encoding='utf-8', newline='') as csv_file:
      
      csv_writer = csv.writer(csv_file)
      
      csv_writer.writerow(['# Background', None])
      csv_writer.writerow(['m', m_fit])
      csv_writer.writerow(['b', b_fit])
      csv_writer.writerow([])
      
      if h_fit is not None:
        csv_writer.writerow(['# Peak', None])
        csv_writer.writerow(['h', h_fit])
        csv_writer.writerow(['mu', mu_fit])
        csv_writer.writerow(['sigma', sigma_fit])
        csv_writer.writerow(['tau', tau_fit])
        csv_writer.writerow([])
      
      if h1_fit is not None:
        csv_writer.writerow(['# Peak 1', None])
        csv_writer.writerow(['h1', h1_fit])
        csv_writer.writerow(['mu1', mu1_fit])
        csv_writer.writerow(['sigma1', sigma1_fit])
        csv_writer.writerow(['tau1', tau1_fit])
        csv_writer.writerow([])
      
      if h2_fit is not None:
        csv_writer.writerow(['# Peak 2', None])
        csv_writer.writerow(['h2', h2_fit])
        csv_writer.writerow(['mu2', mu2_fit])
        csv_writer.writerow(['sigma2', sigma2_fit])
        csv_writer.writerow(['tau2', tau2_fit])
        csv_writer.writerow([])


if __name__ == '__main__':
  main()

#!/usr/bin/env python3

"""
# fit.py

Do curve fitting to extract aggregate and monomer constituents
from frame vs intensity data for biological samples.

I (@yawnoc) am not a biologist,
so I am just going to call them `peak1` and `peak2`.

---

Copyright 2022 Conway
Licensed under MIT-0, see LICENSE.
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


def peak_function(x, a, mu, sigma, tau):
  """
  Exponentially modified Gaussian for each peak in the intensity
  (see `exponentially-modified-gaussian.md` for a derivation):
          f(x) =
                  A / (2 tau)
                  . exp[1/2 (sigma/tau)^2 - (x - mu)/tau]
                  . erfc[1/sqrt(2) (sigma/tau - (x - mu)/sigma)],
  where
          A is is the area under the curve,
          mu is the mean of the normal distribution,
          sigma is the standard deviation of the normal distribution,
          tau is the time scale of the exponential distribution.
  """
  return (
    a / (2 * tau)
      * np.exp(1/2 * (sigma/tau)**2 - (x - mu)/tau)
      * sp.erfc(1/np.sqrt(2) * (sigma/tau - (x - mu)/sigma))
  )


def one_peak_model(x, m, b, a, mu, sigma, tau):
  return background_function(x, m, b) + peak_function(x, a, mu, sigma, tau)


def two_peak_model(x, m, b, a1, mu1, sigma1, tau1, a2, mu2, sigma2, tau2):
  return (
    background_function(x, m, b)
      + peak_function(x, a1, mu1, sigma1, tau1)
      + peak_function(x, a2, mu2, sigma2, tau2)
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
  
  No rigorous mathematical theory underlies these guesses,
  which are based purely on intuition and experimentation.
  """
  
  a_guess = y_peak * np.sqrt((x_foot2 - x_peak) * (x_peak - x_foot1))
  mu_guess = x_peak - 0.2 * (x_peak - x_foot1)
  sigma_guess = 0.4 * (x_peak - x_foot1)
  tau_guess = 0.33 * (x_foot2 - x_peak)
  
  return a_guess, mu_guess, sigma_guess, tau_guess


def approximate_fraction_error(a1_fit, a2_fit, a1_fit_error, a2_fit_error):
  """
  Approximate the standard error for A1/(A1+A2) and A2/(A1+A2).
  
  Applying the physicist's propagation of errors for f1 = A1/(A1+A2):
          s[f1]^2 ≃ (∂f1/∂A1)^2 s[A1]^2 + (∂f1/∂A2)^2 s[A2]^2
                  = (A2^2 s[A1]^2 + A1^2 s[A2]^2) / (A1 + A2)^4
          s[f1] = sqrt(A2^2 s[A1]^2 + A1^2 s[A2]^2) / (A1 + A2)^2.
  The result is the same for the complement f2 = A2/(A1+A2).
  """
  return (
    np.sqrt(a2_fit**2 * a1_fit_error**2 + a1_fit**2 * a2_fit_error**2)
      /
    (a1_fit + a2_fit)**2
  )


LOCATION_MARKERS_STYLE = 'rx'

BACKGROUND_CURVE_COLOUR = '#2CA02C'
PEAK1_CURVE_COLOUR = PEAK_CURVE_COLOUR = '#9467BD'
PEAK2_CURVE_COLOUR = '#D62728'

BACKGROUND_LINE_STYLE = 'dotted'


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
      
      a_guess = mu_guess = sigma_guess = tau_guess = None
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
      a1_guess, mu1_guess, sigma1_guess, tau1_guess = \
              heuristic_peak_parameter_guesses(
                x_peak1_foot1,
                x_peak1, y_peak1,
                x_peak1_foot2,
              )
      y_peak1_fit_guess = \
              peak_function(
                x_data,
                a1_guess, mu1_guess, sigma1_guess, tau1_guess,
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
      a2_guess, mu2_guess, sigma2_guess, tau2_guess = \
              heuristic_peak_parameter_guesses(
                x_peak2_foot1,
                x_peak2, y_peak2,
                x_peak2_foot2,
              )
      y_peak2_fit_guess = \
              peak_function(
                x_data,
                a2_guess, mu2_guess, sigma2_guess, tau2_guess,
              )
      y_peak2_fit_guess_with_background = \
              y_peak2_fit_guess + y_background_fit_guess
      
      ################################
      # 2-peak fit
      ################################
      fit_parameters, fit_covariances = \
              opt.curve_fit(
                two_peak_model,
                x_data, y_data,
                (
                  m_guess, b_guess,
                  a1_guess, mu1_guess, sigma1_guess, tau1_guess,
                  a2_guess, mu2_guess, sigma2_guess, tau2_guess,
                )
              )
      (
        m_fit, b_fit,
        a1_fit, mu1_fit, sigma1_fit, tau1_fit,
        a2_fit, mu2_fit, sigma2_fit, tau2_fit,
      ) = \
              fit_parameters
      (
        m_fit_error, b_fit_error,
        a1_fit_error, mu1_fit_error, sigma1_fit_error, tau1_fit_error,
        a2_fit_error, mu2_fit_error, sigma2_fit_error, tau2_fit_error,
      ) = \
              np.sqrt(np.diag(fit_covariances))
      y_fit = \
              two_peak_model(
                x_data,
                m_fit, b_fit,
                a1_fit, mu1_fit, sigma1_fit, tau1_fit,
                a2_fit, mu2_fit, sigma2_fit, tau2_fit,
              )
      
      y_background_fit = background_function(x_data, m_fit, b_fit)
      y_peak1_fit = \
              peak_function(x_data, a1_fit, mu1_fit, sigma1_fit, tau1_fit)
      y_peak1_fit_with_background = y_peak1_fit + y_background_fit
      y_peak2_fit = \
              peak_function(x_data, a2_fit, mu2_fit, sigma2_fit, tau2_fit)
      y_peak2_fit_with_background = y_peak2_fit + y_background_fit
      
      a_fit = mu_fit = sigma_fit = tau_fit = None
      a_fit_error = mu_fit_error = sigma_fit_error = tau_fit_error = None
      y_peak_fit_with_background = None
    
    except ValueError: # 1 peak
      
      x_foot1, y_foot1, x_peak, y_peak, x_foot2, y_foot2 = \
              heuristic_peak_valley_locations(x_data, y_data, m_guess, b_guess)
      peak_valley_x_locations = [x_foot1, x_peak, x_foot2]
      peak_valley_y_locations = [y_foot1, y_peak, y_foot2]
      
      ################################
      # Peak guess
      ################################
      a_guess, mu_guess, sigma_guess, tau_guess = \
              heuristic_peak_parameter_guesses(
                x_foot1,
                x_peak, y_peak,
                x_foot2,
              )
      y_peak_fit_guess = \
              peak_function(x_data, a_guess, mu_guess, sigma_guess, tau_guess)
      y_peak_fit_guess_with_background = \
              y_peak_fit_guess + y_background_fit_guess
      
      a1_guess = mu1_guess = sigma1_guess = tau1_guess = None
      y_peak1_fit_guess_with_background = None
      
      a2_guess = mu2_guess = sigma2_guess = tau2_guess = None
      y_peak2_fit_guess_with_background = None
      
      ################################
      # 1-peak fit
      ################################
      fit_parameters, fit_covariances = \
              opt.curve_fit(
                one_peak_model,
                x_data, y_data,
                (m_guess, b_guess, a_guess, mu_guess, sigma_guess, tau_guess),
              )
      m_fit, b_fit, a_fit, mu_fit, sigma_fit, tau_fit = fit_parameters
      (
        m_fit_error, b_fit_error,
        a_fit_error, mu_fit_error, sigma_fit_error, tau_fit_error,
      ) = \
              np.sqrt(np.diag(fit_covariances))
      y_fit = \
              one_peak_model(
                x_data,
                m_fit, b_fit,
                a_fit, mu_fit, sigma_fit, tau_fit,
              )
      
      y_background_fit = background_function(x_data, m_fit, b_fit)
      y_peak_fit = peak_function(x_data, a_fit, mu_fit, sigma_fit, tau_fit)
      y_peak_fit_with_background = y_peak_fit + y_background_fit
      
      a1_fit = mu1_fit = sigma1_fit = tau1_fit = None
      a2_fit = mu2_fit = sigma2_fit = tau2_fit = None
      a1_fit_error = mu1_fit_error = sigma1_fit_error = tau1_fit_error = None
      a2_fit_error = mu2_fit_error = sigma2_fit_error = tau2_fit_error = None
      y_peak1_fit_with_background = y_peak2_fit_with_background = None
    
    ################################
    # Guess plots
    ################################
    figure, axes = plt.subplots()
    axes.plot(x_data, y_data, label='data')
    axes.plot(
      x_data, y_background_fit_guess,
      BACKGROUND_CURVE_COLOUR,
      label='\n'.join(
        [
          'background guess',
          f'  m={m_guess:.4}',
          f'  b={b_guess:.4}',
        ]
      ),
      linestyle=BACKGROUND_LINE_STYLE,
    )
    axes.plot(
      peak_valley_x_locations, peak_valley_y_locations,
      LOCATION_MARKERS_STYLE,
    )
    if y_peak_fit_guess_with_background is not None:
      axes.plot(
        x_data, y_peak_fit_guess_with_background,
        PEAK_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak guess',
            f'  A={a_guess:.4}',
            f'  μ={mu_guess:.4}',
            f'  σ={sigma_guess:.4}',
            f'  τ={tau_guess:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    if y_peak1_fit_guess_with_background is not None:
      axes.plot(
        x_data, y_peak1_fit_guess_with_background,
        PEAK1_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak1 guess',
            f'  A1={a1_guess:.4}',
            f'  μ1={mu1_guess:.4}',
            f'  σ1={sigma1_guess:.4}',
            f'  τ1={tau1_guess:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    if y_peak2_fit_guess_with_background is not None:
      axes.plot(
        x_data, y_peak2_fit_guess_with_background,
        PEAK2_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak2 guess',
            f'  A2={a2_guess:.4}',
            f'  μ2={mu2_guess:.4}',
            f'  σ2={sigma2_guess:.4}',
            f'  τ2={tau2_guess:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    axes.set(
      title=f'{file_name} (guess)',
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
      BACKGROUND_CURVE_COLOUR,
      label='\n'.join(
        [
          'background fit',
          f'  m={m_fit:.4}',
          f'  b={b_fit:.4}',
        ]
      ),
      linestyle=BACKGROUND_LINE_STYLE,
    )
    if y_peak_fit_with_background is not None:
      axes.plot(
        x_data, y_peak_fit_with_background,
        PEAK_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak fit',
            f'  A={a_fit:.4}',
            f'  μ={mu_fit:.4}',
            f'  σ={sigma_fit:.4}',
            f'  τ={tau_fit:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    if y_peak1_fit_with_background is not None:
      axes.plot(
        x_data, y_peak1_fit_with_background,
        PEAK1_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak1 fit',
            f'  A1={a1_fit:.4}',
            f'  μ1={mu1_fit:.4}',
            f'  σ1={sigma1_fit:.4}',
            f'  τ1={tau1_fit:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    if y_peak2_fit_with_background is not None:
      axes.plot(
        x_data, y_peak2_fit_with_background,
        PEAK2_CURVE_COLOUR,
        label='\n'.join(
          [
            'peak2 fit',
            f'  A2={a2_fit:.4}',
            f'  μ2={mu2_fit:.4}',
            f'  σ2={sigma2_fit:.4}',
            f'  τ2={tau2_fit:.4}',
          ]
        ),
        linestyle=BACKGROUND_LINE_STYLE,
      )
    axes.set(
      title=f'{file_name} (fit)',
      xlabel=f'Normalised frame number [{int(frame_min)}, {int(frame_max)}]',
      ylabel=f'Normalised intensity [{intensity_min:.3}, {intensity_max:.3}]',
    )
    axes.legend()
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.pdf'))
    plt.savefig(os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.png'))
    
    tsv_file_name = os.path.join(OUTPUT_DIRECTORY, f'fit-{file_name}.tsv')
    with open(tsv_file_name, 'w', encoding='utf-8', newline='') as tsv_file:
      
      tsv_writer = csv.writer(tsv_file, delimiter='\t')
      
      tsv_writer.writerow(['# Unscaled ranges'])
      tsv_writer.writerow(['frame_min', int(frame_min)])
      tsv_writer.writerow(['frame_max', int(frame_max)])
      tsv_writer.writerow(['intensity_min', intensity_min])
      tsv_writer.writerow(['intensity_max', intensity_max])
      tsv_writer.writerow([])
      
      tsv_writer.writerow(
        [
          '# Background parameter',
          'Estimate',
          'Standard error',
          'Relative standard error',
        ]
      )
      tsv_writer.writerow(['m', m_fit, m_fit_error, m_fit_error/m_fit])
      tsv_writer.writerow(['b', b_fit, b_fit_error, b_fit_error/b_fit])
      tsv_writer.writerow([])
      
      if a_fit is not None:
        tsv_writer.writerow(
          [
            '# Peak parameter',
            'Estimate',
            'Standard error',
            'Relative standard error',
          ]
        )
        tsv_writer.writerow(['A', a_fit, a_fit_error, a_fit_error/a_fit])
        tsv_writer.writerow(['mu', mu_fit, mu_fit_error, mu_fit_error/mu_fit])
        tsv_writer.writerow(
          ['sigma', sigma_fit, sigma_fit_error, sigma_fit_error/sigma_fit]
        )
        tsv_writer.writerow(
          ['tau', tau_fit, tau_fit_error, tau_fit_error/tau_fit]
        )
        tsv_writer.writerow([])
      
      if a1_fit is not None:
        tsv_writer.writerow(
          [
            '# Peak 1 parameter',
            'Estimate',
            'Standard error',
            'Relative standard error',
          ]
        )
        tsv_writer.writerow(['A1', a1_fit, a1_fit_error, a1_fit_error/a1_fit])
        tsv_writer.writerow(
          ['mu1', mu1_fit, mu1_fit_error, mu1_fit_error/mu1_fit]
        )
        tsv_writer.writerow(
          ['sigma1', sigma1_fit, sigma1_fit_error, sigma1_fit_error/sigma1_fit]
        )
        tsv_writer.writerow(
          ['tau1', tau1_fit, tau1_fit_error, tau1_fit_error/tau1_fit]
        )
        tsv_writer.writerow([])
      
      if a2_fit is not None:
        tsv_writer.writerow(
          [
            '# Peak 2 parameter',
            'Estimate',
            'Standard Error',
          ]
        )
        tsv_writer.writerow(['A2', a2_fit, a2_fit_error, a2_fit_error/a2_fit])
        tsv_writer.writerow(
          ['mu2', mu2_fit, mu2_fit_error, mu2_fit_error/mu2_fit]
        )
        tsv_writer.writerow(
          ['sigma2', sigma2_fit, sigma2_fit_error, sigma2_fit_error/sigma2_fit]
        )
        tsv_writer.writerow(
          ['tau2', tau2_fit, tau2_fit_error, tau2_fit_error/tau2_fit]
        )
        tsv_writer.writerow([])
      
      if a1_fit is not None and a2_fit is not None:
        tsv_writer.writerow(
          [
            '# Area fraction',
            'Estimate',
            'Standard error (approximate)',
            'Relative standard error (approximate)',
          ]
        )
        a1_fraction = a1_fit / (a1_fit + a2_fit)
        a2_fraction = a2_fit / (a1_fit + a2_fit)
        a1_fraction_error = a2_fraction_error = \
                approximate_fraction_error(
                  a1_fit, a2_fit,
                  a1_fit_error, a2_fit_error,
                )
        tsv_writer.writerow(
          [
            'A1/(A1+A2)',
            a1_fraction,
            a1_fraction_error,
            a1_fraction_error/a1_fraction,
          ]
        )
        tsv_writer.writerow(
          [
            'A2/(A1+A2)',
            a2_fraction,
            a2_fraction_error,
            a2_fraction_error/a2_fraction,
          ]
        )
        tsv_writer.writerow([])


if __name__ == '__main__':
  main()

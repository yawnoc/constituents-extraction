#!/usr/bin/env python3

"""
# fit.py

Do fitting to extract Monomer and Aggregate constituents
for {frame number} vs {average intensity} data.
"""


import os

import matplotlib.pyplot as plt
import numpy as np
import scipy.optimize as opt
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


def peak_model(x, a, k, mu, sigma):
  """
  Exponentially modified Gaussian for each peak of normalised intensity.
  
  Here, a is an amplitude parameter with have tacked on.
  The exponentially modified Gaussian distribution has probability density
    f(x, K) = 1/(2K) exp(1/(2K^2) - x/K) erfc(-(x - 1/K) / sqrt(2)),
  and transforms to the Wikipedia version according to K = 1/(sigma lambda).
  See <https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.exponnorm.html>
  and <https://en.wikipedia.org/wiki/Exponentially_modified_Gaussian_distribution>.
  """
  return a * st.exponnorm.pdf(x, k, mu, sigma)


def full_model(x, m, b, a1, k1, mu1, sigma1, a2, k2, mu2, sigma2):
  return (
    background_model(x, m, b)
    + peak_model(x, a1, k1, mu1, sigma1)
    + peak_model(x, a2, k2, mu2, sigma2)
  )


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


def de_normalise(intensity_min, intensity_max, y_data):
  """
  De-normalise frame and intensity data.
  """
  intensity_data = intensity_min + (intensity_max - intensity_min) * y_data
  return intensity_data


def fit(x_data, y_data):
  
  m_bounds = [0, 0.2]
  b_bounds = [0, 0.05]
  a_bounds = [0, 1]
  k_bounds = [0.1, 4]
  mu_bounds = [0.1, 4]
  sigma_bounds = [0.1, 4]
  parameter_bounds = [
    m_bounds, b_bounds,
    a_bounds, k_bounds, mu_bounds, sigma_bounds,
    a_bounds, k_bounds, mu_bounds, sigma_bounds,
  ]
  
  def sum_of_squares_error(parameters):
    y_model = full_model(x_data, *parameters)
    return np.sum((y_data - y_model) ** 2)
  
  initial_parameters = \
          opt.differential_evolution(
            sum_of_squares_error,
            parameter_bounds,
            seed=0,
          ).x
  fitted_parameters, _ = \
          opt.curve_fit(full_model, x_data, y_data, initial_parameters)
  
  return fitted_parameters


def make_plot(file_name, data):
  
  frame_data = data[:, 0]
  intensity_data = data[:, 1]
  
  frame_min, frame_max, intensity_min, intensity_max, x_data, y_data = \
          normalise(frame_data, intensity_data)
  fitted_parameters = fit(x_data, y_data)
  y_data_fitted = full_model(x_data, *fitted_parameters)
  intensity_data_fitted = \
          de_normalise(intensity_min, intensity_max, y_data_fitted)
  
  figure, axes = plt.subplots()
  axes.plot(frame_data, intensity_data)
  axes.plot(frame_data, intensity_data_fitted)
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

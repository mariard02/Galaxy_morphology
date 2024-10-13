# Galaxy Surface Brightness Profile Fitting

This project provides a Python class `galaxy` for analyzing the radial surface brightness profile of galaxies. It fits Sersic models to observational data, offering tools for both single-component and two-component Sersic profile fitting, with visualization and error handling. It was created for the practical exercises of the subject "Galaxies and Cosmology: an introduction", taught at the University of Geneva.

## Features

- Load and visualize surface brightness (SB) profiles from observational data.
- Perform single-component Sersic profile fitting using `astropy`'s `Sersic1D` model.
- Perform two-component Sersic profile fitting for galaxies with more complex structures (e.g., bulge + disk).
- Generate plots comparing observational data with best-fit Sersic models.
- Residual plots for evaluating the goodness of fit.

## Requirements

- Python 3.x
- NumPy
- Matplotlib
- Astropy

## Data Format

The input data for each galaxy should be stored in a text file containing three columns:

- Arcsecond: Angular radius of the galaxy in arcseconds.
- Surface Brightness (SB): Observed surface brightness at that radius (in arbitrary units or counts per arcsecond squared).
- Surface Brightness Error (SB_err): Error in the surface brightness measurement.
- Each file should have no headers and should be formatted as follows:

arcsecond SB SB_err
0.10      21.2    0.3
0.15      20.8    0.2
...

## Usage

1. Initialization

To create a galaxy object, initialize it with the galaxy's name, distance (in megaparsecs), and the path to the data file:

```python
from galaxy_class import galaxy

# Example: create a galaxy object
gal = galaxy(name="NGC 1234", distance=10.5, path="data/ngc1234.txt")
```

2. Display basic information

The galaxy_info method will print the basic properties of the galaxy:

```python
gal.galaxy_info()
# Output:
# Catalogue name: NGC 1234
# Distance = 10.5 Mpc
```

3. Plot the radial profile

Use the galaxy_plot method to visualize the galaxy's radial surface brightness profile:

```python
gal.galaxy_plot()
```

4. Single-Component Sersic Fitting

The galaxy_fit1 method fits a single Sersic profile to the data. You can provide initial guesses for the amplitude (a), effective radius (r), and Sersic index (nm):

```python
gal.galaxy_fit1(a=0.5, r=1, nm=4)
```

This will perform the fit and show a plot comparing the observed profile with the best-fit model, as well as residuals.

5. Two-Component Sersic Fitting

If the galaxy has more complex features (e.g., a bulge and disk), the galaxy_fit2 method fits a two-component Sersic profile. You can provide initial guesses for the amplitude, effective radius, and Sersic index for each component:

```python
gal.galaxy_fit2(a1=0.17, r1=2.6, n1=1, a2=0.1, r2=3.4, n2=4.2)
```

This will also show the best-fit model and the residuals.

## Example

Here is a complete example of how to use the class:

```python
from galaxy_class import galaxy

# Create a galaxy object
gal = galaxy(name="NGC 1234", distance=10.5, path="data/ngc1234.txt")

# Display basic information
gal.galaxy_info()

# Plot the radial surface brightness profile
gal.galaxy_plot()

# Fit a single Sersic model
gal.galaxy_fit1(a=0.5, r=1, nm=4)

# Fit a two-component Sersic model
gal.galaxy_fit2(a1=0.17, r1=2.6, n1=1, a2=0.1, r2=3.4, n2=4.2)

```


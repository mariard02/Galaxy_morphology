import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D

class galaxy(object):
    """
    Class representing a galaxy. This class handles loading observational data,
    fitting radial surface brightness profiles with Sersic models, and visualizing the results.
    """
    def __init__(self, name, distance, path):
        """
        Initialize the galaxy object with its name, distance, and the path to the data file.
        Data file should contain the surface brightness profile in arcseconds, brightness, and errors.

        Parameters:
        name (str): Name of the galaxy.
        distance (float): Distance to the galaxy in Mpc.
        path (str): Path to the file containing galaxy surface brightness data.
        """
        self.name = name  # Name of the galaxy
        self.distance = distance  # Distance to the galaxy in Mpc
        self.path = path  # Path to the data file

        # Load the data: arcsec (angular radius), SB (surface brightness), and SB_err (errors in surface brightness)
        self.arcsec, self.SB, self.SB_err = np.loadtxt(self.path, usecols=[0, 1, 2], unpack=True)

        # Convert angular radius in arcseconds to physical radius in kiloparsecs (kpc)
        self.radius = self.distance * self.arcsec / 206

    def galaxy_info(self):
        """
        Prints basic information about the galaxy, including its name and distance.
        """
        print("Catalogue name: {}".format(self.name))
        print("Distance = {} Mpc".format(self.distance))

    def galaxy_plot(self):
        """
        Plots the radial surface brightness profile of the galaxy with error bars.
        The plot uses a logarithmic scale for surface brightness.
        """
        fig = plt.figure(figsize=(8, 6))
        plt.errorbar(self.radius, self.SB, yerr=self.SB_err, lw=2, markersize=4, fmt='o', color='deepskyblue', capsize=3)
        plt.xlabel(r'Radius [kpc]', fontsize=12)
        plt.ylabel(r'Surface Brightness', fontsize=12)
        plt.xticks(size=12)
        plt.yticks(size=12)
        plt.semilogy()  # Use a logarithmic scale for the y-axis (SB)
        plt.title("Radial profile for {}".format(self.name))
        plt.show()

    def galaxy_fit1(self, a=0.5, r=1, nm=4):
        """
        Fits a single-component Sersic model to the galaxy's radial surface brightness profile.
        
        Parameters:
        a (float): Initial guess for the amplitude of the Sersic function.
        r (float): Initial guess for the effective radius of the Sersic function.
        nm (float): Initial guess for the Sersic index (n parameter, controls the shape of the profile).
        """
        # Initialize a single Sersic model with initial parameter guesses
        Sersic_init = Sersic1D(amplitude=a, r_eff=r, n=nm)

        # Define the fitting algorithm using the least-squares fitter (Levenberg-Marquardt)
        fit = fitting.LevMarLSQFitter()

        # Perform the fit, considering uncertainties (weights = 1/SB_err)
        best_Sersic = fit(Sersic_init, self.radius, self.SB, weights=1/self.SB_err, maxiter=500000)

        # Print the best-fit parameters
        print("")
        print(best_Sersic)

        # Print the uncertainties in the fitted parameters (diagonal of the covariance matrix)
        cov_diag = np.sqrt(np.diag(fit.fit_info['param_cov']))
        print('error =', cov_diag)

        # Create a plot to compare data with the best-fit model
        fig, ax = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 7))

        fig.suptitle("Best fit model for the radial profile of {}".format(self.name))

        # Plot the best-fit Sersic model and the data
        ax[0].plot(self.radius, best_Sersic(self.radius), label='Best-fit Sersic Model', color='grey', lw=2)
        ax[0].errorbar(self.radius, self.SB, yerr=self.SB_err, fmt='o', label='Data', color='deepskyblue', markersize=4, capsize=3)
        ax[0].set_xlabel('Radius (kpc)')
        ax[0].set_ylabel('Surface Brightness [counts/arcsec^2]')
        ax[0].legend()
        ax[0].semilogy()  # Use a logarithmic scale for the surface brightness

        # Plot the residuals (difference between data and the best-fit model)
        ax[1].scatter(self.radius, (self.SB - best_Sersic(self.radius))/self.SB, color='deepskyblue')
        ax[1].set_xlabel('Radius (kpc)')
        ax[1].set_ylabel('Residuals')

        plt.tight_layout()
        plt.show()

    def galaxy_fit2(self, a1=0.17, r1=2.6, n1=1, a2=0.1, r2=3.4, n2=4.2):
        """
        Fits a two-component Sersic model to the galaxy's radial surface brightness profile.
        This is useful for galaxies that may have a bulge and disk component.

        Parameters:
        a1, r1, n1 (float): Initial guesses for the amplitude, effective radius, and index for the first Sersic component.
        a2, r2, n2 (float): Initial guesses for the amplitude, effective radius, and index for the second Sersic component.
        """
        # Filter out invalid data (e.g., NaN or infinite values) in radius, SB, or SB_err
        mask = np.isfinite(self.radius) & np.isfinite(self.SB) & np.isfinite(self.SB_err)
        
        # Apply the mask to the data
        radius_filtered = self.radius[mask]
        SB_filtered = self.SB[mask]
        SB_err_filtered = self.SB_err[mask]

        # Initialize two Sersic models with initial guesses for each component
        Sersic1_init = Sersic1D(amplitude=a1, r_eff=r1, n=n1)  # First Sersic component
        Sersic2_init = Sersic1D(amplitude=a2, r_eff=r2, n=n2)  # Second Sersic component

        # Sum of the two Sersic components
        Sersic_total_init = Sersic1_init + Sersic2_init

        # Define the fitting algorithm (Levenberg-Marquardt least-squares fitter)
        fit = fitting.LevMarLSQFitter()

        # Perform the fit on the filtered data, considering uncertainties (weights = 1/SB_err)
        try:
            best_Sersic_total = fit(Sersic_total_init, radius_filtered, SB_filtered, weights=1/SB_err_filtered, maxiter=5000)
        except Exception as e:
            print(f"Error during fitting: {e}")
            return

        # Print the best-fit parameters for the two-component Sersic model
        print("")
        print(best_Sersic_total)

        # Print the uncertainties in the fitted parameters (diagonal of the covariance matrix)
        cov_diag = np.sqrt(np.diag(fit.fit_info['param_cov']))
        print('error =', cov_diag)

        # Create a plot to compare the data with the best-fit two-component Sersic model
        fig, ax = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, figsize=(7, 7))

        fig.suptitle("Best fit model for the radial profile of {}".format(self.name))

        # Plot the best-fit two-component Sersic model and the data
        ax[0].plot(radius_filtered, best_Sersic_total(radius_filtered), label='Best-fit Sersic Model (two components)', color='grey', lw=2)
        ax[0].errorbar(radius_filtered, SB_filtered, yerr=SB_err_filtered, fmt='o', label='Data', color='deepskyblue', markersize=4, capsize=3)

        ax[0].set_xlabel('Radius (kpc)')
        ax[0].set_ylabel('Surface Brightness [counts/arcsec^2]')
        ax[0].legend()
        ax[0].semilogy()  # Use a logarithmic scale for the surface brightness

        # Plot the residuals (difference between data and the best-fit model)
        residuals = (SB_filtered - best_Sersic_total(radius_filtered)) / SB_filtered
        ax[1].scatter(radius_filtered, residuals, color='deepskyblue')
        ax[1].set_xlabel('Radius (kpc)')
        ax[1].set_ylabel('Residuals')

        plt.tight_layout()
        plt.show()

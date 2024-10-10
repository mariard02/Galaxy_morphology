import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import models, fitting
from astropy.modeling.models import Sersic1D

class galaxy(object):
	def __init__(self, name, distance, path):
		self.name = name
		self.distance = distance
		self.path = path

		self.arcsec, self.SB, self.SB_err = np.loadtxt(self.path, usecols=[0,1,2],unpack=True)

		self.radius = self.distance * self.arcsec / 206

	def galaxy_info(self):
		print("Catalogue name: {}".format(self.name))
		print("Distance = {} Mpc".format(self.distance))

	def galaxy_plot(self):
		fig = plt.figure(figsize = (8, 6))
		plt.errorbar(self.radius, self.SB, yerr=self.SB_err, lw=2, markersize=4, fmt='o', color='deepskyblue', capsize=3)
		plt.xlabel(r'radius [kpc]', fontsize=12) 
		plt.ylabel(r'SB', fontsize=12) 
		plt.xticks(size=12) 
		plt.yticks(size=12)
		plt.semilogy()
		plt.title("Radial profile for {}".format(self.name))

	def galaxy_fit1(self, a = 0.5, r = 1, nm = 4):
		# Initiate function with initial guess for the parameters

		Sersic_init = Sersic1D(amplitude = a, r_eff = r, n = nm)

		# Fitting function
		fit = fitting.LevMarLSQFitter()

		# Perform the fit (note that uncertainties are also considered)
		best_Sersic = fit(Sersic_init, self.radius, self.SB, weights = 1/self.SB_err, maxiter = 500000)

		print("")
		print(best_Sersic)

		cov_diag = np.sqrt(np.diag(fit.fit_info['param_cov']))
		print('error =', cov_diag)

		fig, ax = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, figsize = (7, 7))

		fig.suptitle("Best fit model for the radial profile of {}".format(self.name))

		ax[0].plot(self.radius, best_Sersic(self.radius), label='Best-fit Sersic Model', color='grey', lw=2)
		ax[0].errorbar(self.radius, self.SB, yerr=self.SB_err, fmt='o', label='Data', color='deepskyblue', markersize=4, capsize=3)
		ax[0].set_xlabel('Radius (Mpc)')
		ax[0].set_ylabel('Surface Brightness [counts/arcsec^2]')
		ax[0].legend()
		ax[0].semilogy()

		ax[1].scatter(self.radius, (self.SB - best_Sersic(self.radius))/self.SB, color = 'deepskyblue')
		ax[1].set_xlabel('Radius (kpc)')
		ax[1].set_ylabel('Residuals')
		#ax[1].semilogy()

	def galaxy_fit2(self, a1 = 0.17, r1 = 2.6, n1 = 1, a2 = 0.1, r2 = 3.4, n2 = 4.2):
	    # Filtrar los datos para eliminar NaN o inf en radius, SB o SB_err
	    mask = np.isfinite(self.radius) & np.isfinite(self.SB) & np.isfinite(self.SB_err)
	    
	    # Aplicar el filtro a los datos
	    radius_filtered = self.radius[mask]
	    SB_filtered = self.SB[mask]
	    SB_err_filtered = self.SB_err[mask]

	    # Iniciar las funciones de Sérsic con suposiciones iniciales para los parámetros
	    Sersic1_init = Sersic1D(amplitude=a1, r_eff=r1, n=n1)  # Primer componente Sérsic
	    Sersic2_init = Sersic1D(amplitude=a2, r_eff=r2, n=n2)  # Segundo componente Sérsic

	    # Suma de dos funciones de Sérsic
	    Sersic_total_init = Sersic1_init + Sersic2_init

	    # Función de ajuste
	    fit = fitting.LevMarLSQFitter()

	    # Realizar el ajuste con los datos filtrados
	    try:
	        best_Sersic_total = fit(Sersic_total_init, radius_filtered, SB_filtered, weights=1/SB_err_filtered, maxiter=5000)
	    except Exception as e:
	        print(f"Error durante el ajuste: {e}")
	        return

	    print("")
	    print(best_Sersic_total)

	    # Cálculo de las incertidumbres (diagonal de la matriz de covarianza)
	    cov_diag = np.sqrt(np.diag(fit.fit_info['param_cov']))
	    print('error =', cov_diag)

	    # Gráficos de los resultados
	    fig, ax = plt.subplots(2, gridspec_kw={'height_ratios': [3, 1]}, figsize = (7, 7))

	    fig.suptitle("Best fit model for the radial profile of {}".format(self.name))

	    # Graficar el modelo ajustado (suma de los dos perfiles de Sérsic)
	    ax[0].plot(radius_filtered, best_Sersic_total(radius_filtered), label='Best-fit Sersic Model (two components)', color='grey', lw=2)
	    
	    # Graficar los datos observacionales con barras de error
	    ax[0].errorbar(radius_filtered, SB_filtered, yerr=SB_err_filtered, fmt='o', label='Data', color='deepskyblue', markersize=4, capsize=3)

	    ax[0].set_xlabel('Radius (kpc)')
	    ax[0].set_ylabel('Surface Brightness [counts/arcsec^2]')
	    ax[0].legend()
	    ax[0].semilogy()

	    # Graficar los residuales
	    residuals = (SB_filtered - best_Sersic_total(radius_filtered)) / SB_filtered
	    ax[1].scatter(radius_filtered, residuals, color='deepskyblue')
	    ax[1].set_xlabel('Radius (Mpc)')
	    ax[1].set_ylabel('Residuals')

	    plt.tight_layout()
	    #plt.show()


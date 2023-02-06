from palmiche.umbrella import ic50
from palmiche.utils import tools
import numpy as np
from scipy.integrate import quad


class dG_corrector:
    def __init__(self,BoundRegion, cylinder_radium, cylinder_force_cte_restraint, orientation_theta_naught, absolute_temperature, PmfXvgPath, integrator = "simpson"):
        self.BoundRegion = BoundRegion
        self.cylinder_radium = cylinder_radium
        self.cylinder_force_cte_restraint = cylinder_force_cte_restraint
        self.orientation_theta_naught = orientation_theta_naught
        self.absolute_temperature = absolute_temperature
        self.PmfXvgPath = PmfXvgPath
        self.integrator = integrator

    def cylinder_correction(self):
        calculated_Ki = ic50.main(
            BoundRegion = self.BoundRegion,
            radium = self.cylinder_radium,
            sigma = 2*tools.get_sigma(self.absolute_temperature,self.cylinder_force_cte_restraint),
            absolute_temperature = self.absolute_temperature,
            PmfXvgPath = self.PmfXvgPath,
            integrator = self.integrator
        )
        delta_G_correction = -tools.CTE.R*self.absolute_temperature*np.log(calculated_Ki)
        return delta_G_correction

    def _I_theta_naught1_integrand(self, theta):
        return np.exp(-0.5*tools.beta(self.absolute_temperature)*(theta - self.orientation_theta_naught)**2)*np.sin(theta)

    def _I_theta_naught2_integrand(self, theta):
        return (theta - self.orientation_theta_naught)**2 * self._I_theta_naught1_integrand(theta)


    def orientation_correction(self):
        I_theta_naught1 = quad(self._I_theta_naught1_integrand,self.orientation_theta_naught, np.pi)[0]
        I_theta_naught2 = quad(self._I_theta_naught2_integrand,self.orientation_theta_naught, np.pi)[0]

        A_theta_naught = 1 / (2*np.pi*(I_theta_naught1 - np.cos(self.orientation_theta_naught) + 1))
        deltaS_corr = tools.CTE.Kb*np.log(4*np.pi) + tools.CTE.Kb*np.log(A_theta_naught) - A_theta_naught*np.pi*I_theta_naught2/self.absolute_temperature
        return deltaS_corr

    def __call__(self):
        dg_corrected = self.cylinder_correction() - self.absolute_temperature*self.orientation_correction()
        Ki_corrected = np.exp(-dg_corrected / (tools.CTE.R*self.absolute_temperature))
        return dg_corrected, Ki_corrected

if __name__ == '__main__':
    import os, glob
    from palmiche.utils.tools import get_sigma
    umbrellas = sorted(glob.glob('/home/ale/mnt/smaug/MD/NEW/docking_min_equi/umbrella_iteration/umbrella_Q'))
    ki_experimental = 66E-9
    for umbrella in umbrellas:
        # try:
        corrector = dG_corrector(
            BoundRegion=2,
            cylinder_radium=0.4,
            cylinder_force_cte_restraint=500,
            orientation_theta_naught=np.pi/4,
            absolute_temperature=303.15,
            PmfXvgPath = os.path.join(umbrella, '7e27/BH246/windows/bsResult_coord0_selected.xvg'),
            integrator='simpson',
        )
        dg_corrected, Ki_corrected = corrector()
        # Experimental - theoretical model
        delta_delta_G = -tools.CTE.R*303.15*(np.log(ki_experimental/Ki_corrected))
        print(f"{os.path.basename(umbrella):>50}: dg_corrected = {dg_corrected:.2f}, Ki_corrected = {Ki_corrected}, delta_delta_G = {delta_delta_G}")

        # except Exception as e:
        #     print(e)



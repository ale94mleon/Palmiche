#!/usr/bin/env python3
# -*- coding: utf-8 -*-

from scipy import integrate
import numpy as np
from palmiche.utils import xvg, tools
import matplotlib.pyplot as plt


def main(BoundRegion, radium, sigma, absolute_temperature, PmfXvgPath, integrator = "simpson"):
    """Get the IC50 from a umbrella simulation.
    The method is based in impose equal probabilities to be bound and unbound states.
    For that it is necessary the definition of the Bound region (which is done with BoundRegion)
    In the general case we have PMF = f(x) and we have the following interval:
    x0-----x1-----x2-----x3-----L----->inffinite
    Where x0 is the ligand inside the protein and inffinite is the ligand in the bulk.
    If we deffine the bound region in the interval [x1;x2]. Then the probability to be bound will be;
    Pbound = P_x1_x2 = I_x1_x2(exp(-E(x)/KbT)dx)/Z
    Where Z is the partion function of the system
    With the same analysis we could get the probability to be in all the rest of the {X}, that is suppouse that the system is unbounded. Therefore:
    Punbound = P_x0_x1 + P_x2_x3 + P_x3_xL
    where:
    P_x0_x1 = I_x0_x1(exp(-E(x)/KbT)dx)/Z
    P_x2_x3 = I_x2_x3(exp(-E(x)/KbT)dx)/Z
    P_x3_L = I_x3_L(exp(-E(x)/KbT)dx)/Z

    Generally in around x3 the PMF shows a flatt region, so we assume that beyond x3 the energy of the system will remain constant and equal to E(x3).
    For numerical reason we set thi value as the average of the last five points.
    if that is so:
    P_x3_L = I_x3_L(exp(-E(x)/KbT)dx)/Z = I_x3_L(exp(-E(x3)/KbT)dx)/Z = exp(-E(x3))(L-x3)/Z
    Finally imposed the equality:
    Pbound = Punbound
    P_x1_x2 = P_x0_x1 + P_x2_x3 + P_x3_L
    I_x1_x2(exp(-E(x)/KbT)dx)/Z = I_x0_x1(exp(-E(x)/KbT)dx)/Z + I_x2_x3(exp(-E(x)/KbT)dx)/Z + exp(-E(x3))(L-x3)/Z
    I_x1_x2(exp(-E(x)/KbT)dx) = I_x0_x1(exp(-E(x)/KbT)dx) + I_x2_x3(exp(-E(x)/KbT)dx) + exp(-E(x3))(L-x3)
    exp(-E(x3))(L-x3) = I_x1_x2(exp(-E(x)/KbT)dx) - I_x0_x1(exp(-E(x)/KbT)dx) - I_x2_x3(exp(-E(x)/KbT)dx)

    L = [I_x1_x2(exp(-E(x)/KbT)dx) - I_x0_x1(exp(-E(x)/KbT)dx) - I_x2_x3(exp(-E(x)/KbT)dx)] / exp(-E(x3)) + x3

    Finally, to get the IC50 because the cylinder restrain used on the umbrella simulation:
    V = np.pi*(radium + 2*sigma)**2*L*1E-24
    1E-24 is to convert from nm to dm. This implies that all the lenght units must be in [nm], Energy in [kJ/mol] and temperature in [K].

    IC50 = 1 / (Nav*V) [mol/L]
    There are several possibilities that are handled on the code:
        It is not necessary to extrapolate the data. If not, only the simulation data will be used
        The minimum definition could be from the first point to some distance or in between.

    Args:
        BoundRegion (float, or list of two floats): The definition of the bounded region
        radium (float): The radium used for the flat bottomed potential in [nm]
        sigma (float): The sigma of the flat bottomed potential (you could get it using the function tools.get_sigma()) in [nm]
        absolute_temperature (float): The simulation temperature in [K]
        PmfXvgPath (string path-like): The path to the xvg file with the PMF data (reaction coordinate[nm], energy [kJ/mol])
        integrator (str, optional): Type of integrator to be used. Defaults to "simpson". Valid methods are: 'simpson', 'trapezoid', 'trapz'

    Raises:
        ValueError: [description]

    Returns:
        [type]: [description]
    """
    # Selection of the integration method
    if integrator == "simpson":
        integrator = integrate.simpson
    elif integrator == "trapezoid":
        integrator = integrate.trapezoid
    elif integrator == "trapz":
        integrator = np.trapz
    else:
        raise ValueError(f"{integrator} is not a valid integration method. Valid methods are: 'simpson', 'trapezoid', 'trapz'")

    Nav = 6.02214076E23
    data = xvg.XVG(PmfXvgPath).data
    data = data[data[:,0] < 4.5]
    data = data[data[:,0] > 0.00]
    data = np.array(sorted(data,key=lambda x:x[0])) # I am not able to use, the sort method of a numpy array in this case
    xmax = np.max(data[:,0])
    print(data[:,1].max()-data[:,1].min())

    try:
        x1 = BoundRegion[0]
        x2 = BoundRegion[1]
    except:
        x1 = None
        x2 = BoundRegion

    if x1:
        data_x0_x1 = data[(data[:,0] <= x1)]
        data_x1_x2 = data[(x1 < data[:,0]) and (data[:,0] <= x2)]
        data_x2_x3 = data[(data[:,0] > x2)]
        raise ValueError('Not coded yet')
        # Esto no lo he porgramado todavia
        # Ahora tengo dos intervlos el minimo esta en el medio
    else:
        # The minimum is deffined from the first point of the data to X2
        data_x0_x2 = data[(data[:,0] <= x2)]
        data_x2_x3 = data[(data[:,0] > x2)]


        # This is just to trnaslate the minimum and see the impact
        #to_add = np.array(len(data_x0_x2)*[(0,30)])
        #data_x0_x2 += to_add
        new_data = np.concatenate((data_x0_x2, data_x2_x3), axis = 0)
        fig, ax = plt.subplots(3,1,figsize = (16,9))
        ax[2].plot(data[:,0],data[:,1])
        ax[0].plot(new_data[:,0],new_data[:,1])
        ax[1].plot(new_data[:,0],tools.BoltzmannFactor(new_data[:,1], absolute_temperature))
        fig.savefig("pmf_boltzman_factor.pdf")

        # Take the last five energies to average, and get the Boltzmann Factor of the last point
        BoltzmannFactor_x3 = tools.BoltzmannFactor(np.mean(data[-5:,1]),absolute_temperature)

        P_x0_x2 = integrator(tools.BoltzmannFactor(data_x0_x2[:,1],absolute_temperature),data_x0_x2[:,0])
        P_x2_x3 = integrator(tools.BoltzmannFactor(data_x2_x3[:,1],absolute_temperature),data_x2_x3[:,0])
        #El problema aqui es que tengo que jugar con la data original, y cuando se me acaben los puntos entonces empezar a interpolar con un dX y empezar a adicionar puntos a la data, hasta que converja la tolerancia o hasta que se me acaben los puntos
        #Probablemente sea mejor correr sobre la data, probar una vez si con el maximo se obtiene convergencia (el signo sera imporrtante) sino ya ir directo a la suma de puntos
        if P_x2_x3 > P_x0_x2:
            L = data_x2_x3[1,0]
            error = np.abs(integrator(tools.BoltzmannFactor(data_x2_x3[:1,1],absolute_temperature),data_x2_x3[:1,0]) - P_x0_x2)
            for i in range(2, len(data_x2_x3)):
                tmp_error = np.abs(integrator(tools.BoltzmannFactor(data_x2_x3[:i,1],absolute_temperature),data_x2_x3[:i,0]) - P_x0_x2)
                if tmp_error < error:
                    error = tmp_error
                    L = data_x2_x3[i,0]

        else:
            # n = N/Nav
            # C = n/V = N/(Nav*V)
            L = xmax + (P_x0_x2 - P_x2_x3) / BoltzmannFactor_x3

        V = np.pi*(radium + 2*sigma)**2*L*1E-24
        IC50 = 1 / (Nav*V)
        deltaG = max(data[:,1]) - min(data[:,1])
        # print("Error: ",np.exp(-6/(8.31E-3*303.15)))
        # print(f"Raw delatG = {deltaG}; Raw IC50 = {np.exp(-deltaG/(8.31E-3*303.15))}, After Correction deltaG = {-8.31E-3*303.15*np.log(IC50)}")
        with open('ic50.txt', "w") as f:
            f.write(str(IC50))
    return IC50






if __name__ == '__main__':
    import os, glob
    from palmiche.utils.tools import get_sigma
    root = '/Users/klimt/mnt/HHDs/Backup1/MD/NEW/docking_min_equi/umbrella_iteration'
    umbrellas = sorted(glob.glob(f'{root}/umbrella_Q*'))

    for umbrella in umbrellas:
        try:
            ic50 = main(2, 0.4, 2*get_sigma(303.15,500), 303.15, os.path.join(umbrella, '7e27/cis1_bh267_m_vacumm_Cell_682_SP_param/windows/bsResult_coord0_selected.xvg'))#cis1_bh267_m_vacumm_Cell_682_SP_param #sys_MMV007839_Cell_891_SP_param
            print(f"{os.path.basename(umbrella):>50}: IC50 = {ic50}")

        except Exception as e:
            print(e)


            pass
    # with open('ic50.txt', "w") as f:
    #     f.write(str(ic50))

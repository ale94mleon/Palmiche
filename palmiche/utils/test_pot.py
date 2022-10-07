from palmiche.utils import xvg, potentials
import matplotlib.pyplot as plt
import numpy as np
orient_restrain_init = 30
pullf = ("pull_pullf.xvg")
pullx = ("pull_pullx.xvg")
# Here I will get the info of the COM and the orientation restrain
t = xvg.XVG(pullf).data[:,0]
plt.subplot(2, 1, 1)
plt.ylabel("V [kJ/mol]")
for coord in range(11, xvg.XVG(pullf).data.shape[1]):
    f = np.pi/180*xvg.XVG(pullf).data[:,coord]
    x = xvg.XVG(pullx).data[:,coord]    
    pot = potentials.POT("flat_bottom", f, x, ref = orient_restrain_init).pot
    plt.plot(t, pot, label = str(coord))
    #plt.ylim(0,2000)
plt.legend(loc="upper right")

plt.subplot(2, 1, 2)
plt.xlabel("Time")
plt.ylabel("COM [deg]")
for coord in range(11, xvg.XVG(pullf).data.shape[1]):
    x = xvg.XVG(pullx).data[:,coord]    
    plt.plot(t, x, label = str(coord))
    #plt.ylim(0,60)
plt.legend(loc="upper right")   
plt.show()

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import simpson
import os

def plot_ks(wave, a_eff, C_abs_sim, C_sca_sim, tau_NH,A,beta,rho_ratio,m_H):
    
    """
FUNCIÓN QUE CALCULA KAPPA
    """
    
    k_ext =tau_NH/(rho_ratio*m_H)

    x=np.log(a_eff*1e-4)

    # Creamos array 2D: filas = radios, columnas = longitudes de onda
    C_abs_array = np.array([C_abs_sim[a] for a in a_eff])
  
    C_sca_array = np.array([C_sca_sim[a] for a in a_eff])
    

    factor = np.exp(x * (1 - beta))[:, np.newaxis]  

    # tau/N_H por radio
    integrando_abs = C_abs_array*factor
    
    integrando_sca = C_sca_array*factor

    tau_NH_abs = A*simpson(integrando_abs, x,axis=0)
    
    tau_NH_sca= A*simpson(integrando_sca, x,axis=0)
    
    k_abs = tau_NH_abs/(rho_ratio*m_H)
    
    k_sca = tau_NH_sca/(rho_ratio*m_H)
    
    plt.figure()
    plt.plot(wave,k_ext,color='Blue',label='Ext',zorder=3)
    plt.plot(wave,k_abs,color='Orange',label='Abs',zorder=2)
    plt.plot(wave,k_sca,color='Green',label='Sca')
    plt.xlabel(r"$\lambda\;(\mu m)$")
    plt.ylabel(r"$\kappa\;(cm^{2}/g)$")
    plt.yscale("log")
    plt.xlim([2, 30])
    plt.ylim([1,1e4])
    plt.legend()
    plt.title(r'Contribución de las distintas $\kappa$')



    # nombre_archivo = f"kappa.pdf"
    # ruta_graficos = "/home/atcsol/Escritorio/graficos_pdf/main"
    # os.makedirs(ruta_graficos, exist_ok=True)
    # ruta_completa = os.path.join(ruta_graficos, nombre_archivo)
    # plt.savefig(ruta_completa)
    # plt.close()

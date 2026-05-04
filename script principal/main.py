import numpy as np
import matplotlib.pyplot as plt
import os
import main_funcs as mf
from C_sim import calc_C_sim
from ks import plot_ks

# ============================================================
# CARPETA DONDE GUARDAREMOS LOS RESULTADOS
# ============================================================

ruta_graficos = "/home/atcsol/Escritorio/graficos_pdf/main"
os.makedirs(ruta_graficos, exist_ok=True)

# ============================================================
# CARGAMOS DATOS DE LA PRIMERA GEOMETRÍA
# ============================================================

ruta_qtable = "/home/atcsol/Escritorio/TFG/gfortran/Ejemplos/examples_exp/SILICATOS_Y_AGUA/qtable"
data = np.loadtxt(ruta_qtable, skiprows=15)

# ============================================================
# EXTRAEMOS DATOS
# ============================================================

a_eff = np.unique(data[:, 0])

wave = data[:, 1]

Q_abs = {}

Q_sca = {}

for a in a_eff:
    
    mask = np.isclose(data[:, 0], a)
    
    Q_abs[a] = data[mask, 3]
    
    Q_sca[a] = data[mask, 4]

wave_plt = wave[mask]

# ============================================================
# CALCULAMOS LAS SECCIONES EFICACES DE ABS Y DISP
# ============================================================

C_abs_sim = {a: calc_C_sim({a: Q_abs[a]}, np.array([a]))[0] for a in a_eff}

C_sca_sim = {a: calc_C_sim({a: Q_sca[a]}, np.array([a]))[0] for a in a_eff}

# ============================================================
# CALCULAMOS TAU/N_H Y KAPPA
# ============================================================

flag=1

tau_NH,beta,A,rho_ratio,_,m_H=mf.calculo_tau_NH(a_eff,C_abs_sim,C_sca_sim,flag)

plt.plot(wave_plt,tau_NH)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau/N_{H}\;(cm^{2}/H)$")
plt.yscale("log")
plt.xlim([2, 30])
plt.title('Opacidad (Silicatos+Agua) por átomo de Hidrógeno')
plt.savefig(os.path.join(ruta_graficos, "tau_NH silicatos.pdf"))
plt.close()

plot_ks(wave_plt, a_eff, C_abs_sim, C_sca_sim, tau_NH,A,beta,rho_ratio,m_H)

# ============================================================
# CARGAMOS DATOS DE LA PRIMERA GEOMETRÍA
# ============================================================

ruta_qtable2 = "/home/atcsol/Escritorio/TFG/gfortran/Ejemplos/examples_exp/CARBONO_Y_AGUA/qtable"
data_2 = np.loadtxt(ruta_qtable2, skiprows=15)

# ============================================================
# EXTRAEMOS DATOS
# ============================================================


# Diccionarios por a_eff
Q_abs_2 = {}

Q_sca_2 = {}

for a in a_eff:
    
    mask = np.isclose(data_2[:, 0], a)
    
    Q_abs_2[a] = data_2[mask, 3]
    
    Q_sca_2[a] = data_2[mask, 4]

# ============================================================
# CALCULAMOS LAS SECCIONES EFICACES DE ABS Y DISP
# ============================================================

C_abs_sim_2 = {a: calc_C_sim({a: Q_abs_2[a]}, np.array([a]))[0] for a in a_eff}

C_sca_sim_2 = {a: calc_C_sim({a: Q_sca_2[a]}, np.array([a]))[0] for a in a_eff}

# ============================================================
# CALCULAMOS TAU/N_H Y KAPPA
# ============================================================

flag_2=2

tau_NH_2,_,A_2,_,rho_ratio_2,_=mf.calculo_tau_NH(a_eff,C_abs_sim_2,C_sca_sim_2,flag_2)

plt.plot(wave_plt,tau_NH_2)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau/N_{H}\;(cm^{2}/H)$")
plt.yscale("log")
plt.xlim([2, 30])
plt.title('Opacidad (Carbono+Agua) por átomo de Hidrógeno')
plt.savefig(os.path.join(ruta_graficos, "tau_NH carbono.pdf"))
plt.close()

tau_NH_t=tau_NH+tau_NH_2

plt.plot(wave_plt,tau_NH_t)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau/N_{H}\;(cm^{2}/H)$")
plt.yscale("log")
plt.xlim([2, 30])
plt.title('Opacidad (Silicatos+Grafito+Agua) por átomo de Hidrógeno')
plt.savefig(os.path.join(ruta_graficos, "tau_NH carbono y silicatos.pdf"))
plt.close()

plot_ks(wave_plt, a_eff, C_abs_sim_2, C_sca_sim_2, tau_NH_2,A_2,beta,rho_ratio_2,m_H)

# ============================================================
# CARGAMOS DATOS DE NGC 3256-S 
# ============================================================

ruta_galaxia = "/home/atcsol/Escritorio/TFG/Espectros galaxias/ngc3256_S_2c.txt"
data_galaxia = np.loadtxt(ruta_galaxia, comments=';')

# ============================================================
# EXTRAEMOS DATOS Y GRAFICAMOS
# ============================================================

wave_galaxia=data_galaxia[:,0]

flujo_galaxia=data_galaxia[:,1]

plt.figure()
plt.plot(wave_galaxia,flujo_galaxia)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$f_{\nu}\;(mJy)$")
plt.yscale("log")
plt.title('NGC 3256-S')
plt.xlim([2,28.63])
plt.savefig(os.path.join(ruta_graficos, "flujo galaxia NGC.pdf"))
plt.close()

# ============================================================
# ELIMINAMOS EL CONTINUO
# ============================================================

wave_ajuste= [2,2.5, 4.5, 5.5, 8.0, 12.7,15,18,20,25,27.5]

flujo_spline,flujo_sel,wave_sel=mf.calculo_spline(wave_ajuste,wave_galaxia,flujo_galaxia)

plt.figure()
plt.plot(wave_galaxia,flujo_galaxia,label='Espectro',zorder=1)
plt.plot(wave_galaxia,flujo_spline,label='Spline',zorder=2)
plt.scatter(wave_sel,flujo_sel,s=10,color='black',label='Puntos elegidos',zorder=3)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$f_{\nu}\;(mJy)$")
plt.yscale("log")
plt.title('Ajuste spline')
plt.legend()
plt.xlim([2, 28.63])
plt.savefig(os.path.join(ruta_graficos, "ajuste spline.pdf"))
plt.close()

tau=-np.log(flujo_galaxia/flujo_spline)

plt.figure()
plt.plot(wave_galaxia,tau)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau$")
plt.xlim([2,28.63])
plt.title('Profundidad óptica en NGC 3256-S')
plt.savefig(os.path.join(ruta_graficos, "profundidad óptica NGC.pdf"))
plt.close()

# ============================================================
# LIMPIAMOS LA LÍNEAS DE EMISIÓN
# ============================================================

tau_masked=mf.limpiar_continuo(wave_galaxia, tau)

plt.figure()
plt.plot(wave_galaxia,tau_masked)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau$")
plt.xlim([2,28.63])
plt.title('Profundidad óptica en NGC 3256-S')
plt.savefig(os.path.join(ruta_graficos, "profundidad óptica NGC limpia.pdf"))
plt.close()

# ============================================================
# OBTENEMOS N_H
# ============================================================

tau_T,NH,std_err_NH,r_value_agua,r_value_sil=mf.calculo_NH(wave_plt,tau_NH_t,wave_galaxia,
                                                           tau_masked)

plt.figure()
plt.plot(wave_galaxia,tau_masked,label=r'$\tau$ galaxia')
plt.plot(wave_plt,tau_T,label=r'$\tau$ simulación')
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau$")
plt.xlim([2, 28.63])
plt.legend()
plt.title('Profundidad óptica en NGC 3256-S')
plt.savefig(os.path.join(ruta_graficos, "tau galaxia vs tau simulación.pdf"))
plt.close()

# ============================================================
# OBTENEMOS Y CORREGIMOS EL REDSHIFT
# ============================================================

A_fit, lambda0_fit, sigma_fit, C_fit=mf.calculo_gaussiana(wave_galaxia,flujo_galaxia)

lambda_rest=9.665

z = lambda0_fit / 9.665 - 1

wave_galaxia_corr=wave_galaxia/(1+z)

wave_ajuste_2= [2,2.7,3.5,5.3,7.2,14.5,24,27.5,30]

flujo_spline_2,flujo_sel_2,wave_sel_2=mf.calculo_spline(wave_ajuste_2,wave_plt,tau_T)

plt.figure()
plt.plot(wave_plt,tau_T,label='Espectro',zorder=1)
plt.plot(wave_plt,flujo_spline_2,label='Spline',zorder=2)
plt.scatter(wave_sel_2,flujo_sel_2,s=10,color='black',label='Puntos elegidos',zorder=3)
plt.xlabel(r"$\lambda\;(\mu m)$")
plt.ylabel(r"$\tau$")
plt.title('Ajuste spline')
plt.legend()
plt.xlim([2, 30])
plt.savefig(os.path.join(ruta_graficos, "ajuste spline simulación.pdf"))
plt.close()

# ============================================================
# CSUSTRAEMOS EL CONTINUO A LA SIMULACIÓN Y GRAFICAMOS RESULTADO
# ============================================================

tau_T_2=tau_T-flujo_spline_2

plt.figure()
plt.plot(wave_galaxia_corr,tau_masked,label=r'$\tau$ galaxia')
plt.plot(wave_plt,tau_T_2,label=r'$\tau$ simulación')
plt.xlabel(r"$\lambda_{rest}\;(\mu m)$")
plt.ylabel(r"$\tau$")
plt.xlim([2,28.63])
plt.legend()
plt.title('Profundidad óptica en NGC 3256-S con redshift corregido')
plt.savefig(os.path.join(ruta_graficos, "zoom silicatos.pdf"))
# plt.savefig(os.path.join(ruta_graficos, "tau galaxia vs tau simulación (z corregido).pdf"))
plt.close()

import numpy as np
from scipy.integrate import simpson
from scipy.interpolate import interp1d
from scipy.stats import linregress
from scipy.optimize import curve_fit
from scipy.interpolate import CubicSpline

def calculo_tau_NH(a_eff,C_abs_sim,C_sca_sim,flag): 
    
    """
FUNCIÓN PARA CALCULAR LA OPACIADA DE UNA GEOMETRÍA
    """
    
    beta = 3.5
    
    rho_solid = 2.5  # g/cm^3 #Jones 2013 (SILICATOS)

    rho_ratio = 5.8e-3
    
    rho_solid_2 = 1.3  # g/cm^3 #Jones 2013 (CARBONO)

    rho_ratio_2 = 0.6e-3

    rho_solid_agua = 0.918 # g/cm^3 #Jones 2013
    
    m_H=1.6724e-24;

    a_max_cm = a_eff[-1] * 1e-4

    a_min_cm = a_eff[0] * 1e-4
    
    x=np.log(a_eff*1e-4)
    
    if flag==1:
        
        C_ext_sim = {a: C_abs_sim[a] + C_sca_sim[a] for a in a_eff}
        
        A = (3*(4-beta))*m_H / (4*np.pi*(0.79*rho_solid+0.21*rho_solid_agua)*(a_max_cm**(4-beta) - a_min_cm**(4-beta))) * (rho_ratio)
        
        C_ext_array = np.array([C_ext_sim[a] for a in a_eff])

        factor = np.exp(x * (1 - beta))[:, np.newaxis]  

        integrando = C_ext_array*factor

        tau_NH = A*simpson(integrando, x,axis=0)
        
    if flag==2:
        
        C_ext_sim_2 = {a: C_abs_sim[a] + C_sca_sim[a] for a in a_eff}
        
        A = (3*(4-beta))*m_H / (4*np.pi*(0.63*rho_solid_2+0.37*rho_solid_agua)*(a_max_cm**(4-beta) - a_min_cm**(4-beta))) * (rho_ratio_2)

        C_ext_array_2 = np.array([C_ext_sim_2[a] for a in a_eff])
        
        factor = np.exp(x * (1 - beta))[:, np.newaxis]  

        integrando_2 = C_ext_array_2*factor

        tau_NH = A*simpson(integrando_2, x,axis=0)
        
    return(tau_NH,beta,A,rho_ratio,rho_ratio_2,m_H)

def calculo_spline(wave_ajuste,wave,flujo):
    
    """
FUNCIÓN PARA OBTENER EL SPLINE PARA SUSTRAER EL CONTINUO
    """
    
    indices = [np.argmin(np.abs(wave- l)) for l in wave_ajuste]

    wave_sel = wave[indices]

    flujo_sel = flujo[indices]

    spline = CubicSpline(wave_sel, flujo_sel)

    flujo_spline = spline(wave)
    
    return(flujo_spline,flujo_sel,wave_sel)

def limpiar_continuo(wave_galaxia, tau):
    
    """
FUNCIÓN PARA ELIMINAR LAS LÍNEAS DE EMISIÓN
    """
    
    delta_v = 236.9426752
    
    c = 300000.0
    
    lambda_0_list = [2.185, 2.65, 7.05, 9.75, 12.93, 15.7, 17.2, 18.885]
   
    N = 3
    
    tau_masked = tau.copy()
    
    for lam0 in lambda_0_list:
        
        delta_lambda = (delta_v / c) * lam0
        
        low = lam0 - N * delta_lambda
        
        high = lam0 + N * delta_lambda
        
        region = (wave_galaxia >= low) & (wave_galaxia <= high)
        
        tau_masked[region] = np.nan
        
    lambda_0_list_2 = [9.075,10.61,12.39]
   
    N_2 = 1
    
    for lam0 in lambda_0_list_2:
        
        delta_lambda = (delta_v / c) * lam0
        
        low = lam0 - N_2 * delta_lambda
        
        high = lam0 + N_2 * delta_lambda
        
        region = (wave_galaxia >= low) & (wave_galaxia <= high)
        
        tau_masked[region] = np.nan
        
    
    return (tau_masked)

def calculo_NH(wave_plt,tau_NH_t,wave_galaxia,tau_masked):
    
    """
FUNCIÓN PARA OBTENER N_H Y TAU
    """
    
    interp = interp1d(wave_plt, tau_NH_t, bounds_error=False, fill_value="extrapolate")

    tau_NH_T = interp(wave_galaxia)

    mask_fit_agua = (~np.isnan(tau_masked)) & (wave_galaxia >= 2.75) & (wave_galaxia <= 4)

    NH_agua, _, r_value_agua, _, std_err_agua = linregress(tau_NH_T[mask_fit_agua],tau_masked[mask_fit_agua])

    mask_fit_sil = (~np.isnan(tau_masked)) & (wave_galaxia >= 7.5) & (wave_galaxia <= 12.7)

    NH_sil, _, r_value_sil, _, std_err_sil = linregress(tau_NH_T[mask_fit_sil],tau_masked[mask_fit_sil])

    NH=(NH_agua+NH_sil)/2

    std_err_NH=(std_err_agua+std_err_sil)/2

    tau_T=tau_NH_t*NH
    
    return(tau_T,NH,std_err_NH,r_value_agua,r_value_sil)

def gaussiana(x, A, x0, sigma, C):
    
    """
FUNCIÓN QUE CALCULA UNA GAUSSIANA
    """
    
    return (A * np.exp(-(x - x0)**2 / (2 * sigma**2)) + C)

def calculo_gaussiana(wave_galaxia,flujo_galaxia):
    
    """
FUNCIÓN QUE CALCULA LA GAUSSIANA PARA CORREGIR EL RESHIFT
    """
    
    mask_3=(wave_galaxia>=9.7)&(wave_galaxia<=10.4)

    wave_gauss=wave_galaxia[mask_3]

    flujo_gauss=flujo_galaxia[mask_3]
    
    A0 = max(flujo_gauss) - np.median(flujo_gauss)
    
    x0_0 = wave_gauss[np.argmax(flujo_gauss)]
    
    sigma0 = 0.05                                  
    
    C0 = np.median(flujo_gauss)                    

    p0 = [A0, x0_0, sigma0, C0]

    popt, pcov = curve_fit(gaussiana, wave_gauss, flujo_gauss, p0=p0)

    A_fit, lambda0_fit, sigma_fit, C_fit = popt
    
    return(A_fit, lambda0_fit, sigma_fit, C_fit)

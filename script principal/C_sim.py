import numpy as np

def calc_C_sim(Q_dict, radios):
    
    """
FUNCIÓN QUE CALCULA LA SECCIÓN EFICAZ DE ABSORCIÓN Y DISPERSIÓN
    """
    
    C_sim = []

    for a in radios:
        
        C_sim.append(Q_dict[a]*np.pi*(a*1e-4)**2)

    return np.array(C_sim)

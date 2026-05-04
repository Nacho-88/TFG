import numpy as np
import random

"""
SCRIPT PARA ELEGIR DE MANERA ALEATORIA LOS DIPOLOS A SUSTITUIR POR 
EL SEGUNDO MATERIAL 
    """

# ============================================================
# DEFINIMOS LOS PARÁMETROS PRINCIPALES
# ============================================================

input_file="target.out"

output_file="target_modificado.out"

radio_esfera=8

radio_inclusion=4

num_inclusiones=3

material_nuevo=(2, 2, 2)  
 
material_original=(1, 1, 1)

# ============================================================
# EXTRAEMOS DATOS
# ============================================================

with open(input_file, "r") as f:
    
    lines = f.readlines()

header=[]

data_lines=[]

for line in lines:
    
    parts=line.split()
    
    if len(parts)==7 and parts[0].isdigit():
        
        data_lines.append(parts)
        
    else:
        
        header.append(line)

dipolos=[]

for row in data_lines:
    
    JA=int(row[0])
    IX=int(row[1])
    IY=int(row[2])
    IZ=int(row[3])
    ICX=int(row[4])
    ICY=int(row[5])
    ICZ=int(row[6])
    
    dipolos.append([JA,IX,IY,IZ,ICX,ICY,ICZ])

# ============================================================
# CALCULAMOS EL CENTRO DE LA ESFERA
# ============================================================

coords=np.array([[d[1],d[2],d[3]] for d in dipolos])

centro=coords.mean(axis=0)

# ============================================================
# SELECCIONAMOS LOS DIPOLOS VÁLIDOS PARA CAMBIAR
# ============================================================

validos=[]

for d in dipolos:
    
    pos=np.array([d[1],d[2],d[3]])
    
    if np.linalg.norm(pos-centro)>=4:
        
        validos.append(d)

# ============================================================
# CAMBIAMOS DIPOLOS SIN SOLAPAR INCLUSIONES
# ============================================================

for _ in range(num_inclusiones):

    while True:
        
        centro_sel=random.choice(validos)
        
        if (centro_sel[4],centro_sel[5],centro_sel[6])!=material_nuevo:
            
            break

    cx,cy,cz=centro_sel[1],centro_sel[2],centro_sel[3]

    for d in dipolos:
        
        x,y,z=d[1],d[2],d[3]
        
        dist=np.sqrt((x-cx)**2 +(y-cy)**2+(z-cz)**2)
        
        if dist<=radio_inclusion:
            
            d[4]=material_nuevo[0]
            
            d[5]=material_nuevo[1]
            
            d[6]=material_nuevo[2]

# ============================================================
# GUARDAMOS RESULTADOS CON EL MISMO FORMATO QUE EL INPUT
# ============================================================

with open(output_file, "w") as f:
    
    for line in header:
        
        f.write(line)
        
    for d in dipolos:
        
        f.write(f"{d[0]:5d}{d[1]:5d}{d[2]:5d}{d[3]:5d}{d[4]:3d}{d[5]:3d}{d[6]:3d}\n")

print("Archivo generado:", output_file)
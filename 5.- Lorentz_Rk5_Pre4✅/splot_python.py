import numpy as np
import matplotlib.pyplot as plt

print("Cargando datos...")
t, x, y, z = np.loadtxt('datos.dat', unpack=True)

# 2. CONFIGURACIÓN DE LA GRÁFICA 3D
print("Generando gráfica...")
fig = plt.figure(figsize=(10, 8))
ax = fig.add_subplot(111, projection='3d')

# 3. GRAFICADO (Líneas y puntos pequeños)
# linestyle='-' (Línea continua)
# marker='o' (Puntos)
# markersize=0.5 (Puntos muy pequeños)
# linewidth=0.3 (Línea fina para que los 1,000,000 puntos no saturen la imagen)
ax.plot(x, y, z, color='blue', linestyle='-', marker='o', 
        markersize=0.5, linewidth=0.3, alpha=0.8)

# 4. ESTÉTICA
ax.set_title('Atractor de Lorenz', fontsize=14, fontweight='bold')
ax.set_xlabel('Eje X (x)')
ax.set_ylabel('Eje Y (y)')
ax.set_zlabel('Eje Z (z)')

# Cambiamos el color de fondo para que resalte más
ax.set_facecolor('white')

# 5. GUARDADO VECTORIAL Y VISUALIZACIÓN
# Guardamos en .pdf o .svg (Formatos vectoriales puros, ideal para artículos/tesis)
plt.savefig('atractor_lorenz.png', format='png', dpi=300, bbox_inches='tight')
print("Gráfica guardada como 'atractor_lorenz.pdf'")

# Mostramos la ventana interactiva
plt.show()

import numpy as np
import matplotlib.pyplot as plt

def LU(n, a, e, c, b):
  w = np.zeros(n)
  u = np.zeros(n)
  y = np.zeros(n)
  x = np.zeros(n)

  w[0] = a[0]
  for i in range(n-1):
    u[i] = e[i]/w[i]
    w[i+1] = a[i+1] - c[i+1]*u[i]

  y[0] = b[0]/w[0]

  for i in range(1,n):
    y[i] = (b[i] - c[i]*y[i-1])/w[i]

  x[n-1] = y[n-1]

  for i in range(n-2,-1,-1):
    x[i] = y[i] - u[i]*x[i+1]

  return x

#Definición de Constantes

L = 3 # Longitud de la barra
alpha = 0.05
h = 1e-2
k = 1e-2
N = int(L / h) - 1
M = 5000
gamma = (alpha*k)/h**2
# Condiciones iniciales

x = np.linspace(0, L, N+2)
u = np.zeros((N+2,M))
u[:,0] =  10 * np.exp(-2 * (x - 0.5)**2) # Temperatura inicial

# Condiciones de frontera
u[0,:] = 10
u[N+1,:] = 5

a = np.zeros(N)
e = np.zeros(N)
c = np.zeros(N)
b = np.zeros(N)

for i in range(N):
  a[i] = 2+2*gamma
  e[i] = -gamma
  c[i] = -gamma

#creo que estos es C-N
for j in range(M-1): # Pasos en el tiempo
  for i in range(1,N+1): # Pasos en la posición
    b[i-1] = gamma*(u[i+1,j]+u[i-1,j])+(2-2*gamma)*u[i,j]

  b[0] = b[0] + gamma*u[0,j+1]
  b[N-1] = b[N-1] + gamma*u[N+1,j+1]
  u[1:N+1,j+1] = LU(N, a, e, c, b)


# Graficación

v = np.zeros((M,N+1))
for i in range(N+1):
  for j in range(M):
      v[j,i]=u[i,j]

plt.imshow(v,aspect = 'auto',cmap = 'jet',origin = 'lower',extent = [0,L,0,h*M], interpolation= 'bilinear')
plt.colorbar(label = 'Temperatura')
plt.xlabel('Posición (x)',fontsize=22)
plt.ylabel('Tiempo (t)',fontsize=22)
#plt.title('Temperatura en la varilla',fontsize=24)
plt.show()
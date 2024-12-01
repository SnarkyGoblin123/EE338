import numpy as np
Wc=1.022
N=24
pi = np.pi
p = []
for k in range(0,24):
    if (k<N/2):
        pole = (Wc*np.cos(pi/2 + pi/(2*N) + k*pi/N), Wc*np.sin(pi/2 + pi/(2*N) + k*pi/N))
        p.append(pole)
    else:
        x = N-1-k
        pole = (Wc*np.cos(pi/2 + pi/(2*N) + x*pi/N), Wc*np.sin(pi/2 + pi/(2*N) + x*pi/N))
        p.append(pole)
print(len(p))
p = np.round(p, 4)
for k in range(24):
    print('$p_{',k,'}$ = ', p[k][0], ' + j', p[k][1], '\\\\')


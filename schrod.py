from NanoTCAD_ViDES import *
from numpy.linalg import eig

hbar = 1.0545718e-34; m0 = 9.1e-31; m_eff = 1; q = 1.6e-19;

def Schrod1d(m, x, U, mode=1):
    N = len(x)
    dx = x[1]-x[0]
    t = (hbar**2)/(2*m*(dx**2))/q
    H = 2*t*eye(N) + diag(U) - t*diag(ones(N-1), 1) - t*diag(ones(N-1), -1)
    eigenValues, eigenVectors = eig(H)
    eigenVectors = eigenVectors/sqrt(dx)
    idx = np.argsort(eigenValues)
    eigenValues = eigenValues[idx]
    eigenVectors = eigenVectors[:,idx]
    return eigenValues, eigenVectors

# print(trapz(eigenVectors[:,0]**2, x))
# plot(eigenVectors[:,0:5])
# plot(eigenValues)
# show()

m = m0*m_eff
N = 30; 
L = 30e-9;
x = linspace(0, L, N)

U = ones(N)*q
U[int(.25*N):int(.75*N)] = 0
eigenValues, eigenVectors = Schrod1d(m, x, U)

# plot(eigenValues, 'ro')
print(trapz(x, eigenVectors[:,0]**2))
plot(x, eigenVectors[:,0])
show()

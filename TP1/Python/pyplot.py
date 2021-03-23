#Example from https://pythonprogramming.net/loading-file-data-matplotlib-tutorial/
import matplotlib.pyplot as plt
import numpy as np

x, y = np.loadtxt('example-10000.data', delimiter=' ', unpack=True)
plt.plot(x, y, label='Sequential - Tab with 10 000 int')

plt.xlabel('Tab size')
plt.ylabel('Time')
plt.title('malloc() + rempliAleatoire() + testPrime()')
plt.legend()
plt.show()

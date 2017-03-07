#!/usr/bin/python

import matplotlib.pyplot as plt
plt.plot([1, 2,4], [.530, 1.272, 1.185], 'r--', [1, 2, 4],[.54, 1.540, 1.524],'b--', [1,2,4],[.54, 1.547, 1.548],'g--')
plt.title('serial qr vs omp-qr non-square matrices')
plt.ylabel('GF/s red= 1000x500, blue = 10000x500, green = 20000x500')
plt.xlabel('number of threads (1 = serial version)')
plt.xlim([0,5])
plt.show()

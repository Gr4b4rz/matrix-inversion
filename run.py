from subprocess import run
from timeit import timeit
from functools import partial
import matplotlib.pyplot as plt
import numpy as np


scope = range(20, 400, 20)
labels = list(map(str, scope))
# func1 = lambda x: timeit(partial(run, ["build/matrix-inversion", "{}".format(x)]), number=5) / 5
#
# # TODO dir to be
# func2 = lambda x: timeit(partial(run, ["build1/matrix-inversion", "{}".format(x)]), number=5) / 5
#
# func3 = lambda x: timeit(partial(run, ["build2/matrix-inversion", "{}".format(x)]), number=5) / 5

func1 = lambda x: timeit(partial(run, ["mpirun", "-n", "1", "build_mpi/matrix-inversion", "{}".format(x)]), number=5) / 5
func2 = lambda x: timeit(partial(run, ["mpirun", "-n", "2", "build_mpi/matrix-inversion", "{}".format(x)]), number=5) / 5
func3 = lambda x: timeit(partial(run, ["mpirun", "-n", "4", "build_mpi/matrix-inversion", "{}".format(x)]), number=5) / 5

a = list(map(func1, [i for i in scope]))
b = list(map(func2, [i for i in scope]))
c = list(map(func3, [i for i in scope]))

width = np.min(np.diff(scope))/3

x = np.arange(len(scope))  # the label locations
width = 0.18  # the width of the bars

fig, ax = plt.subplots()

rects1 = ax.bar(x - width, a, width, label='OpenMpi with 1 cores')
rects2 = ax.bar(x, b, width, label='OpenMpi with 2 cores')
rects3 = ax.bar(x + width, c, width, label='OpenMpi with 4 cores')

ax.set_ylabel('Time Elapsed (s)')
ax.set_title('Matrix inversion mean time')
ax.set_xticks(x)
ax.set_xlabel('Matrix Dimension')
ax.set_xticklabels(labels)
ax.legend()


fig.tight_layout()

plt.show()

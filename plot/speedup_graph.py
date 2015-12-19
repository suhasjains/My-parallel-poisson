from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
plt.rcParams['font.size'] = 18


time = np.array([ 287.373, 96.2159, 49.6662, 9.99871 ]);
cores = np.array([ 24, 72, 144, 720 ]);
speedup = time[0] / time;
ideal = cores / cores[0];

plt.subplot(121)
plt.xlabel('Number of cores')
plt.ylabel('Time taken in seconds')
plt.title('Time taken vs the number of cores for 72 million cells')
plt.xlim([cores[0],cores[3]])
plt.ylim([time[0],time[3]])
plt.yticks(time)
plt.xticks(cores)
plt.plot(cores, time,'b-o', clip_on=False)



plt.subplot(122)
plt.xlabel('Number of cores')
plt.ylabel('Normalized speedup')
plt.title('Speed up vs the number of cores for 72 million cells')
plt.rcParams.update({'axes.titlesize': 'small'})
plt.rcParams.update({'axes.labelsize': 'small'})
plt.xlim([cores[0],cores[3]])
plt.ylim([ideal[0],ideal[3]])
plt.yticks(speedup)
plt.xticks(cores)
plt.plot(cores, speedup,'b-o', clip_on=False, label='Actual speedup')
plt.plot(cores, ideal, 'r-o', clip_on=False, label='Ideal speedup')



#legend
legend = plt.legend(loc='upper center', shadow=True)

# The frame is matplotlib.patches.Rectangle instance surrounding the legend.
frame = legend.get_frame()
frame.set_facecolor('0.90')

# Set the fontsize
for label in legend.get_texts():
    label.set_fontsize('large')

for label in legend.get_lines():
    label.set_linewidth(1.5)  # the legend line width



plt.rcParams.update({'axes.titlesize': 'large'})
#plt.grid(True)
#plt.savefig("speedup.eps")
plt.suptitle("Parallel poisson solver on SahasraT (Cray XC40)");
plt.show()

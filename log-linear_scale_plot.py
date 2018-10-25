
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

# Numbers from -50 to 50, with 0.1 as step
xdomain = np.arange(-50,50, 0.1)
x=[0,1,2,3,4,5,6,7,8,9,10,11,20,30,40,50]
y=[0,1,2,3,4,5,6,7,8,9,10,11,20,30,40,50]

ax = plt.subplot(111)
ax.plot(x, y)
ax.set_xscale('log')
ax.set_xlim((0,550))
ax.spines['left'].set_visible(False)
ax.yaxis.set_ticks_position('right')
ax.yaxis.set_visible(False)

divider = make_axes_locatable(ax)
axLin = divider.append_axes("left", size=2.0, pad=0, sharey=ax)
axLin.set_xscale('linear')
axLin.set_xlim((-10, 0))
axLin.plot(x,y)
axLin.spines['right'].set_visible(False)
axLin.yaxis.set_ticks_position('left')
plt.setp(axLin.get_xticklabels(), visible=True)

plt.title('Linear left, log right')
plt.show()
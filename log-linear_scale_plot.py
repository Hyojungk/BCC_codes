
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
import numpy as np

x=[0,1,2,3,4,5,6,7,8,9,10,11,20,30,40,50]
y=[0,1,2,3,4,5,6,7,8,9,10,11,20,30,40,50]

axlog = plt.subplot(111)
axlog.plot(x, y)
axlog.set_xscale('log')
axlog.set_xlim((0,550))
axlog.spines['left'].set_visible(False)
axlog.yaxis.set_ticks_position('right')
axlog.yaxis.set_visible(False)

divider = make_axes_locatable(axlog)
axLin = divider.append_axes("left", size=2.0, pad=0, sharey=axlog)
axLin.set_xscale('linear')
axLin.set_xlim((-10, 0))
axLin.plot(x,y)
axLin.spines['right'].set_visible(False)
axLin.yaxis.set_ticks_position('left')
plt.setp(axLin.get_xticklabels(), visible=True)

plt.show()

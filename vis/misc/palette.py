
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.colors import rgb2hex 

fig = plt.figure(figsize=(2,8))
cmap = mpl.cm.jet

N = 31

for i in range(N+1):
	col = cmap(1.*i/N)
	colname = rgb2hex(col).upper()
	plt.fill_between([0,1],[i,i],[i+1,i+1], color=col)
	plt.annotate(text=colname, xy=(0.05,i+.5), color='0.25', ha='left', va='center')


plt.xlim([0,1])
plt.ylim([N+1,0])
plt.xticks([])
plt.yticks([])
plt.tight_layout()

plt.savefig('jet.pdf', dpi=500)



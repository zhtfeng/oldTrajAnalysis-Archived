import matplotlib
import matplotlib.pyplot as plt
import numpy as np


labels = ['methyl','H','Vinyl', 'phenyl', 'cyclopropyl', 'ethyl', 'allyl','Mom', '$\mathit{i}$-Pr','CH2-CPP','cyclobutyl','acyl' ]
# men_means = [-321.7,-312.3,-300.3,-293.4,-289.0,-278.6,-271.1,-267.6,-259.0,-252.3,-247.0,-242.0]
men_means = [-329.9,-320.3,-307.2,-301.5,-296.7,-286.4,-278.9,-275.4,-266.7,-260.3,-254.9,-249.1]
men_means.reverse()


x = np.arange(len(labels))  # the label locations
width = 0.35  # the width of the bars

fig, ax = plt.subplots(dpi=600)
rects1 = ax.bar(x, men_means, width,edgecolor='black',color = 'w')


# Add some text for labels, title and custom x-axis tick labels, etc.
ax.set_ylabel('Hydride Affinity (kcal/mol)',fontsize=12,fontname='Times New Roman')
ax.set_xticks(x)
ax.set_xticklabels(men_means,fontsize=8,fontname='Times New Roman')
ax.set_ylim((-348,-228))
ax.tick_params(axis='y', labelsize=8)
for tick in ax.get_yticklabels():
    tick.set_fontname("Times New Roman")
ax.xaxis.set_ticks_position('top')


def autolabel(rects):
    """Attach a text label above each bar in *rects*, displaying its height."""
    for counter,rect in enumerate(rects):
        height = rect.get_height()
        print(rect)
        ax.annotate('{}'.format(labels[counter]),
                    xy=(rect.get_x() + rect.get_width()/2, height-7),
                    xytext=(0, 1),  # 3 points vertical offset
                    textcoords="offset points",
                    ha='center', va='bottom',fontsize=8,fontname='Times New Roman')


#autolabel(rects1)
fig.tight_layout()
plt.savefig('HA_enthalpy.png')
plt.show()
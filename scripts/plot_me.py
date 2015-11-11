from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.lines as mlines
length = [0.3, 1.0, 10.0]
coverage = [75, 150, 225, 300, 450]

#01
#bk_times = [[0.263, 1.514, 9.884, 33.417, 832.0], [0.722, 3.548, 16.169, 31.022], [5.234, 28.809, 547.086]]
#cl_times = [[0.265, 1.499, 9.849, 34.205, 457.801], [0.692, 3.47, 16.917, 33.248], [5.168, 28.349, 2090.666]]

#001
#bk_times = [[0.478, 3.072, 73.329, 13437.074], [0.924, 6.367, 43.502, 380.746], [21.943, 241.138, 59727.867]]
#cl_times = [[0.481, 3.328, 64.014,4200.011], [0.945, 6.516, 44.341, 392.988], [18.766, 211.106, 11371.128]]

#005
bk_times = [[0.407, 2.71, 17.297, 188.512, 2968.672], [0.686, 4.647, 25.251, 113.334, 2202.791], [5.767, 35.582, 135.259, 1558.178]]
cl_times = [[0.397, 2.673, 16.142, 123.531, 1022.705], [0.694, 4.713, 24.734, 117.308, 1516.212], [5.692, 34.159, 132.162, 966.585]]

MAXX = 500

for i in range(0, 3):
    f, ax = plt.subplots()

    bkl = len(bk_times[i])
    cll = len(bk_times[i])

    bk = np.polyfit(coverage[0:bkl], np.log(bk_times[i]), 1)
    cl = np.polyfit(coverage[0:cll], np.log(cl_times[i]), 1)

    bks = 'exp(' + "{0:.3f}".format(bk[0]) + 'x + ' + "{0:.3f}".format(bk[1]) + ')'
    cls = 'exp(' + "{0:.3f}".format(cl[0]) + 'x + ' + "{0:.3f}".format(cl[1]) + ')'

    bkfn = np.poly1d(bk)
    clfn = np.poly1d(cl)

    x = np.linspace(0, MAXX, 100)

    plt.plot(coverage[0:bkl], bk_times[i], 'or')
    plt.plot(coverage[0:cll], cl_times[i], 'ob')

    plt.plot(x, np.exp(bkfn(x)), '-r')
    plt.plot(x, np.exp(clfn(x)), '-b')

    ax.set_yscale('log')
    ax.yaxis.set_label_text('time [in sec]')
    ax.xaxis.set_label_text('coverage')

    blue_line = mlines.Line2D([],[], color='blue', label='CLEVER\n'+ cls)
    red_line = mlines.Line2D([],[], color='red', label='Bron-Kerbosch\n' + bks)

    plt.axis([0, 500, -100, 1000000])
    plt.legend(handles=[red_line, blue_line], bbox_to_anchor=(0.5, 1))

plt.show()

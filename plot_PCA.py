import glob
import numpy as np

#Very lazy code for drawing HIST plot and PP plot

fig = plt.figure(figsize = (32,50))

for i in range(0,16):
    if i==1:
        Qnull_file = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qnull_PC1.npy"
        Qsvd_file  = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qsvd_PC1.npy"
    else:
        Qnull_file = glob.glob("/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qnull_PC{}*".format(i))[0]
        Qsvd_file  = glob.glob("/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qsvd_PC{}*".format(i))[0]
    Qnull_sort = np.load(Qnull_file)
    Qsvd_sort  = np.load(Qsvd_file)
    
    PC = "PC"+str(i)
    xmin = np.min(Qnull_sort)
    xmax = np.max(Qnull_sort)
    bordercolor = '#333333'
    borderwidth = 3
    axis_font_size = 25
    label_font_size = 20
    legend_font_size = 20
    if i == 15:
        j=i-1
    else:
        j=i
    ax = fig.add_subplot(6,3,i+1)
    for tick in ax.yaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
        tick.label1.set_fontweight('light')

    ax.tick_params(axis='both', which = 'major',
                       length = 10, width = borderwidth, pad=10,
                       labelsize = label_font_size+10,
                       color = bordercolor,
                       labelcolor = bordercolor,
                       bottom = True, top = False, left = True, right = False)
    ax.hist(Qsvd_sort, bins = np.linspace(xmin, xmax, 60), alpha = 0.3, label = "All SNPs")#,normed=True)
    ax.hist(Qnull_sort, bins = np.linspace(xmin, xmax, 60), alpha = 0.3, label = r"$χ^2$")#,normed=True)
    ax.legend()
    ax.legend(prop={'size': label_font_size})
    ax.set_xlabel(r"Q-values",{'size': axis_font_size, 'color': bordercolor},labelpad = 25)
    ax.set_ylabel(r"Density",{'size': axis_font_size, 'color': bordercolor},labelpad = 25)
    ax.set_title("{}".format(PC), fontsize=20)
    ax.set_ylim([-12,1950])
    for side, border in ax.spines.items():
        border.set_linewidth(borderwidth)
        border.set_color(bordercolor)
    
#np.save("sorted_Qsvd_0109_0013.npy", Qsvd_sort)
plt.tight_layout()
plt.savefig("/home/anubhavk/Desktop/Internship_Qstat/results/PCA_corrected_Analysis_hist_usingQ.png")

fig = plt.figure(figsize = (32,50))

for i in range(0,16):
    if i==1:
        Qnull_file = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qnull_PC1.npy"
        Qsvd_file  = "/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qsvd_PC1.npy"
    else:
        Qnull_file = glob.glob("/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qnull_PC{}*".format(i))[0]
        Qsvd_file  = glob.glob("/home/anubhavk/Desktop/Internship_Qstat/codebase/lib/usingQ/Qsvd_PC{}*".format(i))[0]
    Qnull_sort = np.load(Qnull_file)
    Qsvd_sort  = np.load(Qsvd_file)
    PC = "PC"+str(i)
    xmin = np.min(Qnull_sort)
    xmax = np.max(Qnull_sort)
    bordercolor = '#333333'
    borderwidth = 3
    axis_font_size = 25
    label_font_size = 20
    legend_font_size = 20
    ax1 = fig.add_subplot(6,3,i+1)
    for tick in ax1.xaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
        tick.label1.set_fontweight('light')
    for tick in ax1.yaxis.get_major_ticks():
        tick.label1.set_fontsize(15)
        tick.label1.set_fontweight('light')

    ax1.tick_params(axis='both', which = 'major',
                       length = 10, width = borderwidth, pad=10,
                       labelsize = label_font_size+10,
                       color = bordercolor,
                       labelcolor = bordercolor,
                       bottom = True, top = False, left = True, right = False)


    ax1.plot([xmin, xmax], [xmin, xmax], lw = 2, color='black', ls = 'dashed')
    ax1.scatter(Qnull_sort, Qsvd_sort,color="#C10020")
    ax1.set_xlabel(r"expected $\bf{Q}$",{'size': axis_font_size, 'color': bordercolor},labelpad = 25)
    ax1.set_ylabel(r"Observed $\bf{Q}$",{'size': axis_font_size, 'color': bordercolor},labelpad = 25)
    ax1.set_ylim([95,570])
    ax1.set_title("{}".format(PC), fontsize=20)
    for side, border in ax1.spines.items():
        border.set_linewidth(borderwidth)
        border.set_color(bordercolor)
    
plt.tight_layout()
plt.savefig("/home/anubhavk/Desktop/Internship_Qstat/results/PCA_corrected_Analysis_pp_usingQ.png")




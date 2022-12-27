import demes
import demesdraw
import matplotlib.pyplot as plt
import sys

demeFile = sys.argv[1]
demeGraphPDF = sys.argv[2]

# demeFile = "/Users/xudaxing/Documents/deleteriousMutation/001_analysis/007_vmap2_1062_spelt/021_Simulation/003_twoWay_proportion0.1_genetration100/001_deme/two_way_admixture_proportion0.1_genetration100.yaml"
# demeGraphPDF = "/Users/xudaxing/Desktop/temp.pdf"

graph = demes.load(demeFile)
demesdraw.tubes(graph)
# plt.show()
plt.savefig(fname=demeGraphPDF,format="pdf")

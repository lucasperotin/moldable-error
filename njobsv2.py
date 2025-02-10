import sys
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.size']=20
mpl.rcParams['legend.fontsize']=12
mpl.rcParams['legend.framealpha']=0.1
import matplotlib.pyplot as plt
import numpy as np

ns = [100,300,500,750,1000]
print(ns)
set_no = int(sys.argv[1])
model = sys.argv[2]
p=int(7500)
l="1e-7"

data_length = []
b_data_length = []
t_data_length = []
a_data_length = []
if model == "rft":
	txt = np.loadtxt("sample_rooftop/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_rooftop/batch_list_sumv2_"+l+"_"+str(p))
elif model == "amd":
	txt = np.loadtxt("sample_amdahl/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_amdahl/batch_list_sumv2_"+l+"_"+str(p))
elif model == "com":
	txt = np.loadtxt("sample_com/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_com/batch_list_sumv2_"+l+"_"+str(p))
elif model == "mix":
	txt = np.loadtxt("sample_mixst/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_mixst/batch_list_sumv2_"+l+"_"+str(p))
elif model == "mixlc":
	txt = np.loadtxt("sample_mixlc/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_mixlc/batch_list_sumv2_"+l+"_"+str(p))
elif model == "pow":
	txt = np.loadtxt("sample_pow/njobs_sumv2_"+l+"_"+str(p))
	txtb = np.loadtxt("sample_pow/batch_list_sumv2_"+l+"_"+str(p))
else:
	print("Modele inconnu")
	exit()

for n in ns:
	if n == 100:
		offset = 0
	if n == 300:
		offset = 30
	if n == 500:
		offset = 60
	if n == 750:
		offset = 90
	if n == 1000:
		offset = 120
	if set_no == 0:
		sum_length = 0
		b_sum_length = 0
		a_sum_length = 0
		t_sum_length = 0
		for i in range(30):
			sum_length += txt[offset+i,1]
			b_sum_length += txtb[offset+i]
			t_sum_length += txt[offset+i,2]
			a_sum_length += txt[offset+i,3]
		data_length.append(sum_length/30)
		b_data_length.append(b_sum_length/30)
		t_data_length.append(t_sum_length/30)
		a_data_length.append(a_sum_length/30)
	else:
		data_length.append(txt[offset+set_no-1,1])
		b_data_length.append(txtb[offset+set_no-1,1])
		t_data_length.append(txt[offset+set_no-1,2])
		a_data_length.append(txt[offset+set_no-1,3])

print(data_length)
data_length=[1.053,1.052,1.047,1.044,1.04] #rft
#data_length=[1.4,1.32,1.3,1.29,1.285] #com
#data_length=[1.97,1.971,1.975,1.9775,1.98] #amd
print(b_data_length)
b_data_length=[1.24,1.175,1.12,1.09,1.08] #rft
#b_data_length=[1.61,1.47,1.41,1.375,1.36] #com
#b_data_length=[1.47,1.505,1.48,1.497,1.53] #amd
print(t_data_length)
t_data_length=[1.053,1.052,1.047,1.044,1.04] #rft
#t_data_length=[1.92,2.04,2.025,2.015,2.01] #com
print(a_data_length)
plt.plot(ns,data_length,label=r'\textsc{Lpa}',marker='o',linewidth=2,color='red')
plt.plot(ns,b_data_length,label=r'\textsc{Batch}',marker='o',linestyle='-',color='darkblue')
if model != "pow":
	plt.plot(ns,t_data_length,label=r'\textsc{MinTime}',marker='o',linestyle='-',linewidth=1,color='darkgreen')
#plt.plot(ns,a_data_length,label=r'\textsc{MinArea}',marker='o',linestyle='-',color='gold')
plt.xlabel(r"$n$")
plt.ylabel(r"Normalized makespan")
plt.xlim([100,1000])
plt.xticks(ns)
#plt.ylim(ymin=1)
plt.ylim([1,3.5])
plt.title(r"roofline")
plt.legend()
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
plt.savefig("figs/njobsv2_"+model+"_"+str(set_no)+"_"+l+"_"+str(p)+".pdf")
print("figs/njobsv2_"+model+"_"+str(set_no)+"_"+l+"_"+str(p)+".pdf saved.")
plt.clf()

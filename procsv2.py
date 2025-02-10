import sys
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.size']=20
mpl.rcParams['legend.fontsize']=12
mpl.rcParams['legend.framealpha']=0.1
import matplotlib.pyplot as plt
import numpy as np

procs = [1000,2000,3000,4000,5000,6250,7500,8750,10000,11250,12500,13750,15000]
#procs = [1000,2000,3000,4000,5000,7500,12500,15000]
print(procs)
set_no = int(sys.argv[1])
model = sys.argv[2]
l=str(sys.argv[3])
big=int(sys.argv[4])


data_length = []
b_data_length = []
t_data_length = []
a_data_length = []
for p in procs:
	if model == "rft":
		txt = np.loadtxt("sample_rooftop/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_rooftop/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_rooftop/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_rooftop/area_sumv2_"+l+"_"+str(p))
	elif model == "amd":
		txt = np.loadtxt("sample_amdahl/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_amdahl/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_amdahl/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_amdahl/area_sumv2_"+l+"_"+str(p))
	elif model == "com":
		txt = np.loadtxt("sample_com/list_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_com/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_com/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_com/area_sum_"+l+"_"+str(p))
	elif model == "mix":
		txt = np.loadtxt("sample_mixst/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_mixst/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_mixst/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_mixst/area_sumv2_"+l+"_"+str(p))
	elif model == "mixlc":
		txt = np.loadtxt("sample_mixlc/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_mixlc/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_mixlc/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_mixlc/area_sumv2_"+l+"_"+str(p))
	elif model == "pow":
		txt = np.loadtxt("sample_pow/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_pow/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_pow/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_pow/area_sumv2_"+l+"_"+str(p))
	else:
		print("Modele inconnu")
		exit()
	if set_no == 0:
		sum_length = 0
		b_sum_length = 0
		a_sum_length = 0
		t_sum_length = 0
		if big:
			for i in range(30,60):
				sum_length += txt[i-30]
				b_sum_length += txtb[i-30]
				t_sum_length += txtT[i-30]
				a_sum_length += txtA[i-30]
				
		else:
			for i in range(30):
				sum_length += txt[i]
				b_sum_length += txtb[i]
				t_sum_length += txtT[i]
				a_sum_length += txtA[i]
		data_length.append(sum_length/30)
		b_data_length.append(b_sum_length/30)
		t_data_length.append(t_sum_length/30)
		a_data_length.append(a_sum_length/30)
	else:
		if big:
			data_length.append(txt[30+set_no-1,1])
			b_data_length.append(txtb[30+set_no-1,1])
			t_data_length.append(txtT[30+set_no-1,1])
			a_data_length.append(txtA[30+set_no-1,1])
		else:
			data_length.append(txt[set_no-1,1])
			b_data_length.append(txtb[set_no-1,1])
			t_data_length.append(txtT[set_no-1,1])
			a_data_length.append(txtA[set_no-1,1])

print(data_length)
data_length=[1.04,1.08,1.117,1.146,1.148,1.083,1.045,1.05,1.028,1.026,1.027,1.015,1.013] #rft
#data_length=[1.375,1.275,1.26,1.265,1.27,1.278,1.284,1.292,1.3,1.31,1.32,1.33,1.34] #com
#data_length=[1.7,2.05,2.05,2.01,1.99,1.98,1.975,1.972,1.97,1.97,1.97,1.97,1.97] #amd
print(b_data_length)
b_data_length=[1.12,1.141,1.144,1.166,1.168,1.124,1.118,1.132,1.141,1.167,1.195,1.215,1.24] #rft
#b_data_length=[1.51,1.39,1.28,1.3,1.32,1.34,1.36,1.38,1.4,1.42,1.44,1.46,1.48] #com
#b_data_length=[1.78,1.755,1.57,1.51,1.49,1.482,1.48,1.479,1.479,1.479,1.479,1.479,1.479] #amd
print(t_data_length)
t_data_length=[1.04,1.08,1.117,1.146,1.148,1.083,1.045,1.05,1.028,1.026,1.027,1.015,1.013] #rft
#t_data_length=[2.1,2.17,2.06,2.05,2.04,2.032,2.033,2.035,2.037,2.04,2.044,2.047,2.05] #com
print(a_data_length)
plt.plot(procs,data_length,label=r'\textsc{Lpa}',marker='o',color='red',linewidth=2)
plt.plot(procs,b_data_length,label=r'\textsc{Batch}',marker='o',linestyle='--',color='darkblue')
if model != "pow":
	plt.plot(procs,t_data_length,label=r'\textsc{MinTime}',marker='o',linestyle=':',color='darkgreen',linewidth=1)
#plt.plot(procs,a_data_length,label=r'\textsc{MinArea}',marker='o',linestyle='-.',color='gold')
plt.xlabel(r"$P$")
plt.ylabel(r"Normalized makespan")
plt.xlim([1000,15000])
plt.ylim([1,2.7])
#plt.yscale("log")
plt.title(r"roofline")
plt.legend()
plt.gcf().subplots_adjust(bottom=0.45)
plt.tight_layout()
if big:
	plt.savefig("figs/procsv2_500_"+model+"_"+str(set_no)+"_"+l+".pdf")
	print("figs/procsv2_500_"+model+"_"+str(set_no)+"_"+l+".pdf saved.")
else:
	plt.savefig("figs/procsv2_100_"+model+"_"+str(set_no)+"_"+l+".pdf")
	print("figs/procsv2_100_"+model+"_"+str(set_no)+"_"+l+".pdf saved.")
plt.clf()

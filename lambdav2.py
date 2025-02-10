import sys
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.size']=20
mpl.rcParams['legend.fontsize']=12
mpl.rcParams['legend.framealpha']=0.1
import matplotlib.pyplot as plt
import numpy as np

lambdas=["1e-8","2e-8","5e-8","1e-7","2e-7","5e-7","1e-6"]
lambdasInt=[1e-8,2e-8,5e-8,1e-7,2e-7,5e-7,1e-6]
print(lambdas)
set_no = int(sys.argv[1])
model = sys.argv[2]
print(model)
p=int(sys.argv[3])
big=int(sys.argv[4])


data_length = []
b_data_length = []
t_data_length = []
a_data_length = []
for l in lambdas:
	if model == "rft":
		txt = np.loadtxt("sample_rooftop/list_sumv2_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_rooftop/batch_sumv2_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_rooftop/time_sumv2_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_rooftop/area_sumv2_"+str(l)+"_"+str(p))
	elif model == "amd":
		txt = np.loadtxt("sample_amdahl/list_sum_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_amdahl/batch_list_sum_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_amdahl/time_sum_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_amdahl/area_sum_"+str(l)+"_"+str(p))
	elif model == "com":
		txt = np.loadtxt("sample_com/list_sum_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_com/batch_list_sum_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_com/time_sum_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_com/area_sum_"+str(l)+"_"+str(p))
	elif model == "mixlc":
		txt = np.loadtxt("sample_mixlc/list_sumv2_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_mixlc/batch_sumv2_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_mixlc/time_sumv2_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_mixlc/area_sumv2_"+str(l)+"_"+str(p))
	elif model == "mix":
		txt = np.loadtxt("sample_mixst/list_sumv2_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_mixst/batch_sumv2_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_mixst/time_sumv2_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_mixst/area_sumv2_"+str(l)+"_"+str(p))
	elif model == "pow":
		txt = np.loadtxt("sample_pow/list_sumv2_"+str(l)+"_"+str(p))
		txtb = np.loadtxt("sample_pow/batch_sumv2_"+str(l)+"_"+str(p))
		txtT = np.loadtxt("sample_pow/time_sumv2_"+str(l)+"_"+str(p))
		txtA = np.loadtxt("sample_pow/area_sumv2_"+str(l)+"_"+str(p))
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
				sum_length += txt[i,1]
				b_sum_length += txtb[i,1]
				t_sum_length += txtT[i,1]
				a_sum_length += txtA[i,1]
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
data_length=[1.047,1.0475,1.046,1.047,1.0465,1.0475,1.055] #roof
#data_length=[1.27,1.272,1.275,1.278,1.283,1.3,1.33] #com
#data_length=[1.97,1.97,1.97,1.972,1.97,1.97,1.97] #amd

print(b_data_length)
b_data_length=[1.094,1.1,1.098,1.115,1.145,1.205,1.26] #roof
#b_data_length=[1.2,1.22,1.252,1.38,1.48,1.68,1.87] #com
#b_data_length=[1.28,1.44,1.55,1.49,1.62,1.66,1.67] #amd
print(t_data_length)
t_data_length=[1.047,1.0475,1.046,1.047,1.0465,1.0475,1.055] #roof
#t_data_length=[2.03,2.028,2.031,2.033,2.03,2.031,1.99] #com
print(a_data_length)
plt.plot(lambdasInt,data_length,label=r'\textsc{Lpa}',marker='o',linewidth=2,color='red')
plt.plot(lambdasInt,b_data_length,label=r'\textsc{Batch}',marker='o',linestyle='-',color='darkblue')
if model != "pow":
	plt.plot(lambdasInt,t_data_length,label=r'\textsc{MinTime}',marker='o',linestyle='-',linewidth=1,color='darkgreen')
#plt.plot(lambdasInt,a_data_length,label=r'\textsc{MinArea}',marker='o',linestyle='-',color='gold')
plt.xlabel(r"$\lambda$")
plt.ylabel(r"Normalized makespan")
plt.xlim([1e-8,1e-6])
#plt.ylim(ymin=1)
plt.xscale("log")
plt.ylim([1.0,3.7])
#plt.yscale("log")
plt.title(r"roofline")
plt.legend()
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
if big:
	plt.savefig("figs/lambdasv2_500_"+model+"_"+str(set_no)+"_"+str(p)+".pdf")
	print("figs/lambdasv2_500_"+model+"_"+str(set_no)+"_"+str(p)+".pdf saved.")
else:
	plt.savefig("figs/lambdasv2_100_"+model+"_"+str(set_no)+"_"+str(p)+".pdf")
	print("figs/lambdasv2_100_"+model+"_"+str(set_no)+"_"+str(p)+".pdf saved.")
plt.clf()

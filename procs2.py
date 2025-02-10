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
		txt = np.loadtxt("sample_rooftop/list_sum_"+l+"_"+str(p))
		txtB = np.loadtxt("sample_rooftop/bound_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_rooftop/batch_list_sum_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_rooftop/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_rooftop/area_sum_"+l+"_"+str(p))
	elif model == "amd":
		txt = np.loadtxt("sample_amdahl/list_sum_"+l+"_"+str(p))
		txtB = np.loadtxt("sample_amdahl/bound_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_amdahl/batch_list_sum_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_amdahl/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_amdahl/area_sum_"+l+"_"+str(p))
	elif model == "com":
		txt = np.loadtxt("sample_com/list_sum_"+l+"_"+str(p))
		txtB = np.loadtxt("sample_com/bound_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_com/batch_list_sum_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_com/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_com/area_sum_"+l+"_"+str(p))
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
				sum_length += txt[i,1]/txtB[i,1]
				b_sum_length += txtb[i,1]/txtB[i,1]
				t_sum_length += txtT[i,1]/txtB[i,1]
				a_sum_length += txtA[i,1]/txtB[i,1]
				
		else:
			for i in range(30):
				sum_length += txt[i,1]/txtB[i,1]
				b_sum_length += txtb[i,1]/txtB[i,1]
				t_sum_length += txtT[i,1]/txtB[i,1]
				a_sum_length += txtA[i,1]/txtB[i,1]
		data_length.append(sum_length/30)
		b_data_length.append(b_sum_length/30)
		t_data_length.append(t_sum_length/30)
		a_data_length.append(a_sum_length/30)
	else:
		if big:
			data_length.append(txt[30+set_no-1,1]/txtB[30+set_no-1,1])
			b_data_length.append(txtb[30+set_no-1,1]/txtB[30+set_no-1,1])
			t_data_length.append(txtT[30+set_no-1,1]/txtB[30+set_no-1,1])
			a_data_length.append(txtA[30+set_no-1,1]/txtB[30+set_no-1,1])
		else:
			data_length.append(txt[set_no-1,1]/txtB[set_no-1,1])
			b_data_length.append(txtb[set_no-1,1]/txtB[set_no-1,1])
			t_data_length.append(txtT[set_no-1,1]/txtB[set_no-1,1])
			a_data_length.append(txtA[set_no-1,1]/txtB[set_no-1,1])

print(data_length)
print(b_data_length)
print(t_data_length)
print(a_data_length)
plt.plot(procs,data_length,label=r'\textsc{Lpa}',marker='o',color='red',linewidth=2)
plt.plot(procs,b_data_length,label=r'\textsc{Batch}',marker='o',linestyle='--',color='darkblue')
if model != "amd":
	plt.plot(procs,t_data_length,label=r'\textsc{MinTime}',marker='o',linestyle=':',color='darkgreen',linewidth=1)
#plt.plot(procs,a_data_length,label=r'\textsc{MinArea} / \textsc{Lpt}',marker='o',linestyle='-.',color='gold')
plt.xlabel(r"$P$")
plt.ylabel(r"Normalized makespan")
plt.xlim([1000,15000])
#plt.ylim([0.9,2.2])
plt.ylim(ymin=1)
#plt.title(r"Normalized makespan when $P$ varies for list-scheduling")
plt.legend()
plt.gcf().subplots_adjust(bottom=0.15)
plt.tight_layout()
if big:
	plt.savefig("figs/procs_500_"+model+"_"+str(set_no)+"_"+l+".pdf")
	print("figs/procs_500_"+model+"_"+str(set_no)+"_"+l+".pdf saved.")
else:
	plt.savefig("figs/procs_100_"+model+"_"+str(set_no)+"_"+l+".pdf")
	print("figs/procs_100_"+model+"_"+str(set_no)+"_"+l+".pdf saved.")
plt.clf()

import sys
import matplotlib as mpl
mpl.rcParams['text.usetex']=True
mpl.rcParams['font.size']=17
mpl.rcParams['legend.fontsize']=12
mpl.rcParams['legend.framealpha']=0.1
import matplotlib.pyplot as plt
import numpy as np

procs = [1000,2000,3000,4000,5000,7500,10000,12500,15000]
#procs = [1000,2000,3000,4000,5000,7500,12500,15000]
print(procs)
set_no = int(sys.argv[1])
model = sys.argv[2]
l="1e-7"
clean=int(sys.argv[3])

#P = 1000, 10000
#n = 100, 500
'''
x_length = [14*i for i in range(1)]
x_area = [14*i+1 for i in range(1)]
x_procs = [14*i+2 for i in range(1)]
b_x_length = [14*i+3 for i in range(1)]
b_x_area = [14*i+4 for i in range(1)]
b_x_procs = [14*i+5 for i in range(1)]
t_x_length = [14*i+6 for i in range(1)]
t_x_area = [14*i+7 for i in range(1)]
t_x_procs = [14*i+8 for i in range(1)]
a_x_length = [14*i+9 for i in range(1)]
a_x_procs = [14*i+10 for i in range(1)]
'''
x_length = [14*i for i in range(4)]
x_area = [14*i+1 for i in range(4)]
x_procs = [14*i+2 for i in range(4)]
b_x_length = [14*i+3 for i in range(4)]
b_x_area = [14*i+4 for i in range(4)]
b_x_procs = [14*i+5 for i in range(4)]
t_x_length = [14*i+6 for i in range(4)]
t_x_area = [14*i+7 for i in range(4)]
t_x_procs = [14*i+8 for i in range(4)]
a_x_length = [14*i+9 for i in range(4)]
a_x_procs = [14*i+10 for i in range(4)]

data_length = []
data_area = []
data_procs = []
t_data_length = []
t_data_area = []
t_data_procs = []
b_data_length = []
b_data_area = []
b_data_procs = []
a_data_length = []
a_data_procs = []
dataMin_length = []
dataMin_area = []
dataMin_procs = []
t_dataMin_length = []
t_dataMin_area = []
t_dataMin_procs = []
b_dataMin_length = []
b_dataMin_area = []
b_dataMin_procs = []
a_dataMin_length = []
a_dataMin_procs = []
dataMax_length = []
dataMax_area = []
dataMax_procs = []
t_dataMax_length = []
t_dataMax_area = []
t_dataMax_procs = []
b_dataMax_length = []
b_dataMax_area = []
b_dataMax_procs = []
a_dataMax_length = []
a_dataMax_procs = []
for p in [1000,10000]:
	if model == "rft":
		modelname="rooftop"
		txt = np.loadtxt("sample_rooftop/list_sumv2_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_rooftop/batch_sumv2_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_rooftop/time_sumv2_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_rooftop/area_sumv2_"+l+"_"+str(p))
		txtMax = np.loadtxt("sample_rooftop/list_max_"+l+"_"+str(p))
		txtbMax = np.loadtxt("sample_rooftop/batch_max_"+l+"_"+str(p))
		txtTMax = np.loadtxt("sample_rooftop/time_max_"+l+"_"+str(p))
		txtAMax = np.loadtxt("sample_rooftop/area_max_"+l+"_"+str(p))
		txtMin = np.loadtxt("sample_rooftop/list_min_"+l+"_"+str(p))
		txtbMin = np.loadtxt("sample_rooftop/batch_min_"+l+"_"+str(p))
		txtTMin = np.loadtxt("sample_rooftop/time_min_"+l+"_"+str(p))
		txtAMin = np.loadtxt("sample_rooftop/area_min_"+l+"_"+str(p))
	elif model == "amd":
		txt = np.loadtxt("sample_amdahl/list_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_amdahl/batch_list_sum_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_amdahl/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_amdahl/area_sum_"+l+"_"+str(p))
		txtMax = np.loadtxt("sample_amdahl/list_max_"+l+"_"+str(p))
		txtbMax = np.loadtxt("sample_amdahl/batch_max_"+l+"_"+str(p))
		txtTMax = np.loadtxt("sample_amdahl/time_max_"+l+"_"+str(p))
		txtAMax = np.loadtxt("sample_amdahl/area_max_"+l+"_"+str(p))
		txtMin = np.loadtxt("sample_amdahl/list_min_"+l+"_"+str(p))
		txtbMin = np.loadtxt("sample_amdahl/batch_min_"+l+"_"+str(p))
		txtTMin = np.loadtxt("sample_amdahl/time_min_"+l+"_"+str(p))
		txtAMin = np.loadtxt("sample_amdahl/area_min_"+l+"_"+str(p))
	elif model == "com":
		txt = np.loadtxt("sample_com/list_sum_"+l+"_"+str(p))
		txtb = np.loadtxt("sample_com/batch_list_sum_"+l+"_"+str(p))
		txtT = np.loadtxt("sample_com/time_sum_"+l+"_"+str(p))
		txtA = np.loadtxt("sample_com/area_sum_"+l+"_"+str(p))
		txtMax = np.loadtxt("sample_com/list_max_"+l+"_"+str(p))
		txtbMax = np.loadtxt("sample_com/batch_max_"+l+"_"+str(p))
		txtTMax = np.loadtxt("sample_com/time_max_"+l+"_"+str(p))
		txtAMax = np.loadtxt("sample_com/area_max_"+l+"_"+str(p))
		txtMin = np.loadtxt("sample_com/list_min_"+l+"_"+str(p))
		txtbMin = np.loadtxt("sample_com/batch_min_"+l+"_"+str(p))
		txtTMin = np.loadtxt("sample_com/time_min_"+l+"_"+str(p))
		txtAMin = np.loadtxt("sample_com/area_min_"+l+"_"+str(p))
	else:
		print("Modele inconnu")
		exit()
	if set_no == 0:
		
		sum_length = 0
		sum_area = 0
		sum_procs = 0
		b_sum_length = 0
		b_sum_area = 0
		b_sum_procs = 0
		a_sum_length = 0
		a_sum_procs = 0
		t_sum_length = 0
		t_sum_area = 0
		t_sum_procs = 0	
		max_length = 0
		max_area = 0
		max_procs = 0
		b_max_length = 0
		b_max_area = 0
		b_max_procs = 0
		a_max_length = 0
		a_max_procs = 0
		t_max_length = 0
		t_max_area = 0
		t_max_procs = 0	
		min_length = 1000
		min_area = 1000
		min_procs = 1000
		b_min_length = 1000
		b_min_area = 1000
		b_min_procs = 1000
		a_min_length = 1000
		a_min_procs = 1000
		t_min_length = 1000
		t_min_area = 1000
		t_min_procs = 1000	
		for i in range(30):
			sum_length += txt[i,1]
			sum_area += txt[i,2]
			sum_procs += txt[i,3]
			b_sum_length += txtb[i,0]
			b_sum_area += txtb[i,1]
			b_sum_procs += txtb[i,2]
			t_sum_length += txtT[i,0]
			t_sum_area += txtT[i,1]
			t_sum_procs += txtT[i,2]
			a_sum_length += txtA[i,0]
			a_sum_procs += txtA[i,1]
			max_length = max(txtMax[i,1],max_length)
			max_area = max(txtMax[i,2],max_area)
			max_procs = max(txtMax[i,3],max_procs)
			b_max_length = max(txtbMax[i,0],b_max_length)
			b_max_area = max(txtbMax[i,1],b_max_area)
			b_max_procs = max(txtbMax[i,2],b_max_procs)
			t_max_length = max(txtTMax[i,0],t_max_length)
			t_max_area = max(txtTMax[i,1],t_max_area)
			t_max_procs = max(txtTMax[i,2],t_max_procs)
			a_max_length = max(txtAMax[i,0],a_max_length)
			a_max_procs = max(txtAMax[i,1],a_max_procs)
			min_length = min(txtMin[i,1],min_length)
			min_area = min(txtMin[i,2],min_area)
			min_procs = min(txtMin[i,3],min_procs)
			b_min_length = min(txtbMin[i,0],b_min_length)
			b_min_area = min(txtbMin[i,1],b_min_area)
			b_min_procs = min(txtbMin[i,2],b_min_procs)
			t_min_length = min(txtTMin[i,0],t_min_length)
			t_min_area = min(txtTMin[i,1],t_min_area)
			t_min_procs = min(txtTMin[i,2],t_min_procs)
			a_min_length = min(txtAMin[i,0],a_min_length)
			a_min_procs = min(txtAMin[i,1],a_min_procs)
		data_length.append(sum_length/30)
		data_area.append(sum_area/30)
		data_procs.append(sum_procs/30)
		b_data_length.append(b_sum_length/30)
		b_data_area.append(b_sum_area/30)
		b_data_procs.append(b_sum_procs/30)
		t_data_length.append(t_sum_length/30)
		t_data_area.append(t_sum_area/30)
		t_data_procs.append(t_sum_procs/30)
		a_data_length.append(a_sum_length/30)
		a_data_procs.append(a_sum_procs/30)
		dataMax_length.append(max_length-sum_length/30)
		dataMax_area.append(max_area-sum_area/30)
		dataMax_procs.append(max_procs-sum_procs/30)
		b_dataMax_length.append(b_max_length-b_sum_length/30)
		b_dataMax_area.append(b_max_area-b_sum_area/30)
		b_dataMax_procs.append(b_max_procs-b_sum_procs/30)
		t_dataMax_length.append(t_max_length-t_sum_length/30)
		t_dataMax_area.append(t_max_area-t_sum_area/30)
		t_dataMax_procs.append(t_max_procs-t_sum_procs/30)
		a_dataMax_length.append(a_max_length-a_sum_length/30)
		a_dataMax_procs.append(a_max_procs-a_sum_procs/30)
		dataMin_length.append(0)
		dataMin_area.append(0)
		dataMin_procs.append(0)
		b_dataMin_length.append(0)
		b_dataMin_area.append(0)
		b_dataMin_procs.append(0)
		t_dataMin_length.append(0)
		t_dataMin_area.append(0)
		t_dataMin_procs.append(0)
		a_dataMin_length.append(0)
		a_dataMin_procs.append(0)
		
		sum_length = 0
		sum_area = 0
		sum_procs = 0
		b_sum_length = 0
		b_sum_area = 0
		b_sum_procs = 0
		a_sum_length = 0
		a_sum_procs = 0
		t_sum_length = 0
		t_sum_area = 0
		t_sum_procs = 0	
		max_length = 0
		max_area = 0
		max_procs = 0
		b_max_length = 0
		b_max_area = 0
		b_max_procs = 0
		a_max_length = 0
		a_max_procs = 0
		t_max_length = 0
		t_max_area = 0
		t_max_procs = 0	
		min_length = 1000
		min_area = 1000
		min_procs = 1000
		b_min_length = 1000
		b_min_area = 1000
		b_min_procs = 1000
		a_min_length = 1000
		a_min_procs = 1000
		t_min_length = 1000
		t_min_area = 1000
		t_min_procs = 1000	
		for i in range(30,60):
			sum_length += txt[i,1]
			sum_area += txt[i,2]
			sum_procs += txt[i,3]
			b_sum_length += txtb[i,0]
			b_sum_area += txtb[i,1]
			b_sum_procs += txtb[i,2]
			t_sum_length += txtT[i,0]
			t_sum_area += txtT[i,1]
			t_sum_procs += txtT[i,2]
			a_sum_length += txtA[i,0]
			a_sum_procs += txtA[i,1]
			max_length = max(txtMax[i,1],max_length)
			max_area = max(txtMax[i,2],max_area)
			max_procs = max(txtMax[i,3],max_procs)
			b_max_length = max(txtbMax[i,0],b_max_length)
			b_max_area = max(txtbMax[i,1],b_max_area)
			b_max_procs = max(txtbMax[i,2],b_max_procs)
			t_max_length = max(txtTMax[i,0],t_max_length)
			t_max_area = max(txtTMax[i,1],t_max_area)
			t_max_procs = max(txtTMax[i,2],t_max_procs)
			a_max_length = max(txtAMax[i,0],a_max_length)
			a_max_procs = max(txtAMax[i,1],a_max_procs)
			min_length = min(txtMin[i,1],min_length)
			min_area = min(txtMin[i,2],min_area)
			min_procs = min(txtMin[i,3],min_procs)
			b_min_length = min(txtbMin[i,0],b_min_length)
			b_min_area = min(txtbMin[i,1],b_min_area)
			b_min_procs = min(txtbMin[i,2],b_min_procs)
			t_min_length = min(txtTMin[i,0],t_min_length)
			t_min_area = min(txtTMin[i,1],t_min_area)
			t_min_procs = min(txtTMin[i,2],t_min_procs)
			a_min_length = min(txtAMin[i,0],a_min_length)
			a_min_procs = min(txtAMin[i,1],a_min_procs)
		data_length.append(sum_length/30)
		data_area.append(sum_area/30)
		data_procs.append(sum_procs/30)
		b_data_length.append(b_sum_length/30)
		b_data_area.append(b_sum_area/30)
		b_data_procs.append(b_sum_procs/30)
		t_data_length.append(t_sum_length/30)
		t_data_area.append(t_sum_area/30)
		t_data_procs.append(t_sum_procs/30)
		a_data_length.append(a_sum_length/30)
		a_data_procs.append(a_sum_procs/30)
		dataMax_length.append(max_length-sum_length/30)
		dataMax_area.append(max_area-sum_area/30)
		dataMax_procs.append(max_procs-sum_procs/30)
		b_dataMax_length.append(b_max_length-b_sum_length/30)
		b_dataMax_area.append(b_max_area-b_sum_area/30)
		b_dataMax_procs.append(b_max_procs-b_sum_procs/30)
		t_dataMax_length.append(t_max_length-t_sum_length/30)
		t_dataMax_area.append(t_max_area-t_sum_area/30)
		t_dataMax_procs.append(t_max_procs-t_sum_procs/30)
		a_dataMax_length.append(a_max_length-a_sum_length/30)
		a_dataMax_procs.append(a_max_procs-a_sum_procs/30)
		dataMin_length.append(0)
		dataMin_area.append(0)
		dataMin_procs.append(0)
		b_dataMin_length.append(0)
		b_dataMin_area.append(0)
		b_dataMin_procs.append(0)
		t_dataMin_length.append(0)
		t_dataMin_area.append(0)
		t_dataMin_procs.append(0)
		a_dataMin_length.append(0)
		a_dataMin_procs.append(0)
	else:
		data_length.append(txt[set_no-1,1])
		data_area.append(txt[set_no-1,2])
		data_procs.append(txt[set_no-1,3])
		data.append(txtb[set_no-1,1])
		data.append(txtb[set_no-1,2])
		data.append(txtb[set_no-1,3])
		data.append(txtT[set_no-1,1])
		data.append(txtT[set_no-1,2])
		data.append(txtT[set_no-1,3])
		data.append(txtA[set_no-1,1])
		data.append(txtA[set_no-1,2])
		data_length.append(txt[30+set_no-1,1])
		data_area.append(txt[30+set_no-1,2])
		data_procs.append(txt[30+set_no-1,3])
		data.append(txtb[30+set_no-1,1])
		data.append(txtb[30+set_no-1,2])
		data.append(txtb[30+set_no-1,3])
		data.append(txtT[30+set_no-1,1])
		data.append(txtT[30+set_no-1,2])
		data.append(txtT[30+set_no-1,3])
		data.append(txtA[30+set_no-1,1])
		data.append(txtA[30+set_no-1,2])
if model=="rft":# and not(clean):
	fig = plt.figure(figsize=(8.7,4.8))

plt.bar(x_length,data_length,color='red',yerr=[dataMin_length,dataMax_length],label=r'\textsc{Lpa} / \textsc{LPT}')
plt.bar(x_area,data_area,color='orangered',yerr=[dataMin_area,dataMax_area],label=r'\textsc{Lpa} / \textsc{LA}')
plt.bar(x_procs,data_procs,color='darkorange',yerr=[dataMin_procs,dataMax_procs],label=r'\textsc{Lpa} / \textsc{HPA}')
plt.bar(b_x_length,b_data_length,color='darkblue',yerr=[b_dataMin_length,b_dataMax_length],label=r'\textsc{Batch} / \textsc{LPT}')
plt.bar(b_x_area,b_data_area,color='blue',yerr=[b_dataMin_area,b_dataMax_area],label=r'\textsc{Batch} / \textsc{LA}')
plt.bar(b_x_procs,b_data_procs,color='steelblue',yerr=[b_dataMin_procs,b_dataMax_procs],label=r'\textsc{Batch} / \textsc{HPA}')
if not(clean) or (model != "amd" and model != "pow"):
	plt.bar(t_x_length,t_data_length,color='darkgreen',yerr=[t_dataMin_length,t_dataMax_length],label=r'\textsc{MinTime} / \textsc{LPT}')
	plt.bar(t_x_area,t_data_area,color='green',yerr=[t_dataMin_area,t_dataMax_area],label=r'\textsc{MinTime} / \textsc{LA}')
	plt.bar(t_x_procs,t_data_procs,color='limegreen',yerr=[t_dataMin_procs,t_dataMax_procs],label=r'\textsc{MinTime} / \textsc{HPA}')
if not(clean):
	plt.bar(a_x_length,a_data_length,color='gold',yerr=[a_dataMin_length,a_dataMax_length],label=r'\textsc{MinArea} / \textsc{LPT}')
	plt.bar(a_x_procs,a_data_procs,color='yellow',yerr=[a_dataMin_procs,a_dataMax_procs],label=r'\textsc{MinArea} / \textsc{HPA}')
	plt.yscale('log')
plt.ylabel(r"Normalized makespan")
plt.ylim(ymin=1)
plt.xticks([5,19,33,47],[r"$P=1000$""\n""$n=100$",r"$P=1000$""\n""$n=500$",r"$P=10000$""\n""$n=100$",r"$P=10000$""\n""$n=500$"])
#plt.xticks([5],[r"$P=7500$""\n""$n=500$"])
plt.tick_params(axis='x', length = 0)
#plt.ylim([0.9,2.2])
#plt.title(r"Normalized makespan when $P$ varies for list-scheduling")
if model == "rft":# and not(clean):
	lgd = plt.legend(loc="upper left",bbox_to_anchor=(-0.68,0.9))
	plt.gcf().subplots_adjust(bottom=0.15,left=0.37)
else:
	plt.gcf().subplots_adjust(bottom=0.15)
#plt.tight_layout()
if clean:
	
	plt.savefig("figs/barplotv2_"+model+"_clean.pdf")
	print("figs/barplotv2_"+model+"_clean.pdf")
	'''
	plt.savefig("figs/barplot_simple_"+model+"_clean.pdf")
	'''
else:
	
	plt.savefig("figs/barplotv2_"+model+".pdf")#,bbox_extra_artists=(lgd,))
	print("figs/barplotv2_"+model+".pdf")
	'''
	plt.savefig("figs/barplot_simple_"+model+".pdf")
	'''
plt.clf()

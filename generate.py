import sys
import random
import math

def generate_mold(filename,t):
	f = open(filename,'w')
	for i in range(1000):
		length = random.uniform(5000,4000000) #total work
		param = 0
		if t == "rand":
			choice = random.randint(1,3)
			if choice == 1:
				tt = "com"
			elif choice == 2:
				tt = "rft"
			elif choice == 3:
				tt = "amd"
			elif choice == 4:
				tt = "mix"
		else:
			tt = t
		if tt == "com":
			power = random.randint(0,3)
			coef = random.uniform(1,2)
			param = (2**power)*coef
		elif tt == "amd":
			#draw log
			power = random.randint(-7,-2)
			coef = random.uniform(0,10)
			param = coef*(10**power)
		elif tt == "rft":
			param = random.randint(100,4000)

		if tt=="mix":
			
			param1 = random.randint(100,4000)

			power = random.randint(-7,-2)
			coef = random.uniform(0,10)
			param4 = coef*(10**power)

			param2=length*param4
			length=length*(1-param4)

			power = random.randint(0,3)
			coef = random.uniform(1,2)
			param3 = (2**power)*coef
			
			name = "Task"+str(i)
			f.write(name+" mld "+str(length)+" "+tt+" "+str(param1)+" "+str(param2)+" "+str(param3)+" "+"\n")

		else:
			name = "Task"+str(i)
			f.write(name+" mld "+str(length)+" "+tt+" "+str(param)+"\n")
	f.close()

generate_mold(sys.argv[1],sys.argv[2])


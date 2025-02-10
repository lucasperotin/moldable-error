#include <iostream>
#include <iomanip>
#include <functional>
#include <vector>
#include <math.h>
#include <queue>
#include <time.h>
#include <sys/time.h>
#include <fstream>
#include <string>
#include <limits>

#include <omp.h>

#define DEBUG(X) std::cerr << X << "\n";

template<typename T> void print_queue(T& q)
{
	while (!q.empty())
	{
		std::string s = q.top().getName();
		std::cout << s << " ";
		q.pop();
	}
	std::cout << "\n";
}

struct taskOptions
{
	int rooftopLimit; //Limit for rooftop model
	double amdhalSeq; //Sequential fraction of the job for Amdhal model
	double comComm; //Communication overhead for communication model
        std::tuple<int, double, double> mix; 
};

class Task
{
	public:
		Task(std::string n, double l, int t, long p, std::function<double(double,long,taskOptions)> s) : name(n), length(l), speedup(s), type(t), np(p), failures(0) {}
		Task(std::string n, double l, int t, long p, std::function<double(double,long,taskOptions)> s, taskOptions to) : name(n), length(l), speedup(s), type(t), np(p), failures(0), options(to) {}
		void setName(std::string s) { name = s; }
		std::string getName() const { return name; }
		double time() const;
		double getPriority() { return priority; }
		void setPriority(double p) { priority = p; }
		int getType() const { return type; }
  double getOption() const {
  			if (getType()==1)
  				return options.rooftopLimit;
  			else if (getType()==2)
  				return options.amdhalSeq;
  			else if (getType()==3)
  				return options.comComm;
			//	else if (getType()==4)
			// return 0;
  		}
		long getProcs() const { return np; }
		void setProcs(long n) { if (type) np=n; }
		int getFailures() const { return failures; }	
		int getMaxFailures() const { return maxFailures; }
		void setMaxFailures(int ff) { maxFailures = ff; }
		void fail() { failures++; }//std::cerr << name << " failed\n"; }
		void minimizeTime(int P);
	private:
		std::string name;
		double length;
		int type;
		long np;
		std::function<double(double,long,taskOptions)> speedup;
		double priority;
		int failures;
		int maxFailures;
		taskOptions options;
};

double Task::time() const
{
	if (type) //if moldable
		return speedup(length,np,options);
	else //if rigid
		return length;
}

void Task::minimizeTime(int P)
{
	if (getType() == 1) //rft
		setProcs(std::min(P,options.rooftopLimit));
	if (getType() == 2) //amd
		setProcs(P);
	if (getType() == 3) //com
	{
		setProcs(1);
		double tmptime = time();
		setProcs(std::min(P,(int)floor(sqrt(tmptime/options.comComm))));
	}
	if (getType()==4); //mix
	{
	  setProcs(1);
	  double tmptime=time();
	  int minProc=1;
	  for (int i=2;i<P;i++)
	    {
	      setProcs(i);
	      if (tmptime>time())
		{
		  minProc=i;
		  tmptime=time();
		}
	    }
	  setProcs(minProc);
	}
}

struct Result
{
	double exec_time;
	std::vector<Task> tasks;
};

double BoundAchievable(std::vector<Task>& jobs,int P, double X, double epsilon)
{
    double Z;
    int maxFail;
    maxFail=0;
    for (int i=0;i<jobs.size();i++)
      {
	maxFail=std::max(maxFail,jobs[i].getFailures());
      }
    Z=epsilon*X/(maxFail+1);
    //std::cout << "\n Z " << Z << "\n";
    int X2;
    X2=floor(X/Z)+maxFail;
    double tempValeursmin;
    //std::cout << X2 << " X2 \n";
    double Valeurs0[X2+1];
    double Valeursjmult[X2+1];
    double Valeursjrest[X2+1];
    double Score[X2+1];
    double tempValeurs;
    int tempChoicel;
    double tempSum;
    int pp1,pp2;
    double temptime,temparea1,temparea2;
    std::vector<int> interesting_processors;

    for (int k=0;k<X2+1;k++)
      {
	Score[k]=0;
      }

    for (int i = 0; i < jobs.size(); i++)
    {
      jobs[i].setProcs(1);
      temptime=jobs[i].time();
      interesting_processors={1};
      
      for (int l=2;l<P+1;l++)
	{
	  jobs[i].setProcs(l);
	  if (jobs[i].time()<temptime)
	    {
	      temptime=jobs[i].time();
	      interesting_processors.push_back(l);
	    }
	}
	//DEBUG("liste done")
      int k=0;
      int l=interesting_processors.size()-1;
      pp1=interesting_processors[l];
      jobs[i].setProcs(pp1);
      temptime=jobs[i].time()/Z;
      temparea1=pp1*temptime/P;
      while (k<temptime && k<X2+1)
	{
	  Valeurs0[k]=(X/Z+maxFail)*2;
	  k++;
	}
      pp2=pp1;
      temparea2=temparea1;

	/*for (auto tmpp:interesting_processors)
		DEBUG(tmpp+" ")
	DEBUG("\n")*/

      for (int l=interesting_processors.size()-2;l>-1;l--)
	{
	  pp1=interesting_processors[l];
	  jobs[i].setProcs(pp1);
	  temptime=jobs[i].time()/Z;
	  temparea1=pp1*temptime/P;
	  while(k<temptime and k<X2+1)
	    {
	      Valeurs0[k]=temparea2;
	      k++;
	    }
	  pp2=pp1;
	  temparea2=temparea1;
	}
      while(k<X2+1)
	{
	  Valeurs0[k]=temparea2;
	  k++;
	}
	//DEBUG("1st part")
      for (int k=1;k<X2+1;k++)
	{
	  Valeursjmult[k]=Valeurs0[k];
	  Valeursjrest[k]=0;
	}

      int numberep=jobs[i].getFailures()+1;
      
      
      while(numberep>0)
        {
	  if (numberep%2==1)
	    {
	      for (int k=X2;k>0;k--)
		{
	            tempValeursmin=(X/Z+jobs[i].getFailures()+1)*2;
                    for (int l=1;l<k;l++)
                    {
                         tempValeurs=Valeursjrest[l]+Valeursjmult[k-l];
                         if (tempValeurs<tempValeursmin)
                         {
                             tempValeursmin=tempValeurs;
                         }
                    }
		    Valeursjrest[k]=tempValeursmin;
		}
	        numberep-=1;
	    }
	  
	  else
	  {
            for (int k=X2; k>0; k--)
            {
	        tempValeursmin=(X/Z+jobs[i].getFailures()+1)*2;
                tempChoicel=0;
                for (int l=1;l<k;l++)
                {

                    tempValeurs=Valeursjmult[l]+Valeursjmult[k-l];
                    if (tempValeurs<tempValeursmin)
                    {
                        tempValeursmin=tempValeurs;
                    }
                }
		Valeursjmult[k]=tempValeursmin;
            }
	    numberep/=2;
	  }
        }
	for (int k=1;k<X2+1;k++)
	  {
	    Score[k]+=Valeursjrest[k];
	  }
    }

    
    for (int k=1;k<X2+1;k++)
      {
	if (Score[k]<k && k*Z<X)
	  {
	    return (k*Z);
	  }
      }
   return -1;
}


double lowerBoundOnSet(std::vector<Task>& jobs,int P, double epsilon1,double epsilon2)
{
    double Xmax1=0;
    double Xmax2=0;
    double Xmax=0;
    double Xmin=0;
    double Xmed=0;
    float b;

    for (int i = 0; i < jobs.size(); i++)
	{
        jobs[i].setProcs(1);
        Xmax1=std::max(Xmax1,jobs[i].time()*(jobs[i].getFailures()+1));
        Xmax2=Xmax2+jobs[i].time()*(jobs[i].getFailures()+1);
	}
    Xmax2=Xmax2/P;
    Xmax=std::max(Xmax1,Xmax2);
	//DEBUG("BoundAchievable begin")
    b=BoundAchievable(jobs, P, Xmax, epsilon1);
	//DEBUG("BoundAchievable done")
    Xmax=b;

    while (Xmin==0 or Xmax/Xmin>1+epsilon2)
    {
      //std::cout << "\n Xmin : " << Xmin << "\n Xmax :" << Xmax << "\n " << Xmax/Xmin;
        Xmed=(Xmax+Xmin)/2;
        
        b=BoundAchievable(jobs, P, Xmed, epsilon1);

        if (b==-1)
        {
            Xmin=Xmed;
        }
        else
        {
            Xmax=b;
        }
    }
    Xmax=Xmax/(1+epsilon1)/(1+epsilon2);
    return Xmax;
}

double simpleBound(std::vector<Task> jobs, int p)
{
	double maxTime = 0;
	double totalArea = 0.0;
	for (int i = 0; i < jobs.size(); i++)
	{
		jobs[i].setProcs(1);
		totalArea += jobs[i].time()*(jobs[i].getMaxFailures()+1);
		jobs[i].minimizeTime(p);
		if (jobs[i].time()*(jobs[i].getMaxFailures()+1) > maxTime)
			maxTime = jobs[i].time()*(jobs[i].getMaxFailures()+1);
	}
	return std::max(maxTime, totalArea/(double)p);
}

void printResults(Result r)
{
	std::cout << "Execution time: " << r.exec_time << " seconds.\n";
	for (auto t:r.tasks)
	{
		std::cout << "Task " << t.getName() << " failed " << t.getFailures() << " times.\n";
	}
}

void printResultsSHORT(Result r,int p)
{
	double Xmax=0;
	int nb_fails=0;
	std::vector<Task> jobs;
	
	for(auto t:r.tasks)
	  {
	    jobs.push_back(t);
	    nb_fails+=t.getMaxFailures();
	  }
	Xmax = simpleBound(jobs,p);
	
	std::cout << r.exec_time << " " << Xmax << " " << nb_fails << "\n";
	if (Xmax>r.exec_time)
	  {
	      std::cerr << "The lower bound is higher... from a factor " << Xmax/r.exec_time << "\n";
	  }
}

class Event
{
	public:
		Event(double t, bool tp, Task rt) : time(t), related_task(rt), type(tp) {}
		bool type; //false = start, true = end
		double time;
		Task related_task;
};
auto cmpEvent = [](Event left, Event right) { return left.time > right.time; };

void printEvent(Event e)
{
	std::cerr << e.related_task.getName() << " " << e.time << " " << e.type << "\n";
}
void printEventQueue(std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> eventQueue)
{
	return;
	std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> tmpEventQueue(cmpEvent);
	while (!eventQueue.empty())
	{
		Event e = eventQueue.top();
		eventQueue.pop();
		printEvent(e);
		tmpEventQueue.push(e);
	}
	while (!tmpEventQueue.empty())
	{
		Event e = tmpEventQueue.top();
		tmpEventQueue.pop();
		eventQueue.push(e);
	}
	
}

class Simulator
{
	public:
		Simulator(int s, double l, double v) : seed(s), lambda(l), verif_time(v) {}
		void init() { srand(seed); }
		bool fail(Task t);
		double getLambda() const { return lambda; }
		double getVerif() const { return verif_time; }
		void setSeed(double s) { seed=s; }
		void setLambda(double l) { lambda = l; }
		void setVerif(double v) { verif_time = v; }
		std::function<double(Task)> priority_fun;
		bool batch() { return batchVar>0; }
		bool naive() { return naiveVar>0; }
		bool shelves() { return shelfVar>0; }
		void setBatch(int b) { batchVar=b; }
		void setNaive(int b) { naiveVar=b; }
		void setShelves(int b) { shelfVar=b; }
	private:
		double seed;
		double lambda;
		double verif_time;
		int batchVar, naiveVar, shelfVar;
};

bool Simulator::fail(Task t)
{
	/*t.setProcs(1);
	double randVal = (double)rand()/(double)RAND_MAX;
	double proba = 1 - exp(-lambda * t.time());
	return randVal <= proba;*/
	return t.getFailures()<t.getMaxFailures();
}

bool checkAvailProc(double startTime, Task t, int startProcs, std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> eQ, Simulator* s)
{
	//std::cerr << "CHECK AVAIL PROC(" << startTime << "," << t.getName() << "," << startProcs << ")\n";

//	auto cmpEvent = [](Event left, Event right) { return left.time > right.time; };
//	std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> tmpEventQueue(cmpEvent);
	if (startProcs < t.getProcs())
		return false;
	if (eQ.empty())
		return true;
	Event e = eQ.top();
//	tmpEventQueue.push(e);
	eQ.pop();
	int currentProcs = startProcs;
	while (e.time < startTime + t.time() + s->getVerif())
	{
		if (e.type)
			currentProcs += e.related_task.getProcs();
		else
			currentProcs -= e.related_task.getProcs();
		if (currentProcs < t.getProcs())
		{
/*				while (!tmpEventQueue.empty())
				{
					eventQueue.push(tmpEventQueue.top());
					tmpEventQueue.pop();
				}*/
				return false;
		}
		if (!eQ.empty())
		{
			e = eQ.top();
			eQ.pop();
		} else
			break;
	}/*
	while (!tmpEventQueue.empty())
	{
		eventQueue.push(tmpEventQueue.top());
		tmpEventQueue.pop();
	}*/

	return true;
}

//batch = true => failed tasks are rescheduled after all tasks are done once
//naive = true => always schedule by decreasing priority (no fill)
//shelves = true => schedule tasks at once and then wait for them to be finished
Result simulate(Simulator& s, std::vector<Task> jobs, long nbProcs, std::string mParam)
{
	double cTime = 0;
	s.init();

	Result r;

	//Init the priority queue of jobs
	for (unsigned i=0; i<jobs.size(); i++)
		jobs[i].setPriority(s.priority_fun(jobs[i]));

	auto cmpTask = [](Task left, Task right) { return left.getPriority() < right.getPriority(); };
	std::priority_queue<Task,std::vector<Task>,decltype(cmpTask)> jobQueue(cmpTask);
	std::priority_queue<Task,std::vector<Task>,decltype(cmpTask)> tmpJobQueue(cmpTask);
	std::priority_queue<Task,std::vector<Task>,decltype(cmpTask)> reservedJobQueue(cmpTask); //to keep the list of reservations made (in order to delete them in case a task fails)
	std::priority_queue<Task,std::vector<Task>,decltype(cmpTask)> failedJobs(cmpTask); //for next batch

	for (Task j:jobs)
		jobQueue.push(j);

	long pAvail = nbProcs;

	//Create a priority queue for the events (ends of tasks)
//	auto cmpEvent = [](Event left, Event right) { return left.time > right.time; };
	std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> eventQueue(cmpEvent);
	std::priority_queue<Event,std::vector<Event>,decltype(cmpEvent)> tmpEventQueue(cmpEvent);

	int m;

	if (mParam == "zero")
		m = 0;
	else if (mParam == "one")
		m = 1;
	else if (mParam == "all")
		m = jobQueue.size();

	begin:
	//Schedule the first tasks

	int cpt = 0;
	while (!jobQueue.empty() && cpt < m)
	{
		Task first = jobQueue.top();
		jobQueue.pop();
		//find next event where to put the task
		
		//std::cerr << "Reserving for " << first.getName(); 

		//if (pAvail >= first.getProcs()) //OK to schedule now as others
		if (checkAvailProc(cTime,first,pAvail,eventQueue,&s)) //check until end of task because of the reservations
		{
			//std::cerr << " at cTime " << cTime << "\n";
			pAvail -= first.getProcs();
			Event e(cTime+first.time()+s.getVerif(),true,first);
			eventQueue.push(e);
			cpt++;
		} else { //or in the future
			int pAvailFuture = pAvail;
			while (!eventQueue.empty())
			{
				Event next = eventQueue.top();
				eventQueue.pop();
				tmpEventQueue.push(next);
				if (next.type)
					pAvailFuture += next.related_task.getProcs();
				else
					pAvailFuture -= next.related_task.getProcs();
				//if (pAvailFuture >= first.getProcs())
				if (checkAvailProc(next.time,first,pAvailFuture,eventQueue,&s))
				{
					//std::cerr << " at time " << next.time << "\n";
					Event e(next.time,false,first);
					eventQueue.push(e);
					Event f(next.time+first.time()+s.getVerif(),true,first);
					eventQueue.push(f);
					reservedJobQueue.push(first);
					cpt++; //number of reservations
					break;
				}
			}
		
			while (!tmpEventQueue.empty())
			{
				eventQueue.push(tmpEventQueue.top());
				tmpEventQueue.pop();
			}
		}
	}

	while (!jobQueue.empty())
	{
		Task first = jobQueue.top();
		jobQueue.pop();
		//if (first.getProcs() <= pAvail) //If enough procs available
		if (checkAvailProc(cTime,first,pAvail,eventQueue,&s))
		{
			//std::cerr << "Task " << first.getName() << " scheduled.\n";
			pAvail -= first.getProcs();
			Event e(cTime+first.time()+s.getVerif(),true,first);
			eventQueue.push(e);
		} else { //otherwise needs to be rescheduled later
			//std::cout << "Task " << first.getName() << " not scheduled.\n";
			tmpJobQueue.push(first);
			if (s.naive())
				break; //When naive is true, we stop as soon as we cannot add the next task with highest priority (like Turek paper)
		}
	}
	//Put the non-scheduled jobs back in the queue
	while (!tmpJobQueue.empty())
	{
		jobQueue.push(tmpJobQueue.top());
		tmpJobQueue.pop();
	}

	//MAIN LOOP
	while (!eventQueue.empty())
	{
		Event e = eventQueue.top();
		eventQueue.pop();
		cTime = e.time;
	
		//std::cout << "At time " << cTime << ":\n";
	
		if (e.type) { //It is an ending event so we need to schedule new tasks


		//std::cout << "Task " << e.related_task.getName() << " (execution #" << e.related_task.getFailures()+1 << ") finished at time " << e.time << "\n";
		if (s.fail(e.related_task)) //If it failed
		{
			//std::cerr << "Task " << e.related_task.getName() << " failed.\nBefore cancellation:\n";
			printEventQueue(eventQueue);
			e.related_task.fail();
			e.related_task.setPriority(s.priority_fun(e.related_task)); //USELESS SO FAR BUT WE MAY WANT TO CHANGE THE PRIORITY OF A FAILED TASK
			while (!reservedJobQueue.empty()) //CANCEL RESERVATIONS
			{
				Task job = reservedJobQueue.top();
				reservedJobQueue.pop();
				int nEvents = 0;
				while (!eventQueue.empty() && nEvents<2)
				{
					Event ee = eventQueue.top();
					eventQueue.pop();
					tmpEventQueue.push(ee);
					/*if (e.related_task.getName() != job.getName()) //discard the event if same name
					else
					{*/
					if (nEvents == 0 && e.type) //the task has already started as we didn't find the begin event
						break;
					nEvents++; //max 2 events for the same task
					if (nEvents == 2) //OK WE FOUND THE RESERVATION (OTHERWISE THE JOB IS RUNNING OR WAS ALREADY DONE) so we can put it back to being scheduled
						jobQueue.push(job); //back to being scheduled
					//}
				}
				while (!tmpEventQueue.empty()) //we put back the events saved
				{
					eventQueue.push(tmpEventQueue.top());
					tmpEventQueue.pop();
				}
			}
			//std::cerr << "After cancellation:\n";
			printEventQueue(eventQueue);
			if (s.batch())
				failedJobs.push(e.related_task);
			else
				jobQueue.push(e.related_task);
		} else {
			//std::cerr << e.related_task.getName() << " finished at time " << cTime << ". Remaining: " << jobQueue.size() << " tasks.\n";
			printEventQueue(eventQueue);
			r.tasks.push_back(e.related_task);
		}
		pAvail += e.related_task.getProcs();

		if (!s.shelves() || eventQueue.empty()) //RESCHEDULE ONLY IF NOT SHELVES OR SHELF FINISHED
		{
			//Schedule the next tasks
			if (mParam == "zero")
				m = 0;
			else if (mParam == "one")
				m = 1;
			else if (mParam == "all")
				m = jobQueue.size();

			int cpt = 0;
			while (!jobQueue.empty() && cpt < m)
			{
				Task first = jobQueue.top();
				jobQueue.pop();
				//std::cerr << "---Reserving for " << first.getName();
				//find next event where to put the task
				

				//if (pAvail >= first.getProcs()) //OK to schedule now as others
				if (checkAvailProc(cTime,first,pAvail,eventQueue,&s)) //check until end of task because of the reservations
				{
					//std::cerr << " at cTime " << cTime << "\n";
					pAvail -= first.getProcs();
					Event e(cTime+first.time()+s.getVerif(),true,first);
					eventQueue.push(e);
					cpt++;
				} else { //or in the future
					int pAvailFuture = pAvail;
					while (!eventQueue.empty())
					{
						Event next = eventQueue.top();
						eventQueue.pop();
						tmpEventQueue.push(next);
						if (next.type)
							pAvailFuture += next.related_task.getProcs();
						else
							pAvailFuture -= next.related_task.getProcs();
						//if (pAvailFuture >= first.getProcs())
						if (checkAvailProc(next.time,first,pAvailFuture,eventQueue,&s))
						{
							//std::cerr << " at time " << next.time << "\n";
							Event e(next.time,false,first);
							eventQueue.push(e);
							Event f(next.time+first.time()+s.getVerif(),true,first);
							eventQueue.push(f);
							reservedJobQueue.push(first);
							cpt++;
							break;
						}
					}
				
					while (!tmpEventQueue.empty())
					{
						eventQueue.push(tmpEventQueue.top());
						tmpEventQueue.pop();
					}
				}
			}

			//Other jobs to schedule (backfilling)
			while (!jobQueue.empty())
			{
				Task first = jobQueue.top();
				jobQueue.pop();
				//if (first.getProcs() <= pAvail) //If enough procs available
				if (checkAvailProc(cTime,first,pAvail,eventQueue,&s)) //check with the existing reservations
				{
					//std::cerr << "---Task " << first.getName() << " scheduled.\n";
					pAvail -= first.getProcs();
					Event e(cTime+first.time()+s.getVerif(),true,first);
					eventQueue.push(e);
				} else { //otherwise needs to be rescheduled later
					tmpJobQueue.push(first);
					if (s.naive())
						break;
				}
			}
			//Put the non-scheduled jobs back in the queue
			while (!tmpJobQueue.empty())
			{
				jobQueue.push(tmpJobQueue.top());
				tmpJobQueue.pop();
			}
		}

		} else { //it is a begin event so we account for the processors used
			pAvail -= e.related_task.getProcs();
		}
		
	}
	if (s.batch() && !failedJobs.empty()) //Batch mode and some jobs need to be re-executed
	{
		while (!failedJobs.empty())
		{
			jobQueue.push(failedJobs.top());
			failedJobs.pop();
		}
		goto begin; //start again with list of failed jobs
	}

	//std::cout << "Execution time: " << cTime << "\n";

	r.exec_time = cTime;

	return r;
}

double pTaskLength(Task t)
{
	return t.time();
}

double pTaskArea(Task t)
{
	return t.time()*t.getProcs();
}

double pTaskProcs(Task t)
{
	return t.getProcs();
}

double pTaskLengthS(Task t)
{
	return -t.time();
}

double pTaskAreaS(Task t)
{
	return -(t.time()*t.getProcs());
}

double pTaskProcsS(Task t)
{
	return -t.getProcs();
}

double pTaskRandom(Task t)
{
	return rand()/(double)RAND_MAX;
}

double linearSpeedup(double seq_length, long nb_proc, taskOptions to)
{
	return seq_length/(double)nb_proc;
}

double rooftopSpeedup(double seq_length, long nb_proc,taskOptions to)
{
	return std::max(seq_length/(double)nb_proc,seq_length/(double)to.rooftopLimit);
}

double amdahlSpeedup(double seq_length, long nb_proc, taskOptions to)
{
	return seq_length*(to.amdhalSeq+(1.0-to.amdhalSeq)/(double)nb_proc);
}

double comSpeedup(double seq_length, long nb_proc, taskOptions to)
{
	return seq_length/(double)nb_proc + (nb_proc-1)*to.comComm;
}

double mixSpeedup(double seq_length, long nb_proc, taskOptions to)
{
  int rooflimit=std::get<0>(to.mix);
  double amdval=std::get<1>(to.mix);
  double comval=std::get<2>(to.mix);
  return std::max(seq_length/(double)nb_proc,seq_length/(double)rooflimit)+amdval+(nb_proc-1)*comval;
}

int readInput(std::string filename, std::vector<Task> &v)
{
	//int batch,naive,shelf;
	int shelf;
	//double lambda,verif,length;
	double length;
	std::string name;
	std::string type;
	std::string speedup;
	int proc;

	double lower_bound = 0; //To evaluate the "optimal" makespan

	std::ifstream input(filename,std::ios::in);	
	//input >> lambda >> *p >> verif >> priority >> batch >> naive >> shelf;
	/*s->setLambda(lambda);
	s->setVerif(verif);
	s->setBatch(batch);
	s->setNaive(naive);
	s->setShelves(shelf);*/
	/*
	if (priority == "length")
		s->priority_fun = pTaskLength;
	else if (priority == "procs")
		s->priority_fun = pTaskProcs;
	else if (priority == "area")
		s->priority_fun = pTaskArea;
	else if (priority == "rlength")
		s->priority_fun = pTaskLengthS;
	else if (priority == "rprocs")
		s->priority_fun = pTaskProcsS;
	else if (priority == "rarea")
		s->priority_fun = pTaskAreaS;
	else if (priority == "rand")
		s->priority_fun = pTaskRandom;
	else {
		std::cerr << "Unrecognized priority function: " << priority << ".\n";
		return -1;
	}*/

	while (input >> name >> type >> length)
	{
		int intTmp;
		double doubleTmp;
		if (type != "rgd" && type != "mld") {
			std::cerr << "Undefined type of task.\n";
			return -1;
		}
		if (type == "rgd")
		{
			input >> proc;
			v.push_back(Task(name,length,0,proc,linearSpeedup));
			lower_bound += length*proc;
		}
		else
		{
			input >> speedup;
			if (speedup == "lin")
				v.push_back(Task(name,length,4,0,linearSpeedup));
			else if (speedup == "rft")
			{
				taskOptions to;
				input >> to.rooftopLimit;
				v.push_back(Task(name,length,1,0,rooftopSpeedup,to));
			}
			else if (speedup == "amd")
			{
				taskOptions to;
				input >> to.amdhalSeq;
				v.push_back(Task(name,length,2,0,amdahlSpeedup,to));
			}
			else if (speedup == "com")
			{
				taskOptions to;
				input >> to.comComm;
				v.push_back(Task(name,length,3,0,comSpeedup,to));
			}
			else if (speedup=="mix")
			  {
			    taskOptions to;
			    input >> std::get<0>(to.mix);
			    input >> std::get<1>(to.mix);
			    input >> std::get<2>(to.mix);
			    
			    v.push_back(Task(name,length,4,0,mixSpeedup,to));
			  }
			else
			{
				std::cerr << "Unrecognized speedup function.\n";
				return -1;
			}
		}
	}

	//std::cout << "A lower bound on the optimal makespan: " << lower_bound/(*p) << ".\n";

	input.close();

	return 0;
}

double computeR(double alpha, double beta, int P)
{
	double r = 0;
	if (alpha >= beta)
		r = 2*alpha;
	else
		r = (double)P/(double)(P-1)*alpha + (double)(P-2)/(double)(P-1)*beta;
	return r;
}

void assignProcs(std::vector<Task>& jobs,int P)
{
	for (int i = 0; i < jobs.size(); i++)
	{
		//First compute a_min and t_min
		jobs[i].setProcs(1);
		double tmin = jobs[i].time();
		double amin = tmin;
		double ttmp;
		for (int p = 2; p <= P; p++)
		{
			jobs[i].setProcs(p);
			ttmp = jobs[i].time();
			if (ttmp < tmin)
				tmin = ttmp;
			if (p*ttmp < amin)
				amin = p*ttmp;
		}
		//Then compute r_min and p_j
		jobs[i].setProcs(1);
		double alpha = jobs[i].time()/amin;
		double beta = jobs[i].time()/tmin;
		double rmin = computeR(alpha,beta,P);
		double rtmp, pj = 1;
		for (int p = 2; p <= P; p++)
		{
			jobs[i].setProcs(p);
			beta = jobs[i].time()/tmin;
			alpha = beta*p*tmin/amin;
			rtmp = computeR(alpha, beta,P);
			if (rtmp < rmin)
			{
				pj = p;
				rmin = rtmp;
			}
		}
		//assign p_j procs to job j
		jobs[i].setProcs(pj);
	}
}

void assignProcsTime(std::vector<Task>& jobs, int P)
{
	for (int i = 0; i < jobs.size(); i++)
	{
		jobs[i].minimizeTime(P);
	}
}

void assignProcsArea(std::vector<Task>& jobs, int P)
{
	for (int i = 0; i < jobs.size(); i++)
	{
		jobs[i].setProcs(1);
	}
}

void printAssign(std::vector<Task>& jobs)
{
	for (unsigned i=0; i<jobs.size(); i++)
		std::cout << jobs[i].getProcs() << " ";
	std::cout << "\n";
}

int main(int argc, char** argv)
{
	/*Task job1("test",1000,false,100,linearSpeedup);
	Task job2("bigtask",100000,false,1389,linearSpeedup);
	Task job3("shorttask",100,false,382,linearSpeedup);
	Task job4("rand",58443,false,1000,linearSpeedup);
	Task job5("mamamia",29847,false,2000,linearSpeedup);*/
	struct timeval time;
	gettimeofday(&time, NULL);
	double val = (time.tv_sec * 1000) + (time.tv_usec / 1000);
	Simulator s(val, 0,0);
	std::vector<Task> jobs;

	//./moldable jobs failures nb_iter num_procs priority_fun shelves?
	if (argc < 9)
	{
		std::cout << "Missing arguments: " << 9-argc << "\n";
		std::cout << argv[0] << " job_list failure_list nb_iter nb_procs priority batch naive shelves.\n";
		std::cout << "Priority is among: length, rlength, procs, rprocs, area, rarea, random.\n";
		exit(2);
	}

	//int ok = readInput(argv[1],&s,jobs,&num_procs);
	int ok = readInput(argv[1],jobs);
	std::ifstream failures(argv[2],std::ios::in);
	/*double lambda = atof(argv[2]);
	s.setLambda(lambda);*/
	int nb_iter = stoi(std::string(argv[3]));
	int num_procs = stoi(std::string(argv[4]));
	if (num_procs < 2)
	{
		std::cout << "You need at least 2 processors.\n";
		exit(2);
	}
	std::string priority = std::string(argv[5]);
	if (priority == "length")
		s.priority_fun = pTaskLength;
	else if (priority == "procs")
		s.priority_fun = pTaskProcs;
	else if (priority == "area")
		s.priority_fun = pTaskArea;
	else if (priority == "rlength")
		s.priority_fun = pTaskLengthS;
	else if (priority == "rprocs")
		s.priority_fun = pTaskProcsS;
	else if (priority == "rarea")
		s.priority_fun = pTaskAreaS;
	else if (priority == "rand")
		s.priority_fun = pTaskRandom;
	else {
		std::cerr << "Unrecognized priority function: " << priority << ".\n";
		return -1;
	}
	/*std::string assign_fun = std::string(argv[6]);
	if (assign_fun != "FF" && assign_fun != "EJ" && assign_fun != "FW") {
		std::cerr << "Unrecognized assignment cost function: " << assign_fun << ".\n";
		return -1;
	}*/
	s.setBatch(stoi(std::string(argv[6])));
	s.setNaive(stoi(std::string(argv[7])));
	s.setShelves(stoi(std::string(argv[8])));

	if (ok==0)
	{
		double avg_time = 0;
		for (int i=0; i<nb_iter; i++)
		{
			std::cerr << i << "\n";
			s.setSeed(val+i);
			int nb_fail_tmp = 0;
			for (unsigned j=0; j<jobs.size(); j++)
			{
				failures >> nb_fail_tmp;
				jobs[j].setMaxFailures(nb_fail_tmp);
			}


			assignProcs(jobs,num_procs);
			//DEBUG("assignment done");
			//printAssign(jobs);
			Result r=simulate(s,jobs,num_procs,"zero");
			
			//DEBUG("finished")
			{
				avg_time += r.exec_time;
			}
			printResultsSHORT(r,num_procs);
		}
		avg_time /= nb_iter;
		//std::cout << avg_time << "\n";
	}

	failures.close();

	return 0;
}

/*OLD FUNCTION

void assignProcs(std::vector<Task>& jobs,int P, double lambda, std::string cost_fun)
{
	for (auto& t:jobs)
		t.setProcs(1);
	bool better = true;
	double cost = std::numeric_limits<double>::max();
	std::vector<int> failure_free(jobs.size(),0);
	while (better)
	{
		unsigned int index = jobs.size();
		double min_cost = cost;	
		double cur_cost = 0;
		for (unsigned int i=0; i<jobs.size(); i++)
		{
			if (jobs[i].getProcs() < P)
			{
				jobs[i].setProcs(jobs[i].getProcs()+1);
				if (cost_fun == "FF") //Failure-Free Lower Bound
					cur_cost = FFLB(jobs,P,failure_free);
				else if (cost_fun == "FW") //Failure-Weighted Lower Bound
					cur_cost = FWLB(jobs,P,1,lambda);
				else if (cost_fun == "EJ") //Estimated Jobs Lower Bound
					cur_cost = EJLB(jobs,P,lambda);
				jobs[i].setProcs(jobs[i].getProcs()-1);
				//std::cerr << cur_cost << "\n";
				if (cur_cost <= min_cost)
				{
					min_cost = cur_cost;
					index = i;
				}
			}
		}
		better = (index != jobs.size());
		if (better)
			jobs[index].setProcs(jobs[index].getProcs()+1);
		cost = min_cost;
	}
	
}*/
/*
double FFLB(const std::vector<Task>& jobs,int P, std::vector<int>& failures) //failures are here to be more general as the function is used in FWLB with different failure scenarios
{
	if (failures.size() != jobs.size())
		std::cerr << "Jobs and Failures have different sizes: " << jobs.size() << " " << failures.size() << ".\n";
	double area = 0;
	double lb = 0;
	for (unsigned i=0; i<jobs.size(); i++)
	{
		area += jobs[i].getProcs()*(failures[i]+1)*jobs[i].time();
		if ((failures[i]+1)*jobs[i].time() > lb)
			lb = (failures[i]+1)*(jobs[i].time());
	}
	if (area/(double)P > lb) lb = area/(double)P;
	return lb;
}

double proba(const std::vector<Task>& jobs, std::vector<int>& failures, double lambda)
{
	if (failures.size() != jobs.size())
		std::cerr << "Jobs and Failures have different sizes.\n";
	double p = 1;
	for (int i=0; i<jobs.size(); i++)
	{
		double pi = 1;
		for (int j=0; j<failures[i]; j++)
			pi *= (1-exp(-lambda*jobs[i].getProcs()*jobs[i].time()));
		pi *= exp(-lambda*jobs[i].getProcs()*jobs[i].time());
		p *= pi;
	}
	return p;
}

int nextBitArray(std::vector<int>& t, int k) //no more than k failures in t
{
        //we always increment last bit but we need to know when it should be "reset"
        int nb_assigned = 0;
	int return_value = 1; //to say if we came back to 00...00
        for (unsigned i=0; i<t.size(); i++)
        {
                nb_assigned += t[i];
        }
        if (nb_assigned > k)
        {
                std::cerr << "Non standard bit array for value " << k << ".\n";
        }

        for (int i=(int)t.size()-1; i>=0; i--)
        {
                if (nb_assigned < k)
                {
                        t[i]++;
                        break;
                } else {
                        nb_assigned -= t[i];
                        t[i] = 0;
			if (i==0) //it means t[0] = nb_assigned so we reached the last step
				return_value = 0;
                }
        }

	return return_value;
}

double FWLB(const std::vector<Task>& jobs, int P, int k, double lambda) //k is the max number of faults
{
	std::vector<int> failures(jobs.size(),0);
	double total_proba = 0;
	double total_val = 0;
	do
	{
		double L = FFLB(jobs,P,failures);
		double Q = proba(jobs,failures,lambda);
		total_proba += Q;
		total_val += L*Q;
	} while (nextBitArray(failures,k));
	return total_val / total_proba;
}

double EJLB(const std::vector<Task>& jobs,int P,double lambda)
{
	double area = 0;
	double lb = 0;
	for (unsigned i=0; i<jobs.size(); i++)
	{
		double time = jobs[i].time()*exp(lambda*jobs[i].time()*jobs[i].getProcs());
		area += jobs[i].getProcs()*time;
		if (time > lb)
			lb = time;
	}
	if (area/(double)P > lb) lb = area/(double)P;
	return lb;
}*/

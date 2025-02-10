#include <iostream>
#include <iomanip>
#include <fstream>
#include <vector>
#include <string>
#include <math.h>

int main(int argc, char** argv)
{
	int nb_iter(stoi(std::string(argv[3])));
	std::string filename(argv[1]);
	std::string filename2(argv[2]);
	std::string tmp;
	std::ifstream input(filename, std::ios::in);
	std::ofstream output(filename2,std::ios::out);
	double lambda,length; int procs;
	int nb_tasks=0;
	double avg_area=0;
	int total_failures=0;
	for (int i=0; i<nb_iter; i++)
	{
		//input >> lambda;
		lambda = atof(argv[4]);
		//getline(input,tmp);
		input.ignore();
		while(input >> tmp >> tmp)
		{
			if (tmp == "mld")
			{
				input >> length >> tmp;
				if (tmp=="mix")
				  {
				    input>>tmp >> tmp>> tmp;
				  }
				else
				  {
				    input>>tmp;
				  }
				if (i==0) {
				avg_area += length;
				nb_tasks++; }
				bool failed = true;
				int nb_failures = 0;
				while (failed)
				{
					double r = (double)rand()/(double)RAND_MAX;
					double proba = 1-exp(lambda*length*(-1.0));
					failed = r < proba;
					if (failed)
						nb_failures++;	
				}
				output << nb_failures << " ";
				total_failures += nb_failures;
			} else {
				std::cerr << "Not a moldable task.\n";
				input.close();
				output.close();
				return -1;
			}
		}
		output << "\n";
		input.clear();
		input.seekg(0);
	}	
	std::cout << "Average number of failures : "<< (double)total_failures/(double)nb_iter << "\n";
	std::cout << "Average length = " << avg_area/nb_tasks << ".\n";
	input.close();
	output.close();
}

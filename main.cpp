//
//  main.cpp
//  transitionMatrix
//
//  Created by Will Pearse on 31/10/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//

#include <iostream>
#include <stdlib.h>
#include "communitySimulator.h"
#include "communityData.h"

using namespace std;
using namespace boost::numeric;

int main (int argc, const char * argv[])
{
    if(argc == 2)
    {
        //Only works for the first community...
        cout << "Loading dataset: " << argv[1] << endl;
        DataSet data(argv[1]);
        data.communities[0].optimise();
        cout << "First community -> second community event matrix:" << endl;
        data.communities[0].print_event_matrix(0,8);
        cout << "Overall transition matrix:" << endl;
        data.communities[0].print_transition_matrix(8);
    } else
    {
        if(argc == 4)
        {
            cout << "Randomising some data for you..." << endl;
            ublas::matrix<double> test_t_m(3,5);
            test_t_m(0,0)=0.8;test_t_m(0,1)=0.1;test_t_m(0,2)=0.1;test_t_m(0,3)=0.0;test_t_m(0,4)=0.0;
            test_t_m(1,0)=0.1;test_t_m(1,1)=0.8;test_t_m(1,2)=0.1;test_t_m(1,3)=0.0;test_t_m(1,4)=0.0;
            test_t_m(2,0)=0.1;test_t_m(2,1)=0.1;test_t_m(2,2)=0.8;test_t_m(2,3)=0.0;test_t_m(2,4)=0.0;
            vector<string> test_names(3);
            test_names[0]="Will";test_names[1]="Fern";test_names[2]="Andre";
            Community data(test_t_m, test_names, atoi(argv[1]), atoi(argv[2]), "Test Community", atoi(argv[3]));
            data.optimise();
            cout << "Real transition matrix used to generate the data was:" << endl;
            data.print_real_transition_matrix(8);
            cout << "First community -> second community event matrix:" << endl;
            data.print_event_matrix(0,8);
            cout << "Overall transition matrix:" << endl;
            data.print_transition_matrix(8);
        }
        else
            cout << "Whoops! Please specify a file, or the number of communities, how many individuals you have to start, and a random seed for testing." << endl;
    }
                         
    return 0;
}

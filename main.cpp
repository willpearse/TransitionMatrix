//
//  main.cpp
//  transitionMatrix
//
//  Created by Will Pearse on 31/10/2011.
//  Copyright 2011 Imperial College London. All rights reserved.
//

#include <iostream>
#include <iomanip>
#include <fstream>
#include "communitySimulator.h"
#include "communityData.h"

using namespace std;
using namespace boost::numeric;

int main (int argc, const char * argv[])
{
    DataSet test("test.txt");
    //DataSet test("/Users/will/Library/Developer/Xcode/DerivedData/transitionMatrix-gxhgkoqjhuhbegctfntobfttbhgb/Build/Products/Debug/test.csv");
    test.communities[0].optimise();
    test.communities[0].print_community(0);
    test.communities[0].print_event_matrix(0);
    test.communities[0].print_transition_matrix();
    return 0;
}

//
//  main.cpp
//  scan1
//
//  Created by 孟令凯 on 2021/3/4.
//

#include <iostream>
#include "string"
#include <fstream>
#include "Graph.hpp"
using namespace std;

int main(int argc, const char * argv[]) {
    
    if(argc<5) return 0;
    // insert code here...
    clock_t startTime,endTime;
    startTime = clock();//计时开始
    
    std::cout << "Hello, World!\n";
//    string str = "/Users/milk/test_data/wiki-Talk.txt";
    string str = argv[1];
//    cout<<argv[0];
    string f = "/Users/milk/test_data/index/";
    
    
//    double a =0.4;//atof(argv[1]);
    Graph graph(str,f);
    graph.creatIndex(atof(argv[2]),atof(argv[3]),atoi(argv[4]));
//    graph.creatIndex(0.4,0.4,5);
    return 0;
}


//
//  Graph.hpp
//  scan1
//
//  Created by 孟令凯 on 2021/3/4.
//

#ifndef Graph_hpp
#define Graph_hpp
#include <iostream>
#include <stdio.h>
#include "string"
#include <ext/hash_map>
#include <vector>
#include <list>
#include <fstream>
#include <queue>

using namespace std;
using namespace __gnu_cxx;

class Graph{
    string str;//图数据
    string f;//索引数据
    int n,m,m2;
    hash_map<int,vector<int>> out_edges;
    hash_map<int,vector<int>> in_edges;
    
    vector<vector<vector<double>>> core_order1;
    vector<vector<vector<double>>> core_order2;
    
    vector<vector<vector<double>>> neighbor_order1;
    vector<vector<vector<double>>> neighbor_order2;
    
//    int *neiNum;
    vector<int> neiNum;
    
   
    
public:
    Graph(string str, string f);
    void creatIndex(double parameter1,double parameter2, int parameterNum);
    void readIndex();
    void buildNeighborOrder();
    double getSimi1(int u, int v);
    double getSimi2(int u, int v);
    void getNeiNum();
    bool BinarySearch(vector<int> a, int value, int n);
    void writeNOrder();
    void writeCOrder();
    void writeCOrder2(int parameterNum);
    void query(double parameter1,double parameter2, int parameterNum);
    void query2(double parameter1,double parameter2, int parameterNum);
    vector<int> getAllNei(double parameter1,double parameter2,int v);
    vector<int> getAllCore(double parameter1,double parameter2,int corei);
    vector<int> getcoreNum(vector<vector<double>> order1, vector<vector<double>> order2, double parameter1,double parameter2);
    
};

#endif /* Graph_hpp */


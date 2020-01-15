#include <iostream>
#include <fstream>
#include <string>
#include <vector>

using namespace std; 

// Read data from file
// Expected format: parameter = value, all spaces and tabs are ignored

void readDataFromFile(ifstream &file, string parameter, double &value);
void readDataFromFile(ifstream &file, string parameter, int &value);

// Read column vector from file
// Expected format:

//# comments
//value11 value21
//value12 value22
//value13 value23
//value14 value24
//value15 value25
//value16 value26
//value17 value27

// columnNumber starts at 0

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<double> &vec);
void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<int> &vec);

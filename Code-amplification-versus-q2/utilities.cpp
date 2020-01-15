#include <iostream>
#include <fstream>
#include <string>
#include <vector>

#include "utilities.h"

using namespace std; 

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Remove tabs and spaces from string

// Function used by readDataFromFile

void removeTabsAndSpaces(string &str)
{
	int N = str.size();
	string str_new = "";
	
	for (int i=0; i<N; i++)
		if (str[i]!=' ' && str[i]!='	')
			str_new = str_new + str[i];
	
	str = str_new;
}

////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Read data from file

// Expected format: parameter = value, all spaces and tabs are ignored

void readDataFromFile(ifstream &file, string parameter, double &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter && line[b]=='=') 
			value = stod(line.substr(b+1, a-b-1));
	}
}

void readDataFromFile(ifstream &file, string parameter, int &value)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	string line = "";
	while (getline(file, line))
	{
		removeTabsAndSpaces(line);
		
		int a = line.size();
		int b = parameter.size();
		
		if (line.substr(0,b) == parameter && line[b]=='=') 
			value = stoi(line.substr(b+1, a-b-1));
	}
}




////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////

// Read column vector

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<double> &vec)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	vec.empty();
	
	string line;
	
	while (getline(file, line))
	{
		char separator = ' ';
		char comment = '#';
		int indexStart = 0;
		int indexStop = line.size();
		int separatorCount = 0;
		
		if (line.size()>0) if(line[0]==comment) continue; // skip comments
		
		for (int i=0; i<line.size(); i++) 
		{
			if (line[i]==separator) 
			{
				separatorCount ++;
				if (separatorCount==columnNumber)   indexStart = i+1;
				if (separatorCount==columnNumber+1) indexStop = i-1;
			}
		}
		
		string value_str = line.substr(indexStart,indexStop-indexStart+1);
		vec.push_back(stod(value_str));
	}
}

void readColumnVectorFromFile(ifstream &file, int columnNumber, vector<int> &vec)
{
	file.clear();
	file.seekg(0, ios_base::beg);
	
	vec.empty();
	
	string line;
	
	while (getline(file, line))
	{
		char separator = ' ';
		char comment = '#';
		int indexStart = 0;
		int indexStop;
		int separatorCount = 0;
		
		if (line.size()>0) if(line[0]==comment) continue; // skip comments
		
		for (int i=0; i<line.size(); i++) 
		{
			if (line[i]==separator) 
			{
				separatorCount ++;
				if (separatorCount==columnNumber)   indexStart = i+1;
				if (separatorCount==columnNumber+1) indexStop = i-1;
			}
		}
		
		string value_str = line.substr(indexStart,indexStop);
		vec.push_back(stoi(value_str));
	}
}


////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////


#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>   // for exit(), srand(), rand()
#include <time.h>
#include<sstream>
#include<math.h>
#include <algorithm>
#include "matrix.h"
using namespace std;
void readHapFile(vector<string> &hapVec, string genoSetFile);
void readTagFile(vector<int> &tagVec, string tagFile);
Matrix hapVec2Matrix(vector<string> hapVec);
void signMatrix (Matrix &resultMatrix);
void signHapMatrix (Matrix &resultMatrix);
double differ(Matrix resultMatrix, Matrix predMatrix);
double differSquareSum(Matrix resultMatrix, Matrix predMatrix);
int differSquareSum4Columns(Matrix resultMatrix, Matrix predMatrix);
void write2TagFile(vector<int> tagVec, string tagFile);
bool isTag(const vector<int>& tagVec, int SNP);
Matrix hap2GenoMatrix(Matrix hapMatrix1, Matrix hapMatrix2);
vector<string> Matrix2HapVec(Matrix x, vector<string> oneVec); 
void write2GenoFile(vector<string> genoVec, string genoFile);

int main(int argc, const char *argv[])
{	//--- argv[1] for genoFile population
	//--- argv[2] for tag file
	//--- argv[3] for tag-restricted geno
	//--- argv[5] for result genoFile
	//--- argv[4] for version;
			//--1 genotype 
			//--2 haplotype 

cout <<"       **************************************" <<endl;
cout <<"       *      MLR-tagging Phasing Package   *" <<endl;
cout <<"       *          version 1.3               *" <<endl;
cout <<"       *  Tag Selection based on Prediction *" <<endl;
cout <<"       *                                    *" <<endl;
cout <<"       *      Alexander Zelikovsky          *" <<endl;
cout <<"       *         Jingwu He                  *" <<endl;
cout <<"       * Department of Computer Science     *" <<endl;
cout <<"       *     Georgia State University       *" <<endl;
cout <<"       *            March 2006              *" <<endl;
cout <<"       **************************************" <<endl;
cout << endl;


	vector<string> oneVec, sVec;
	vector<int> tagVec, newTagVec;
	readHapFile(oneVec, argv[3]);
	readHapFile(sVec, argv[1]);
	readTagFile(newTagVec, argv[2]);
	string version(argv[4]);
	int i,j, SNP, maxSNP;
	double max =0;	
	Matrix coffMatrix, predMatrix, sampleMatrix;
	Matrix tempMatrix, resultMatrix, tagMatrix, transMatrix, oneMatrix;
	double a;
	vector<double> accVec;
	sampleMatrix = hapVec2Matrix(sVec);
	Matrix inputOneMatrix = hapVec2Matrix(oneVec);
	oneMatrix = Matrix(-1, inputOneMatrix.m_nRows, sampleMatrix.m_nCols);
    for(j=0; j<inputOneMatrix.m_nRows; j++)
	for (i=0; i < newTagVec.size(); i++)
	{
		oneMatrix.m_pData[j][newTagVec.at(i)] = inputOneMatrix.m_pData[j][i];
	}
	tempMatrix = sampleMatrix;
	int rank;
	vector<int> basePosition, rowPosition;
	tempMatrix.RREF(rank, basePosition, rowPosition);

	if (newTagVec.size() > rank)
	{
		for (i=0; i < rank; i++)
			tagVec.push_back(newTagVec.at(i));
		tagMatrix = Matrix(-1, sampleMatrix.m_nRows, tagVec.size()+1);
		for (i=0; i < sampleMatrix.m_nRows; i++)
			for (j=0; j < tagVec.size(); j++)
		tagMatrix.m_pData[i][j+1] = sampleMatrix.m_pData[i][tagVec.at(j)]; 
		transMatrix = tagMatrix.GetTransposed();
		tempMatrix = transMatrix * tagMatrix;
		tempMatrix.Invert();
		tempMatrix = tempMatrix * transMatrix;	
		coffMatrix = tempMatrix * sampleMatrix;  		
	}
	else 
	{
		for (i=0; i < newTagVec.size(); i++)
		tagVec.push_back(newTagVec.at(i));
		tagMatrix = Matrix(-1, sampleMatrix.m_nRows, tagVec.size()+1);
		for (i=0; i < sampleMatrix.m_nRows; i++)
			for (j=0; j<tagVec.size(); j++)
				tagMatrix.m_pData[i][j+1] = sampleMatrix.m_pData[i][tagVec.at(j)]; 

		transMatrix = tagMatrix.GetTransposed();  
		tempMatrix = transMatrix * tagMatrix;
		tempMatrix.Invert();
		tempMatrix = tempMatrix * transMatrix;
		coffMatrix = tempMatrix * sampleMatrix;  
	}

	  

		tagMatrix = Matrix(-1, oneMatrix.m_nRows, tagVec.size()+1);
		for (i=0; i < oneMatrix.m_nRows; i++)
			for (j=0; j<tagVec.size(); j++)
		tagMatrix.m_pData[i][j+1] = oneMatrix.m_pData[i][tagVec.at(j)]; 
			resultMatrix = tagMatrix * coffMatrix;

			if (version.find("G") != -1)
				signMatrix(resultMatrix);
			else
				signHapMatrix(resultMatrix);

			for (i=0; i < oneMatrix.m_nRows; i++)
				for (j=0; j < newTagVec.size(); j++)
					resultMatrix.m_pData[i][newTagVec.at(j)] = inputOneMatrix[i][j];



	vector<string> resultVec = Matrix2HapVec(resultMatrix, oneVec);

	string inputS = sVec.at(0);
	string outputS = resultVec.at(0);
	for (j=0; j < newTagVec.size(); j++)
	{
				
		if (inputS.at(j) != outputS.at(newTagVec.at(j)))
		{
			cout << "error" << endl;
		}
				
	}
	
	write2GenoFile(resultVec, argv[5]);


	return 0;	
}

Matrix hap2GenoMatrix(Matrix hapMatrix1, Matrix hapMatrix2)
{
	Matrix genoMatrix = hapMatrix1;
	int j;
	for (j=0; j < hapMatrix1.m_nCols;j++)
	{
		if (hapMatrix1.m_pData[0][j] != hapMatrix2.m_pData[0][j])
			genoMatrix.m_pData[0][j] = 0;
	}


	return genoMatrix;
}


bool isTag(const vector<int>& tagVec, int SNP)
{
	int i;
	for (i=0; i < tagVec.size(); i++)
	{
		if (SNP == tagVec.at(i))
			return true;

	}
	return false;
}

void readHapFile(vector<string>& hapSetVec, string hapSetFile) {
    ifstream infile;
    infile.open(hapSetFile.c_str());
	string line;
    if (!infile) {
        cerr << "Unable to read " <<  hapSetFile << endl;
        exit(1);
    }
	
	while(getline(infile, line))
        hapSetVec.push_back(line);
    infile.close();

	vector<string>::iterator startIterator;
	startIterator = hapSetVec.begin();
	hapSetVec.erase(startIterator, startIterator + 3);

}

Matrix hapVec2Matrix(vector<string> hapVec) 
{    
	int i, j;
	Matrix x;
	x =  Matrix(0.0, hapVec.size(), hapVec.at(0).size());	

		for (i=0; i < hapVec.size(); i++)
		{
			for(j=0; j<hapVec.at(0).size(); j++)
			{
				if(hapVec.at(i).substr(j,1).compare("0") == 0)
					x.m_pData[i][j] = -1;
				else if (hapVec.at(i).substr(j,1).compare("1") == 0)
					x.m_pData[i][j] = 1;
				else if (hapVec.at(i).substr(j,1).compare("2") == 0)
					x.m_pData[i][j] = 0;
				else 
					x.m_pData[i][j] = atof(hapVec.at(i).substr(j,1).c_str());
			}
			
		}

    return x;
}


void readTagFile(vector<int> &tagVec, string tagFile)
{
    ifstream infile;
    infile.open(tagFile.c_str());
	string line;
    if (!infile) {
        cerr << "Unable to read " <<  tagFile << endl;
        exit(1);
    }
	int k=0;
	while(getline(infile, line))
	{  
		if (k>=3)
			tagVec.push_back(atoi(line.c_str()));
		k++;
	}
    infile.close();

}


void signMatrix (Matrix &resultMatrix)
{
	int i, j;
	for (i=0; i<resultMatrix.m_nRows; i++)
	{
		//cout << i << endl;
		for (j=0; j<resultMatrix.m_nCols; j++)
		{
			if (resultMatrix.m_pData[i][j] > 0.5)
				resultMatrix.m_pData[i][j] = 1;
			else if (resultMatrix.m_pData[i][j] <= -0.5)
				resultMatrix.m_pData[i][j] = -1;
			else
				resultMatrix.m_pData[i][j] = 0;
		}	
	}	

}

void signHapMatrix (Matrix &resultMatrix)
{
	int i, j;
	for (i=0; i<resultMatrix.m_nRows; i++)
	{
		//cout << i << endl;
		for (j=0; j<resultMatrix.m_nCols; j++)
		{
			if (resultMatrix.m_pData[i][j] > 0)
			{
				resultMatrix.m_pData[i][j] = 1;
			}
			else if (resultMatrix.m_pData[i][j] <= 0)
			{
				resultMatrix.m_pData[i][j] = -1;
			}
				
		}	
	}	

}

double differ(Matrix resultMatrix, Matrix predMatrix)
{
	int i, j;
	double e=0;
	for (i=0; i<resultMatrix.m_nRows; i++)
	{
		//cout << i << endl;
		for (j=0; j<resultMatrix.m_nCols; j++)
		{

				if (resultMatrix.m_pData[i][j] != predMatrix.m_pData[i][j])	
					e++;

				
		}
	}

	double z = (double)(resultMatrix.m_nRows*resultMatrix.m_nCols)/100.0;

	//cout << e << endl;

	//cout << z << endl;

	double a = 100 - e/z;
	return a;

}


double differSquareSum(Matrix resultMatrix, Matrix predMatrix)
{
	int i, j;
	double e=0;
	for (i=0; i<resultMatrix.m_nRows; i++)
	{
		//cout << i << endl;
		for (j=0; j<resultMatrix.m_nCols; j++)
		{
				if (resultMatrix.m_pData[i][j] != predMatrix.m_pData[i][j])
					
					e = e + (resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j])*(resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j]);
		}
	}

	double z = (double)(resultMatrix.m_nRows*resultMatrix.m_nCols)/100.0;

	//cout << e << endl;

	//cout << z << endl;

	double a = 100 - e/z;
	return a;

}


int differSquareSum4Columns(Matrix resultMatrix, Matrix predMatrix)
{
	int i, j;
	double e=0;
	int maxCol;
	double max =0;
	for (j=0; j<resultMatrix.m_nCols; j++)
	{
		e=0;
		for (i=0; i<resultMatrix.m_nRows; i++)
		{
				if (resultMatrix.m_pData[i][j] != predMatrix.m_pData[i][j])				
					e = e + (resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j])*(resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j]);
		}
		if (e > max)
		{
			max = e;
			maxCol = j;
		}
	}

	if (max == 0)
		return -1;
	return maxCol;

}


void write2TagFile(vector<int> tagVec, string tagFile) {

	int i, j;
	ofstream ostr;
	ostr.open(tagFile.c_str());
	ostr << tagVec.size() << endl;
	ostr << 1 << endl;
	ostr << "result tagFile " << endl;
	for(i=0; i<tagVec.size(); i++) {
		ostr << tagVec.at(i) << endl;
	}	
	ostr.close();
}

vector<string> Matrix2HapVec(Matrix x, vector<string> oneVec) 
{   
	oneVec.clear();
	int i, j;
	string tempString="";
	for(i=0; i<x.m_nRows ; i++){
	   for(j=0; j<x.m_nCols; j++)
	   {
	      if(x.m_pData[i][j] == -1)
					tempString = tempString+'0';
				else if (x.m_pData[i][j] == 1)
					tempString = tempString+'1';
				else 
					tempString =tempString+'2';
	   }
	   oneVec.push_back(tempString);
	   tempString.clear();
	   } 
	/*oneVec.clear();
	string tempString (x.m_nCols,'0'); 
	oneVec.push_back(tempString);
	int i, j;
		for (i=0; i < oneVec.size(); i++)
		{
			for(j=0; j<oneVec.at(0).size(); j++)
			{
				if(x.m_pData[i][j] == -1)
					oneVec.at(i).at(j) = '0';
				else if (x.m_pData[i][j] == 1)
					oneVec.at(i).at(j) = '1';
				else 
					oneVec.at(i).at(j) ='2';
			}
		}*/

    return oneVec;
}

void write2GenoFile(vector<string> genoVec, string genoFile) {

	int i, j;
	ofstream ostr;
	ostr.open(genoFile.c_str());
	ostr << genoVec.size() << endl;
	ostr << genoVec.at(0).size() << endl;
	ostr << "result individual " << endl;
	for(i=0; i<genoVec.size(); i++) {
		ostr << genoVec.at(i) << endl;
	}	
	ostr.close();
}


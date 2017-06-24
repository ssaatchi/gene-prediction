
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
bool sameHap(string a, string b);
Matrix hap2GenoMatrix(Matrix hapMatrix1, Matrix hapMatrix2);
double calculateR(vector<string> twoSNPVec);
double getNewRsq(int col,vector<string> PredVec, vector<string> PredTagsVec, Matrix resultMatrix );
vector<string> form2SNPVec(int j, vector<string> hapVec, vector<string> predHapVec);
double CalHapFreq(Matrix resultMatrix, int i);
//double coveredSNPs(Matrix resultMatrix, Matrix predMatrix, double r2);
vector<int> coveredSNPs(Matrix resultMatrix, Matrix predMatrix, double r2);
vector<string> Matrix2HapVec(Matrix x);
double averageR2(Matrix resultMatrix, Matrix predMatrix);
string itos(int i);
string ftos(double i);
//int gl=0;

int main(int argc, const char *argv[])
{	
	//--- argv[1] for sample genoFile
	//--- argv[2] for R2
	//--- argv[4] for genotype or haplotype

	clock_t first=clock();
	string version(argv[4]);	      	   
	vector<string> oneVec, sVec, PredTagsVec, PredVec;
	vector<int> tagVec, newTagVec;
	readHapFile(sVec, argv[1]);
	int numberOfTags;
	double r2= atof(argv[2]);
	int i, j, SNP, maxSNP;
	double max =0, rsq=0;	
	Matrix coffMatrix, predMatrix, sampleMatrix, PredTags;
	Matrix tempMatrix, resultMatrix, tagMatrix, transMatrix, oneMatrix;
	double a;
	vector<double> accVec;
	sampleMatrix = hapVec2Matrix(sVec);
	//cout << sVec.at(0).size() << endl;
	newTagVec.clear();
	double prev=0; 
	vector<int> v, v1;
	
   
	
	int temp;

	int limit;
	
	v.clear();
	cout<<"1\n";
	//************ larger than rank ***********
	tempMatrix = sampleMatrix;
	cout<<"2\n";
	int rank;
	cout<<"3\n";
	vector<int> basePosition, rowPosition;
	cout<<"4\n";
	//tempMatrix.RREF(rank, basePosition, rowPosition);
	//cout << rank << endl;
	
	
	
	bool flag = false;
	
	int ku;
    while (true)
	{
		cout<<"goes here1\n";
		tagVec.clear();	
		
		if (flag)
		{
		
		  
		  cout<<"come here\n";
			for (i=0; i < limit; i++)
			{
				tagVec.push_back(newTagVec.at(i));
			}

			tagMatrix = Matrix(-1, sampleMatrix.m_nRows, tagVec.size());

			for (i=0; i < sampleMatrix.m_nRows; i++)
				for (j=0; j < tagVec.size(); j++)
				{
					tagMatrix.m_pData[i][j] = sampleMatrix.m_pData[i][tagVec.at(j)-1]; 
				}

			transMatrix = tagMatrix.GetTransposed();
			tempMatrix = transMatrix * tagMatrix;
			tempMatrix.Invert();
			tempMatrix = tempMatrix * transMatrix;	
			coffMatrix = tempMatrix * sampleMatrix;  
			resultMatrix = tagMatrix * coffMatrix;
			for (i=0; i < sampleMatrix.m_nRows; i++)
				for (j=0; j < newTagVec.size(); j++)
				{
					resultMatrix.m_pData[i][tagVec.at(j)] = sampleMatrix[i][tagVec.at(j)];
				}
			maxSNP = differSquareSum4Columns(sampleMatrix, resultMatrix);
			if (maxSNP == -1)
			{	
				for (SNP=0; SNP<sampleMatrix.m_nCols; SNP++)
				{
					if (!isTag(newTagVec, SNP))
					{
						maxSNP = SNP;
					}
				}
			}

			newTagVec.push_back(maxSNP);
		//cout << "Tag #" << newTagVec.size() << ":" << newTagVec.at(newTagVec.size()-1) << endl;
		}
		else 
		{
		    cout<<"goes here2\n";
			for (i=0; i < newTagVec.size(); i++)
			{
				tagVec.push_back(newTagVec.at(i));
			}
			tagVec.push_back(0);
			//cout<<"test1 pass\n";
			for (SNP=0; SNP<sampleMatrix.m_nCols; SNP++)
			{
				accVec.push_back(0);
				if (!isTag(v, SNP))
				{			
				   cout<<"SNP: "<<SNP<<"  goes here3\n ";		
					tagVec.at(tagVec.size()-1) = SNP;
					
					cout<<"tagVec : ";
					for(ku=0; ku<tagVec.size(); ku++)
					   cout<<tagVec.at(ku)<<", ";
					cout<<"\n";
					
					
					tagMatrix = Matrix(-1, sampleMatrix.m_nRows, tagVec.size());
					for (i=0; i < sampleMatrix.m_nRows; i++)
						for (j=0; j<tagVec.size(); j++)
						{
							tagMatrix.m_pData[i][j] = sampleMatrix.m_pData[i][tagVec.at(j)]; 
						}
						
						
					/*	cout<<"tagmatrix\n";
						tagMatrix.Display();
						*/
						
						//cout<<"test2 pass\n";
					transMatrix = tagMatrix.GetTransposed(); 
					tempMatrix = transMatrix * tagMatrix;
					tempMatrix.Invert();
					tempMatrix = tempMatrix * transMatrix;
					coffMatrix = tempMatrix * sampleMatrix;  
					resultMatrix = tagMatrix * coffMatrix;
					
					//cout<<"test3 pass\n";
					
					if (version.find("G") != -1)
					{
				
						signMatrix(resultMatrix);
					}
					else
					{
						signHapMatrix(resultMatrix);
					}
					for (i=0; i < resultMatrix.m_nRows; i++)
						for (j=0; j<tagVec.size(); j++)
						{
							resultMatrix.m_pData[i][tagVec.at(j)] = sampleMatrix.m_pData[i][tagVec.at(j)]; 
						}
					v =  coveredSNPs(sampleMatrix, resultMatrix, r2);
				/*	cout << "    SnP TAg: " << SNP << " covered "<< v.size() << " number of SNPs as follows:" << endl;
			for (int k=0; k < v.size(); k++)
				cout << v.at(k) << ",";
					*/
					if (v.size() > max)
					{
						maxSNP = SNP;
						max = v.size();
					}
				}
			}
			//cout << maxSNP << endl;
			if ((newTagVec.size () >= 1) && (maxSNP == tagVec.at(newTagVec.size()-1)))
			{
				cout<<"goes here4\n";
			for (SNP=0; SNP<sampleMatrix.m_nCols; SNP++)
			{
				if ((!isTag(v, SNP))&&(!isTag(newTagVec, SNP))){
					newTagVec.push_back(SNP);
					v.push_back(SNP);
						cout << "tag: " << SNP << " covered "<< v.size() << " number of SNPs as follows:" << endl;
			for (int k=0; k < v.size(); k++)
				cout << v.at(k) << ",";
			cout << endl;
				}
			
			}
			
					
			break;
				
			}
			else{
			cout<<"goes here5\n";
			if (!isTag(newTagVec, maxSNP))
				newTagVec.push_back(maxSNP);
			
			for (i=0; i < sampleMatrix.m_nRows; i++)
			for (j=0; j<newTagVec.size(); j++)
			{
				tagMatrix.m_pData[i][j] = sampleMatrix.m_pData[i][newTagVec.at(j)]; 
			}
			transMatrix = tagMatrix.GetTransposed(); 
			tempMatrix = transMatrix * tagMatrix;
			tempMatrix.Invert();
			tempMatrix = tempMatrix * transMatrix;
			coffMatrix = tempMatrix * sampleMatrix;  
			resultMatrix = tagMatrix * coffMatrix;
				if (version.find("G") != -1)
					{
				
						signMatrix(resultMatrix);
					}
					else
					{
						signHapMatrix(resultMatrix);
					}
			for (i=0; i < resultMatrix.m_nRows; i++)
				for (j=0; j<newTagVec.size(); j++)
			{
				resultMatrix.m_pData[i][newTagVec.at(j)] = sampleMatrix.m_pData[i][newTagVec.at(j)]; 
			}
			
			PredTagsVec.clear();
			PredVec.clear();
			string tmp=""; 
			
			cout<< "rows :"<<tagMatrix.m_nRows<<"Cols: "<<tagMatrix.m_nCols<<"\n";
			
			for (int i2=0; i2 < tagMatrix.m_nRows; i2++)
			{	tmp="";
			   for (int j2=0; j2<tagMatrix.m_nCols; j2++)
				{	
				     tmp = tmp + ftos(tagMatrix.m_pData[i2][j2]);	
				}
				
				PredTagsVec.push_back(tmp);

			}	
			
			//cout<<"loop1 covered\n";
			
			for (int l=0; l< resultMatrix.m_nCols; l++){
			PredVec.clear();
			   for (int l3=0; l3< resultMatrix.m_nRows; l3++){
					PredVec.push_back(ftos(resultMatrix.m_pData[l3][l]));
					//cout<<"r["<<l3<<"]["<<l<<"]   "<<resultMatrix.m_pData[l3][l]<<"\n";
			   }
			 //  cout<<"in 2 "<<l<<"\n";
			  //rsq=getNewRsq(l, PredVec, PredTagsVec, resultMatrix);  
			}
			   
			
			v =  coveredSNPs(sampleMatrix, resultMatrix, r2);
			cout << "tag: " << maxSNP << " covered "<< v.size() << " number of SNPs as follows:" << endl;
			for (int k=0; k < v.size(); k++)
				cout << v.at(k) << ",";
			cout << endl;
			temp = v.size();
			}

		}
		if (v.size() == sampleMatrix.m_nCols)
			break;
			
	}
	
	cout<<"tagmat \n";
	/*
	
	tagMatrix.Display();
	
	cout<<"predtagsvec \n";
	
	for(int k=0; k<PredTagsVec.size(); k++)
	 cout<<PredTagsVec.at(k)<<"\n";
	*/ 
	
	//v =  coveredSNPs(sampleMatrix, resultMatrix, r2);

	//cout<<"v size = "<<v.size();
	
	
	
	
	write2TagFile(newTagVec, argv[3]);
	
	clock_t second=clock();
	
	
	cout<<"Time taken = "<<second-first<<"\n";
	
	
	return 0;	
}

double getNewRsq(int col,vector<string> PredVec, vector<string> PredTagsVec, Matrix resultMatrix )
{
   int i,j;
   double aa=0,ab=0,bb=0,ba=0, rsq=0;
   	cout<<"inside...1"<<"\n";
   
    for (i = 0; i < PredTagsVec.size(); i++){
                aa=0;
                ab=0;
                bb=0;
                ba=0;
                	cout<<"inside...2"<<"\n";
                for (j = 0; j < PredTagsVec.size(); j++){
                    if (PredVec.at(j) == PredVec.at(0)){
                       if(!sameHap(PredTagsVec.at(i), PredTagsVec.at(j))){
                           cout<<"no loop...\n";
                          /*   ab += CalHapFreq(resultMatrix,j);
                        }else{
                            aa += CalHapFreq(resultMatrix,j);
                        }
                    }else{
                        if(!sameHap(PredTagsVec.at(i), PredTagsVec.at(j))){
                            bb += CalHapFreq(resultMatrix,j);
                        }else{
                            ba += CalHapFreq(resultMatrix,j);*/
                        }
                    }
                }
                //p is snp's freq, q is hap's freq
                	cout<<"inside...3"<<"\n";
               double p = aa+ab;
                double q = ba+aa;
                //round to 5 decimal places.
                rsq = ((aa*bb - ab*ba) * (aa*bb - ab*ba) ) / (p*(1-p)*q*(1-q)) ;
               /* if (rsq > curBestRsq){
                    StringBuffer sb = new StringBuffer();
                    for (int j = 1; j < genos[i].length; j++){
                        sb.append(genos[i][j]);
                    }
                    curBestAllele = new Allele(theBlock,sb.toString());
                    curBestRsq = rsq;*/
                }
          // } 
          //}//del this
            
            cout<<"aa = "<<aa<<"\n";
            cout<<"ab = "<<ab<<"\n";
            cout<<"bb = "<<bb<<"\n";
            cout<<"ba = "<<ba<<"\n";
            
            cout<<"rsq "<<rsq<<"\n";
            return rsq;

}

double CalHapFreq(Matrix resultMatrix, int i){

   int count=0, total;
   total = resultMatrix.m_nRows;
   string resultM="", tmp="";
   double percent;
   //cout<<"  resultMatrix.m_pData : "<<resultMatrix.m_pData[i]<<"\n";
   
   for(int k=0; k<resultMatrix.m_nCols; k++)
      resultM=resultM + ftos(resultMatrix.m_pData[i][k]);
   //cout<<"i:"<<i<<"  result: "<<resultM<<"\n";
   
   for(int j=0; j<resultMatrix.m_nRows; j++){
      tmp="";
     if(j!=i){
      for(int j2=0; j2<resultMatrix.m_nCols; j2++){
             tmp = tmp + ftos(resultMatrix.m_pData[j][j2]);	
		}
		//cout<<"tmp "<<tmp<<"\n";
		if(sameHap(tmp,resultM))
		    count++;}
	}
	
	//cout<<"count "<<count<<"\n";
   percent = (double) count/total;
    
    
    //cout<<"per "<<percent<<"\n";
    return percent;

   }

 bool sameHap(string a, string b) {
       
       //cout<<"a "<<a<<"   b  "<<b<<"\n";
        if(a.length() != b.length()) {
            return false;
        }
        for(int i=1;i<a.size();i++) {
            if(a[i] != b[i]) {
                return false;
            }
        }
        
        return true;
 }
 

vector<int> coveredSNPs(Matrix resultMatrix, Matrix predMatrix, double r2)
{
	vector<int> v;
	int i, j;
	double r;
	int maxCol;
	double max =0;
	double k=1;
	vector<string> twoSNPVec;
	vector<string> predVec =  Matrix2HapVec(predMatrix);
	vector<string> resultVec = Matrix2HapVec(resultMatrix);
	for (j=0; j<resultMatrix.m_nCols; j++)
	{
			twoSNPVec = form2SNPVec(j, predVec, resultVec);
			r = calculateR(twoSNPVec);
			//cout<<j+1<<") "<<r<<"\n";
			if ( r>=r2){
				v.push_back(j);
			}
	
	}
	//cout << v.size() << endl;
	return v;
}




Matrix hap2GenoMatrix(Matrix hapMatrix1, Matrix hapMatrix2)
{
	Matrix genoMatrix = hapMatrix1;
	int j;
	for (j=0; j < hapMatrix1.m_nCols;j++)
	{
		if (hapMatrix1.m_pData[0][j] != hapMatrix2.m_pData[0][j])
		{
			genoMatrix.m_pData[0][j] = 0;
		}
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
	{  
        hapSetVec.push_back(line);
    }
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
		{
			tagVec.push_back(atoi(line.c_str()));
		}
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
			{
				resultMatrix.m_pData[i][j] = 1;
			}
			else if (resultMatrix.m_pData[i][j] <= -0.5)
			{
				resultMatrix.m_pData[i][j] = -1;
			}
			else
			{
				resultMatrix.m_pData[i][j] = 0;
			}
		
				
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
				{
					
					e++;

				}
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
				{
					
					e = e + (resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j])*(resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j]);

				}
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
				{				
					e = e + (resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j])*(resultMatrix.m_pData[i][j] - predMatrix.m_pData[i][j]);
				}
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

/*
double coveredSNPs(Matrix resultMatrix, Matrix predMatrix, double r2)
{
	int i, j;
	double r;
	int maxCol;
	double max =0;
	double k=1;
	vector<string> twoSNPVec;
	vector<string> predVec =  Matrix2HapVec(predMatrix);
	vector<string> resultVec = Matrix2HapVec(resultMatrix);
	for (j=0; j<resultMatrix.m_nCols; j++)
	{
			twoSNPVec = form2SNPVec(j, predVec, resultVec);
			r = calculateR(twoSNPVec);
			if ( r>=r2)
				k++;
	
	}
	return k;
}

*/


double averageR2(Matrix resultMatrix, Matrix predMatrix)
{
	int i, j;
	double r;
	int maxCol;
	double max =0;
	double k=0;
	vector<string> twoSNPVec;
	vector<string> predVec =  Matrix2HapVec(predMatrix);
	vector<string> resultVec = Matrix2HapVec(resultMatrix);
	for (j=0; j<resultMatrix.m_nCols; j++)
	{
			twoSNPVec = form2SNPVec(j, predVec, resultVec);
			k = k+ calculateR(twoSNPVec);
	}
	k= k/(double)resultMatrix.m_nCols;
	return k;
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

double calculateR(vector<string> twoSNPVec)
{
	int n = twoSNPVec.size();
	int i;
	int p00=0, p01=0, p11=0, p10=0, p02=0, p20=0, p21=0, p12=0, p22=0;
	int h00=0, h01=0, h10=0, h11=0; 
	for (i=0; i<n; i++)
	{
		if(twoSNPVec.at(i) == "00")
			p00++;
		else if(twoSNPVec.at(i) == "01")
			p01++;
		else if(twoSNPVec.at(i) == "10")
			p10++;
		else if(twoSNPVec.at(i) == "11")
			p11++;
			
		else if(twoSNPVec.at(i) == "02")
			p02++;
		else if(twoSNPVec.at(i) == "20")
		    p20++;
		else if(twoSNPVec.at(i) == "21")
			p21++;
		else if(twoSNPVec.at(i) == "12")
		    p12++;
		    
		else if(twoSNPVec.at(i) == "22")
		    p22++;
	}
	
	//cout<<"22 vals : "<<p20<<"    "<<p02<<"     "<<p21<<"    "<<p12<<"    "<<p22<<"\n";
	
	h00 = (2 * p00) + p02 + p20 + p22;
	h01 = (2 * p01) + p02 + p21;
	h10 = (2 * p01) + p20 + p12;
	h11 = (2 * p11) + p12 + p21 + p22;

	double r = (double)((h11*h00 - h10*h01) * (h11*h00 - h10*h01)) /(double)((h00+h01)*(h10+h11)*(h00+h10)*(h01+h11));
	
	if (((h00+h01)*(h10+h11)*(h00+h10)*(h01+h11)) == 0)
		r =0;

    //gl++;
    //cout<<gl<<") r= "<<r<<"\n";
    
	return r;
}

vector<string> form2SNPVec(int j, vector<string> hapVec, vector<string> predHapVec)
{

	vector<string> twoSNPVec;
	int i;
	for (i = 0; i < hapVec.size(); i++)
	{
		twoSNPVec.push_back("00");
	}

	for (i = 0; i < twoSNPVec.size(); i++)
	{
		if (hapVec.at(i).at(j) == '1')
			twoSNPVec.at(i).at(0) = '1';
	    else if (hapVec.at(i).at(j) == '2')
	        twoSNPVec.at(i).at(0) = '2';
	        
		if (predHapVec.at(i).at(j) == '1')
			twoSNPVec.at(i).at(1) = '1';
		else if (predHapVec.at(i).at(j) == '2')
			twoSNPVec.at(i).at(1) = '2';
	}
	return twoSNPVec;
}

/*double calculateR(vector<string> twoSNPVec)
{
	int n = twoSNPVec.size();
	int i;
	int p00=0, p01=0, p11=0, p10=0; 
	for (i=0; i<n; i++)
	{
		if(twoSNPVec.at(i) == "00")
			p00++;
		else if(twoSNPVec.at(i) == "01")
			p01++;
		else if(twoSNPVec.at(i) == "10")
			p10++;
		else if(twoSNPVec.at(i) == "11")
			p11++;
	}

	double r = (double)((p11*p00 - p10*p01) * (p11*p00 - p10*p01)) /(double)((p00+p01)*(p10+p11)*(p00+p10)*(p01+p11));
	
	if (((p00+p01)*(p10+p11)*(p00+p10)*(p01+p11)) == 0)
		r =0;

    //gl++;
    //cout<<gl<<") r= "<<r<<"\n";
    
	return r;
}


vector<string> form2SNPVec(int j, vector<string> hapVec, vector<string> predHapVec)
{

	vector<string> twoSNPVec;
	int i;
	for (i = 0; i < hapVec.size(); i++)
	{
		twoSNPVec.push_back("00");
	}

	for (i = 0; i < twoSNPVec.size(); i++)
	{
		if (hapVec.at(i).at(j) == '1')
			twoSNPVec.at(i).at(0) = '1';
		if (predHapVec.at(i).at(j) == '1')
			twoSNPVec.at(i).at(1) = '1';
	}
	return twoSNPVec;
}*/


vector<string> Matrix2HapVec(Matrix x) 
{    
	vector<string> resultVec;
	string tempString (x.m_nCols,'0'); 
	int i, j;
		for (i=0; i < x.m_nRows; i++)
		{
			for(j=0; j<x.m_nCols; j++)
			{
				if(x.m_pData[i][j] == -1)
					tempString.at(j) = '0';
				else if (x.m_pData[i][j] == 1)
					tempString.at(j) = '1';
				else 
					tempString.at(j) ='2';
			}
			resultVec.push_back(tempString);
		}

    return resultVec;
}

//----------------------------------------------------
string itos(int i) // convert int to string
{
  stringstream s;
  s << i;
  return s.str();
}

string ftos(double i) // convert float to string
{
  stringstream s;
  s << i;
  return s.str();
}

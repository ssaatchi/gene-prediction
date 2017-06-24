#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <stdlib.h>
#include <time.h>
#include <sstream>
#include <math.h>

using namespace std;
void write2File(vector<string> data, string baseFile);
string itos(int i);

int main(int argc, char *argv[])
{

  //argv[1] - name of file
  //argv[2] - number of parts. As far as possible try to make it divisible
  
  string Sample = argv[1];
  int kfold=atoi(argv[2]);
  vector<string> data, testdata, partdata;
  string  GenoLine, tagLine, Logline,TTLine, tag2Line;
  
  int l1, l2, i;
  
  
  ifstream sampleip;
  ofstream testop,sampop;
  
 // cout<<"reach 2\n";
  
  sampleip.open(Sample.c_str());
  /*for( int i=0; i<parts; i++){
  sampop.open("sample"+i+".txt", ios::app);
  testop.open("test"+i+".txt", ios::app);
  }*/
  
 // cout<<"reach 3\n";
  
  if (!sampleip) {
        cerr << "Unable to read " <<  Sample << endl;
        exit(1);
    }
    
        
  for(i=0; i<3;i++){
	   getline(sampleip,GenoLine);
	 }
   
   //cout<<"reach 4\n";
   
   /*while(!tag2ip.eof())
   {
   //cout<<"reach 5\n";
        getline(tag2ip, tag2Line);
        int p=atoi(tag2Line.c_str());
        tag2vals.push_back(p);
      //  cout<<"reach 7\n";
   }*/
   
    
    int flag=0;
    int k;    
  while(!sampleip.eof())
  {
		getline(sampleip,GenoLine);
		 data.push_back(GenoLine);
  }
  
  int nolines = data.size();
  nolines--;
  
  for(k=1; k<=kfold; k++)
  {
    int start = (k-1)*nolines/kfold;
    int end = k*nolines/kfold;
    
    cout<<"start = "<<start<<"   end = "<<end<<endl;
      
    for (i=start; i<end; i++)
         testdata.push_back(data.at(i));
    if (start!=0)
       for (i=0; i<start; i++)
         partdata.push_back(data.at(i));
    if (end!=nolines)
       for (i=end; i<nolines; i++)
         partdata.push_back(data.at(i));
     
      write2File(partdata,"part"+itos(k)+".txt");
     
      write2File(testdata,"test"+itos(k)+".txt");
      
      partdata.clear();
      testdata.clear();
   }
 
/*for(k=1;k<=kfold;k++)       
{
for(i=1;i<=2;i++)
      {
 if(i==k)    
     write2File(data,"part"+itos(k)+".txt",(i-1)*nolines/kfold,(i*nolines)/kfold);
      else  
         write2File(data,"test"+itos(k)+".txt",(i-1)*nolines/kfold,(i*nolines)/kfold);
       }
  
}*/
}


void write2File(vector<string> data, string baseFile) {
		ofstream ostr;
		ostr.open(baseFile.c_str());
		
		ostr << data.size() << endl;
	ostr << 1 << endl;
	ostr << "result File " << endl;
		for(int i=0; i< data.size(); i++) {
		ostr << data.at(i) << endl;
		}
		ostr.close();
}

string itos(int i) // convert int to string
{
  stringstream s;
  s << i;
  return s.str();
}

/**
 * Miscellaneous Functions Definitions
 * Supvisor: Dr. Alex Zelikovsky
 * <PRE>
 * misc.cc
 *
 * Revisions:   1.0  April 13, 2004
 * misc
 *
 * </PRE>
 * Author: Jingwu He (PH.D candidate)
 * Address: Computer Science Department, Georgia State University
 * @author <A HREF="mailto:jingwu_he@yahoo.com">Jingwu He</A>
 * @version 1.0, April 13, 2004
 *
 */
#include <iostream>
#include <string>
#include "matrix.h"
using namespace std;

//display error message, exit program if necessary
void ErrorMsg(string str, bool mustexit)
{
	cout << "Error: " << str << '\n';

	if(mustexit) {
		cout << "Exiting program\n";
		exit(1);
	}
	else {
		char ch = 'Q';

		do {
			cout << "Continue program execution? (Y/N): ";
			cin >> ch;
			ch = toupper(ch);
		} while((ch != 'Y') && (ch != 'N'));

		if(ch == 'N') exit(1);
		if(ch == 'Y') return;        //yes, superfluous, but makes it more undertandable
	}
}

//get safe int input
int getint(istream& istr)
{
	int temp = 0;
	string str;

	while(!(istr >> temp)) {
		istr.clear();
		istr >> str;
	}

	return temp;
}


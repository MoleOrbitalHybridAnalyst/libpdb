#include "pdb.h"
#include "utili.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <regex>

#include <iostream>
#include <cstdio>

using namespace std;
using namespace PDB_NS;

PDB::PDB(const string& fname) 
{
   boxlens[0]=boxlens[1]=boxlens[2]=0.0;
   ifstream fs(fname);
   string line;
   size_t line_count = 0;
   bool incomplete = false;
   while(getline(fs, line)) {
      line_count++;
      if(line.substr(0,4)=="ATOM") {
         if(line.size()<minLineSize) {
            cerr << "WARNING: Incomplete line at line "<<line_count<<" in "
               << fname << endl;
            incomplete = true;
         }
         linenumbers.push_back(line_count);
         //TODO atomnames, resnames, ...

         string tmpstr;
         tmpstr = "    ";
         readField(line.substr(12,4),tmpstr);
         atomnames.push_back(tmpstr);
         tmpstr = "    ";
         readField(line.substr(17,4),tmpstr);
         resnames.push_back(tmpstr);
         tmpstr = " ";
         readField(line.substr(21,1),tmpstr);
         chainids.push_back(tmpstr[0]);

         if(incomplete) {
            incomplete = false;
         }

         if(line_count>214079 or line_count<5) {
            cout << atomnames.back() << endl;
            cout << resnames.back() << endl;
            cout << chainids.back() << endl;
            //cout << chainids.back() << endl;
            //cout << ss.str() << endl;
            //cout << tmpstr << endl;
         }
      }
      if(line.substr(0,5)=="CRYST") {
         size_t start = line.find_first_of(" \t\r");
         if(start!=string::npos) {
            //istrstream line_stream(line.substr(start).c_str());
            stringstream ss(line.substr(start));
            ss >> boxlens[0] >> boxlens[1] >> boxlens[2];
         }
          printf("%f %f %f\n",boxlens[0],boxlens[1],boxlens[2]);
         pair<size_t,string> linepair(line_count,line);
         nonatomlines.push_back(linepair);
         continue;
      }
      if(line.substr(0,3)=="END") {
         pair<size_t,string> linepair(line_count,line);
         nonatomlines.push_back(linepair);
         break;
      } else { //not CRYST && not END && not ATOM
         pair<size_t,string> linepair(line_count,line);
         nonatomlines.push_back(linepair);
      }
   }
}

void PDB::eraseSpace(string& str) {
   str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
}


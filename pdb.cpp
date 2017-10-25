#include "pdb.h"
#include "utili.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <exception>

#include <iostream>
#include <cstdio>
#include <regex>

using namespace std;
using namespace PDB_NS;

PDB::PDB(const string& fname) 
{
   boxlens[0]=boxlens[1]=boxlens[2]=0.0;
   ifstream fs(fname);
   string line;
   size_t line_count = 0;
   //bool incomplete = false;
   while(getline(fs, line)) {
      line_count++;
      if(line.substr(0,4)=="ATOM") {
         //if(line.size()<minLineSize) {
         //   // TODO move this err to a check function 
         //   //cerr << "WARNING: Incomplete line at line "<<line_count<<" in "
         //   //   << fname << endl;
         //   incomplete = true;
         //}
         //if(line.size()<errLineSize) {
         //    //throw an err
         //}
         linenumbers.push_back(line_count);
         line.resize(pdbLineSize,' ');
         //TODO atomnames, resnames, ...
         vector<bool> defined;

         string tmpstr;
         tmpstr = " ";
         defined.push_back(readField(line.substr(12,4),tmpstr));
         //} catch(std::out_of_range& exception) {
         //    atomnames.push_back(" "); resnames.push_back(" ");
         //    chainids.push_back(' '); resids.push_back(-1);
         //    xs.push_back(0.0); ys.push_back(0.0); zs.push_back(0.0);
         //    occ.push_back(1.0); tempfs.push_back(0.0);
         //    segnames.push(" "); atomtypes.push(" ");
         //    vector<bool> defined(false,11);
         //    defineds.push_back(defined); continue;
         //}
         atomnames.push_back(tmpstr);
         tmpstr = " ";
         defined.push_back(readField(line.substr(17,4),tmpstr));
         resnames.push_back(tmpstr);
         tmpstr = " ";
         defined.push_back(readField(line.substr(21,1),tmpstr));
         chainids.push_back(tmpstr[0]);

         int resid=-1;
         //if(!readField(line.substr(22,4),resid)) {
         //    cerr << "WARNING: No resid at line "<<line_count<<" in "
         //        << fname << endl;
         //}
         defined.push_back(readField(line.substr(22,4),resid));
         resids.push_back(resid);

         float tmpflt;
         tmpflt = 0.0;
         defined.push_back(readField(line.substr(30,8),tmpflt));
         xs.push_back(tmpflt);
         tmpflt = 0.0;
         defined.push_back(readField(line.substr(38,8),tmpflt));
         ys.push_back(tmpflt);
         tmpflt = 0.0;
         defined.push_back(readField(line.substr(46,8),tmpflt));
         zs.push_back(tmpflt);
         tmpflt = 1.0;
         defined.push_back(readField(line.substr(54,6),tmpflt));
         occs.push_back(tmpflt);
         tmpflt = 0.0;
         defined.push_back(readField(line.substr(60,6),tmpflt));
         tempfs.push_back(tmpflt);

         tmpstr = " ";
         defined.push_back(readField(line.substr(72,4),tmpstr));
         resnames.push_back(tmpstr);
         tmpstr = " ";
         defined.push_back(readField(line.substr(76,2),tmpstr));
         atomtypes.push_back(tmpstr);

         if(line_count>214079 or line_count<5) {
            cout << atomnames.back() << endl;
            cout << resnames.back() << endl;
            cout << chainids.back() << endl;
            //cout << chainids.back() << endl;
            //cout << ss.str() << endl;
            //cout << tmpstr << endl;
         }
         defineds.push_back(defined);
         continue;
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
   fs.close();
}

void write2file(const string& fname) const {
}

//void PDB::eraseSpace(string& str) {
//   str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
//}


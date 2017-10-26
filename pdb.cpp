#include "pdb.h"
#include "utili.h"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <utility>
#include <exception>
#include <cstdio>
#include <iostream>
#include <cstdlib>


#include <regex>
#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace PDB_NS;

//static constexpr size_t atomname = PDBField::atomname;
//static constexpr size_t resname = PDBField::resname;
//static constexpr size_t segname = PDBField::segname;
//static constexpr size_t atomtype = PDBField::atomtype;
//static constexpr size_t chainid = PDBField::chainid;
//static constexpr size_t resid = PDBField::resid;
//static constexpr size_t x = PDBField::x;
//static constexpr size_t y = PDBField::y;
//static constexpr size_t z = PDBField::z;
//static constexpr size_t occ = PDBField::occ;
//static constexpr size_t tempf = PDBField::tempf;

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
         // atomnames, resnames, ...
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
         segnames.push_back(tmpstr);
         tmpstr = " ";
         defined.push_back(readField(line.substr(76,2),tmpstr));
         atomtypes.push_back(tmpstr);

         //@@@ if(line_count>214079 or line_count<5) {
         //@@@    cout << atomnames.back() << endl;
         //@@@    cout << resnames.back() << endl;
         //@@@    cout << chainids.back() << endl;
         //@@@    //cout << chainids.back() << endl;
         //@@@    //cout << ss.str() << endl;
         //@@@    //cout << tmpstr << endl;
         //@@@ }
         defineds.push_back(defined);
         continue;
      }
      if(line.substr(0,5)=="CRYST") {
         size_t start = line.find_first_of(" \t\r\n\f\v");
         if(start!=string::npos) {
            //istrstream line_stream(line.substr(start).c_str());
            //auto t1 = high_resolution_clock::now();
            stringstream ss(line.substr(start));
            ss >> boxlens[0] >> boxlens[1] >> boxlens[2];
            //auto t2 = high_resolution_clock::now();
            //cout << duration_cast<duration<double>>(t2-t1).count() << endl;
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
   nAtoms = atomnames.size();
}

void PDB::centerAlignedPrint4(FILE *fp, const string& s) const 
{
   string tmpstr;
   switch(s.size()) {
      case 1: tmpstr = " " + s + "  "; // _A__
              break;
      case 2: tmpstr = " " + s + " ";  // _AB_
              break;
      case 3: tmpstr = " " + s;        // _ABC
              break;
      case 4: tmpstr = s;
              break;
   }
   fprintf(fp, "%s", tmpstr.c_str());
}

void PDB::write2file(const string& fname) const 
{
   FILE *fp = fopen(fname.c_str(),"w");
   if(!fp) {
      cerr << "libpdb internal error: cannot open "<<fname<<'\n';
      abort();
   }
   auto iter = nonatomlines.begin();
   for(size_t index = 0; index < atomnames.size(); ++index) {
      //wait for nonatomlines
      size_t linenumber = linenumbers[index];
      while(iter->first < linenumber) {
         fprintf(fp, "%s\n", (iter->second).c_str());
         iter++;
      }

      //writting ATOM
      fprintf(fp, "ATOM  ");
      if(index >= 99998) {
         fprintf(fp, "***** ");
      } else {
         fprintf(fp, "%5d ", int(index) + 1);
      }
      centerAlignedPrint4(fp, atomnames[index]);
      fprintf(fp," %-4s%c%4d    %8.3f%8.3f%8.3f%6.2f%6.2f",
            resnames[index].c_str(), chainids[index],
            resids[index], xs[index], ys[index], zs[index],
            occs[index], tempfs[index]);
      fprintf(fp,"      %-4s%2s\n", 
            segnames[index].c_str(), atomtypes[index].c_str());
   }

   //write the left nonatomlines
   while(iter != nonatomlines.end()) {
      fprintf(fp, "%s\n", (iter->second).c_str());
      iter++;
   }
   fclose(fp);
}

vector<pair<size_t,string>> PDB::checkUndefined() const 
{
   vector<pair<size_t,string>> results;
   for(auto iter = defineds.begin(); iter != defineds.end(); ++iter) {
      for(auto iter2 = iter->begin(); iter2 != iter->end(); ++iter2) {
         if(*iter2 == false) {
            pair<size_t,string> result;
            result.first = iter - defineds.begin();
            result.second = 
               transField(static_cast<PDBField>(iter2 - iter->begin()));
            results.push_back(result);
         }
      }
   }
   return results;
}

string PDB::transField(const PDBField& pdbfield) const 
{
   switch(pdbfield) {
      case PDBField::atomname: return "atomname";
      case PDBField::resname: return "resname";
      case PDBField::segname: return "segname";
      case atomtype: return "atomtype";
      case chainid: return "chainid";
      case resid: return "resid";
      case x: return "x";
      case y: return "y";
      case z: return "z";
      case occ: return "occ";
      case tempf: return "tempf";
      default: return "\0";
   }
}

bool PDB::guessOneChainid(const size_t index) 
{
   if(defineds[index][chainid]) {
      return true;
   }
   if(defineds[index][PDBField::segname]) {
      chainids[index] = segnames[index][0];
      defineds[index][chainid] = true;
      return true;
   } 
   return false;
}

bool PDB::guessAllChainids() 
{
   bool f = true;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(!guessOneChainid(index)) f = false;
   }
   return f;
}

bool PDB::guessOneSegname(const size_t index) 
{
   if(defineds[index][PDBField::segname]) {
      return true;
   }
   if(defineds[index][chainid]) {
      segnames[index] = string(1,chainids[index]);
      defineds[index][PDBField::segname] = true;
      return true;
   } 
   return false;
}

bool PDB::guessAllSegnames() 
{
   bool f = true;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(!guessOneSegname(index)) f = false;
   }
   return f;
}

bool PDB::guessOneAtomtype(const size_t index) 
{
   if(defineds[index][atomtype]) {
      return true;
   }
   if(defineds[index][PDBField::atomname]) {
      string atomname = atomnames[index];
      if(atomname == "CLA") {
         atomtypes[index]="CL";
      } else if(atomname == "POT") {
         atomtypes[index]="K";
      } else if(atomname == "SOD") {
         atomtypes[index]="NA";
      } else {
         atomtypes[index] = string(1,atomname[0]);
      }
      defineds[index][chainid] = true;
      return true;
   } 
   return false;
}

bool PDB::guessAllAtomtypes() 
{
   bool f = true;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(!guessOneAtomtype(index)) f = false;
   }
   return f;
}

void PDB::swapFields(const size_t i1, const size_t i2, 
      const vector<PDBField>& fields)
{
   if(i1 >= nAtoms and i2 >= nAtoms) {
      cerr << "libpdb error: index exceeds natom when swapping\n";
      abort();
   }
   for(auto iter = fields.begin(); iter != fields.end(); ++iter) {
      // swap and also defineds!!
      switch(*iter) {
         case x: swap(xs[i1], xs[i2]); 
                 break;
         case y: swap(ys[i1], ys[i2]); 
                 break;
         case z: swap(zs[i1], zs[i2]);
                 break;
         case resid: swap(resids[i1], resids[i2]);
                 break;
         case PDBField::atomname: swap(atomnames[i1],atomnames[i2]);
                 break;
         case PDBField::resname: swap(resnames[i1], resnames[i2]);
                 break;
         case PDBField::segname: swap(segnames[i1], segnames[i2]);
                 break;
         case chainid: swap(chainids[i1], chainids[i2]);
                 break;
         case atomtype: swap(atomtypes[i1], atomtypes[i2]);
                 break;
         case occ: swap(occs[i1], occs[i2]);
                 break;
         case tempf: swap(tempfs[i1], tempfs[i2]);
         default: break;
      }
      swap(defineds[i1][*iter], defineds[i2][*iter]);
   }
}

void PDB::swapFields(const size_t i1, const size_t i2, 
      const PDBField& field)
{
   swapFields(i1, i2, vector<PDBField>(1,field));
}

void PDB::swapCoordinates(const size_t i1, const size_t i2)
{
   PDBField xyz[3] = {x, y, z};
   swapFields(i1, i2, vector<PDBField>(xyz, xyz + 3));
}

size_t PDB::reorderWater(bool guess, bool check, bool reorder,
      const PDBDef& defo, const PDBDef& defh)
{
   if(guess) {
      guessAllChainids();
      guessAllSegnames();
      guessAllAtomtypes();
   }
   if(check) {
      //TODO check whether all the things involving def defined
   }
   if(reorder) { //reorder H and O
      //TODO reorder
   }
   //TODO real part of reorder
   return 0;
}

size_t PDB::reorderWater(const PDBDef& defo, const PDBDef& defh)
{
   return reorderWater(false, false, false, defo, defh);
}



//void PDB::eraseSpace(string& str) {
//   str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
//}


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

//vector<pair<size_t,string>> PDB::checkUndefined() const 
vector<pair<size_t,PDBField>> PDB::checkUndefined() const 
{
   vector<pair<size_t,PDBField>> results;
   for(auto iter = defineds.begin(); iter != defineds.end(); ++iter) {
      for(auto iter2 = iter->begin(); iter2 != iter->end(); ++iter2) {
         if(*iter2 == false) {
            pair<size_t,PDBField> result;
            result.first = iter - defineds.begin();
            result.second = 
               static_cast<PDBField>(iter2 - iter->begin());
               //transField(static_cast<PDBField>(iter2 - iter->begin()));
            results.push_back(result);
         }
      }
   }
   return results;
}

//string PDB::transField(const PDBField& pdbfield) const 
//{
//   switch(pdbfield) {
//      case PDBField::atomname: return "atomname";
//      case PDBField::resname: return "resname";
//      case PDBField::segname: return "segname";
//      case PDBField::atomtype: return "atomtype";
//      case PDBField::chainid: return "chainid";
//      case PDBField::resid: return "resid";
//      case PDBField::x: return "x";
//      case PDBField::y: return "y";
//      case PDBField::z: return "z";
//      case PDBField::occ: return "occ";
//      case PDBField::tempf: return "tempf";
//      default: return "\0";
//   }
//}

bool PDB::guessOneChainid(const size_t index) 
{
   if(defineds[index][static_cast<size_t>(PDBField::chainid)]) {
      return true;
   }
   if(defineds[index][static_cast<size_t>(PDBField::segname)]) {
      chainids[index] = segnames[index][0];
      defineds[index][static_cast<size_t>(PDBField::chainid)] = true;
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
   if(defineds[index][static_cast<size_t>(PDBField::segname)]) {
      return true;
   }
   if(defineds[index][static_cast<size_t>(PDBField::chainid)]) {
      segnames[index] = string(1,chainids[index]);
      defineds[index][static_cast<size_t>(PDBField::segname)] = true;
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
   if(defineds[index][static_cast<size_t>(PDBField::atomtype)]) {
      return true;
   }
   if(defineds[index][static_cast<size_t>(PDBField::atomname)]) {
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
      defineds[index][static_cast<size_t>(PDBField::atomtype)] = true;
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
         case PDBField::x: std::swap(xs[i1], xs[i2]); 
                 break;
         case PDBField::y: std::swap(ys[i1], ys[i2]); 
                 break;
         case PDBField::z: std::swap(zs[i1], zs[i2]);
                 break;
         case PDBField::resid: std::swap(resids[i1], resids[i2]);
                 break;
         case PDBField::atomname: std::swap(atomnames[i1],atomnames[i2]);
                 break;
         case PDBField::resname: std::swap(resnames[i1], resnames[i2]);
                 break;
         case PDBField::segname: std::swap(segnames[i1], segnames[i2]);
                 break;
         case PDBField::chainid: std::swap(chainids[i1], chainids[i2]);
                 break;
         case PDBField::atomtype: std::swap(atomtypes[i1], atomtypes[i2]);
                 break;
         case PDBField::occ: std::swap(occs[i1], occs[i2]);
                 break;
         case PDBField::tempf: std::swap(tempfs[i1], tempfs[i2]);
         default: break;
      }
      std::swap(defineds[i1][static_cast<size_t>(*iter)], 
                           defineds[i2][static_cast<size_t>(*iter)]);
   }
}

void PDB::swapFields(const size_t i1, const size_t i2, 
      const PDBField& field)
{
   swapFields(i1, i2, vector<PDBField>(1,field));
}

void PDB::swapCoordinates(const size_t i1, const size_t i2)
{
   PDBField xyz[3] = {PDBField::x, PDBField::y, PDBField::z};
   swapFields(i1, i2, vector<PDBField>(xyz, xyz + 3));
}

size_t PDB::reorderWater(bool guess, bool check, bool reorder,
      const PDBDef& defo, const PDBDef& defh, const PDBDef& defhyd)
{
   if(guess) {
      guessAllChainids();
      guessAllSegnames();
      guessAllAtomtypes();
   }
   if(check) {
      for(size_t index = 0; index < nAtoms; ++index) {
         auto defstr = defo.getDefstr();
         for(auto iter = defstr.begin(); iter != defstr.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         defstr = defh.getDefstr();
         for(auto iter = defstr.begin(); iter != defstr.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         auto defflt = defo.getDefflt();
         for(auto iter = defflt.begin(); iter != defflt.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         defflt = defh.getDefflt();
         for(auto iter = defflt.begin(); iter != defflt.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         auto defint = defh.getDefint();
         for(auto iter = defint.begin(); iter != defint.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         defint = defo.getDefint();
         for(auto iter = defint.begin(); iter != defint.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         auto defchr = defh.getDefchr();
         for(auto iter = defchr.begin(); iter != defchr.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
         defchr = defo.getDefchr();
         for(auto iter = defchr.begin(); iter != defchr.end(); ++iter) {
            if(!defineds[index][static_cast<size_t>(iter->first)]) {
               cerr<< "libpdb error: undefined "<< transField(iter->first)<<
                  " of atom " << index + 1<<'\n'; abort();
            }
         }
      }
   }
   //Lv1 index vector gives real atom index by dereference once
   vector<size_t> oindexesLv1, hindexesLv1;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(isMatched(index, defo)) oindexesLv1.push_back(index);
      if(isMatched(index, defh)) hindexesLv1.push_back(index);
   }
   size_t hydindex; int count_hyd=0;
   for(auto iter = oindexesLv1.begin(); iter != oindexesLv1.end(); ++iter) {
      if(isMatched(*iter, defhyd)) {
         if(count_hyd > 1) {
            cerr << "libpdb error: more than one hydronium\n"; abort();
         }
         hydindex = *iter; count_hyd++;
      }
   }
   if(!count_hyd) {
      cerr << "libpdb error: no hydronium found\n"; abort();
   }
   //cout << oindexesLv1.size() << ' ' << hindexesLv1.size() << endl;
   if(hindexesLv1.size() != 2*oindexesLv1.size() + 1) {
      cerr << "libpdb error: number of hydrogens and "<<
                     "number of oxygens do not match\n";
      abort();
   }
   if(reorder) { 
      //TODO reorder H and O so that they are in the order of OHHOHH... 
      //and throw H and O to the bottom of the pdb
      //defined allfields for swaping two atoms completely
      PDBField tmparr[11] = {PDBField::atomname,PDBField::resname,
         PDBField::segname,PDBField::atomtype,PDBField::chainid,PDBField::resid,
      PDBField::x,PDBField::y,PDBField::z,PDBField::occ,PDBField::tempf};
      vector<PDBField> allfields(tmparr, tmparr + 11);
      cerr << "libpdb warning: reorder not implememted yet\n";
   }
   //now every O has two following H's
   vector<size_t> hindexesLv1bck = hindexesLv1;
   for(auto io = oindexesLv1.begin(); io != oindexesLv1.end(); ++io) {
      //countj* are Lv2 indexes
      size_t countj1 = size_tMax, countj2 = size_tMax; countj = 0;
      for(auto jh = hindexesLv1.begin(); jh != hindexesLv1.end(); ++jh) {
         float d1, d2, d2_;
         if(*jh == size_tMax) continue;
         if(countj1 == size_tMax) {
            d1 = pbcDistance2(*io, *jh);
            countj1 = countj;
         } else if(countj2 == size_tMax) {
            d2 = pbcDistance2(*io, *jh);
            countj2 = countj;
            if(d1 > d2) {
               std::swap(d1, d2);
               std::swap(countj1, countj2);
            }
         } else {
            d2_ = pbcDistance2(*io, *jh);
            if(d2_ < d2 and d2_ > d1) {
               d2 = d2_; countj2 = countj;
            } else if(d2_ < d1 and d2_ > d2) {
               d1 = d2_; countj1 = countj;
            } else if(d2_ < d1 and d2_ < d2) {
               if(d1 < d2) {
                  d2 = d2_; countj2 = countj;
               } else {
                  d1 = d2_; countj1 = countj;
               }
            }
         }
         int flag = 0;
         for(size_t countj = 0; countj < hindexesLv1bck.size(); countj++) {
            if(hindexesLv1bck[countj] == *io + 1) {
               flag++; hindexesLv1[countj] = size_tMax;
            }
            if(hindexesLv1bck[countj] == *io + 2) {
               flag++; hindexesLv1[countj] = size_tMax;
            }
         }
         if(flag != 2) {
            cerr << "libpdb error: cannot find two hydrogens of atom "
               <<*io+1 << '\n';
            abort();
         }
         swapCoordinates(*io+1, hindexesLv1bck[countj1]);
         swapCoordinates(*io+2, hindexesLv1bck[countj2]);
         countj++;
      }
   }

   return hydindex;
}

size_t PDB::reorderWater(
      const PDBDef& defo, const PDBDef& defh, const PDBDef& defhyd)
{
   return reorderWater(false, false, false, defo, defh, defhyd);
}

array<float,3> PDB::pbcDistance(const array<float,3>& x1,
                           const array<float,3>& x2) const 
{
   array<float,3> diff;
   for(int i=0; i < 3; ++i) diff[i] = pbcDiff(x1[i], x2[i], i);
   return diff;
}

array<float,3> PDB::pbcDistance(size_t i1, size_t i2) const 
{
   array<float,3> x1, x2;
   x1[0] = xs[i1]; x1[1] = ys[i1]; x1[2] = zs[i1];
   x2[0] = xs[i2]; x2[1] = ys[i2]; x2[2] = zs[i2];
   return pbcDistance(x1, x2);
}

bool PDB::isMatched(size_t index, const PDBDef& def) const
{
   if(!isMatched(atomnames[index], def, PDBField::atomname)) return false;
   if(!isMatched(resnames[index], def, PDBField::resname)) return false;
   if(!isMatched(segnames[index], def, PDBField::segname)) return false;
   if(!isMatched(atomtypes[index], def, PDBField::atomtype)) return false;
   if(!isMatched(chainids[index], def, PDBField::chainid)) return false;
   return true;
}

//bool PDB::isMatched(size_t index, const PDBdef& def) const
//{
//   if(index > nAtoms) 
//      return false;
//   bool result = false;
//   while((auto search = def.find(PDBField::atomname)) != def.end()) {
//      if(search->second == atomnames[index]) {
//         result = true; break;
//      }
//   }
//   if(!result) return false;
//   result = false;
//   while((auto search = def.find(PDBField::resname)) != def.end()) {
//      if(search->second == resnames[index]) {
//         result = true; break;
//      }
//   }
//   if(!result) return false;
//   result = false;
//   while((auto search = def.find(PDBField::segname)) != def.end()) {
//      if(search->second == segnames[index]) {
//         result = true; break;
//      }
//   }
//   if(!result) return false;
//   result = false;
//   while((auto search = def.find(PDBField::segname)) != def.end()) {
//      if(search->second == segnames[index]) {
//         result = true; break;
//      }
//   }
//
//   return true;
//
//}


//void PDB::eraseSpace(string& str) {
//   str.erase(remove_if(str.begin(),str.end(),::isspace),str.end());
//}


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


//#include <regex>
//#include <chrono>
//using namespace std::chrono;

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
         //try {
            defined.push_back(readField(line.substr(22,4),resid));
         //}
         //catch(invalid_argument &) {
         //   cerr << "libpdb warning: cannot read resid at line "
         //      << line_count << " in " << fname << endl;
         //   defined.push_back(false);
         //}
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
         //@@@  printf("%f %f %f\n",boxlens[0],boxlens[1],boxlens[2]);
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
      cerr << "libpdb error: index "
         << i1 << ' ' << i2 << " exceeds natom when swapping\n";
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

void PDB::swapFields(const size_t i1, const size_t i2)
{
   PDBField fields[11] = { PDBField::atomname, PDBField::resname, PDBField::segname, PDBField::atomtype, PDBField::chainid, PDBField::resid, PDBField::x, PDBField::y, PDBField::z, PDBField::occ, PDBField::tempf};
   vector<PDBField> fields_vec(fields, fields + 11);
   swapFields(i1, i2, fields_vec);
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
   if(reorder) { 
      //assemble H and O so that they are in the order of OHHOHH... 
      //and throw H and O to the bottom of the pdb
      //defined allfields for swaping two atoms completely
      //PDBField tmparr[11] = {PDBField::atomname,PDBField::resname,
      //   PDBField::segname,PDBField::atomtype,PDBField::chainid,PDBField::resid,
      //PDBField::x,PDBField::y,PDBField::z,PDBField::occ,PDBField::tempf};
      //vector<PDBField> allfields(tmparr, tmparr + 11);
      //cerr << "libpdb warning: assemble water not implememted yet\n";
      if(!assembleWater(false, false, defo, defh)) {
         cerr << "libpdb error: assemble water failed\n"; abort();
      }
   }
   //now every O has two following H's 
   //Lv1 index vector gives real atom index by dereference once
   vector<size_t> oindexesLv1, hindexesLv1;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(isMatched(index, defo)) oindexesLv1.push_back(index);
      if(isMatched(index, defh)) hindexesLv1.push_back(index);
   }
   size_t hydindex = 0; int count_hyd=0;
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
   //modified version of my implementation in python
   //vector<size_t> hindexesLv1bck(hindexesLv1);
   auto starth = hindexesLv1.begin();
   for(auto io = oindexesLv1.begin(); io != oindexesLv1.end(); ++io) {
      //j*Lv2 are Lv2 indexes if defined as j*Lv2[i] = i
      starth += 2;
      size_t j1Lv2, j2Lv2, jLv2 = 0;
      float d1 = pbcDistance2(*io, *io + 1);
      float d2 = pbcDistance2(*io, *io + 2);
      j1Lv2 = starth - hindexesLv1.begin() - 2;
      j2Lv2 = j1Lv2 + 1;
      jLv2 = j2Lv2 + 1;
      if(d1 > d2) {
         std::swap(d1, d2);
         std::swap(j1Lv2, j2Lv2);
      }

      for(auto jh = starth; jh != hindexesLv1.end(); ++jh) {
         float d = pbcDistance2(*io, *jh);
         if(d < d2 and d > d1) {
            d2 = d; j2Lv2 = jLv2;
         } else if(d < d1 and d > d2) {
            d1 = d; j1Lv2 = jLv2;
         } else if(d < d1 and d < d2) {
            if(d1 < d2) {
               d2 = d; j2Lv2 = jLv2;
            } else {
               d1 = d; j1Lv2 = jLv2;
            }
         }
         jLv2++;
      }
      swapCoordinates(*io+1, hindexesLv1[j1Lv2]);
      swapCoordinates(*io+2, hindexesLv1[j2Lv2]);
   }

// slow because of frequent swapping
//   auto starth = hindexesLv1.begin() + 2;
//   for(auto io = oindexesLv1.begin(); io != oindexesLv1.end(); ++io) {
//      float d1 = pbcDistance2(*io, *io + 1);
//      float d2 = pbcDistance2(*io, *io + 2);
//      if(d1 > d2) swapCoordinates(*io + 1, *io + 2);
//      for(auto jh = starth; jh != hindexesLv1.end(); ++jh) {
//         float d = pbcDistance2(*io, *jh);
//         if(d < d2) swapCoordinates(*io + 2, *jh);
//         if(d1 > d2) swapCoordinates(*io + 1, *io + 2);
//      }
//      starth += 2;
//   }

   float mindist2 = 
      boxlens[0]*boxlens[0] + boxlens[1]*boxlens[1] + boxlens[2]*boxlens[2];
   for(auto io = oindexesLv1.begin(); io != oindexesLv1.end(); ++io) {
      //double d = pbcDistance2(*io, hindexesLv1bck.back());
      double d = pbcDistance2(*io, hindexesLv1.back());
      if(d < mindist2) {
         mindist2 = d; hydindex = *io;
      }
   }
   swapCoordinates(hydindex, oindexesLv1.back());
   swapCoordinates(hydindex + 1, oindexesLv1.back() + 1);
   swapCoordinates(hydindex + 2, oindexesLv1.back() + 2);
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

bool PDB::moveTo(const size_t i1, const size_t i2) 
{
   int step;
   if(i1 == i2) {
      return true;
   } else if(i1 > i2) {
      step = -1;
   } else {
      step = 1;
   }
   if(i1 >= nAtoms or i2 >= nAtoms) return false;
   for(size_t walker = i1; walker != i2; walker += step) 
      swapFields(walker, walker + step);
   return true;
}

bool PDB::assembleWater(bool guess, bool check,
      const PDBDef& defo, const PDBDef& defh)//, const PDBDef& defhyd)
{
   if(guess) {
      guessAllChainids();
      guessAllSegnames();
      guessAllAtomtypes();
   }
   //duplicated in reorderWater becaue I am a so noob programmer
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
   // real assembling stuff
   vector<size_t> oindexesLv1, hindexesLv1;
   for(size_t index = 0; index < nAtoms; ++index) {
      if(isMatched(index, defh)) hindexesLv1.push_back(index);
   }
   // first put all H to the bottom
   // start from end so that hindexes is always valid
   size_t pos = nAtoms - 1;
   for(auto ith = hindexesLv1.rbegin(); ith != hindexesLv1.rend(); ++ith) {
         moveTo(*ith, pos); pos--;
   } // hindexes is not valid any more
   for(size_t index = 0; index < nAtoms; ++index) {
      if(isMatched(index, defo)) oindexesLv1.push_back(index);
   }
   if(hindexesLv1.size() != 2*oindexesLv1.size() + 1) {
      cerr << "libpdb error: number of hydrogens and "<<
                     "number of oxygens do not match\n";
      abort();
   }
   for(auto ito = oindexesLv1.rbegin(); ito != oindexesLv1.rend(); ++ito) {
      moveTo(*ito, pos); pos--;
   } // oindexes is not valid any more
   // switch some O and H
   pos = nAtoms - 2 * oindexesLv1.size() - 2; 
   for(size_t dest = nAtoms - 4;pos > nAtoms - 3 * oindexesLv1.size() - 1; dest -= 3) {
      swapFields(pos, dest); pos--;
   }
   pos = nAtoms - 2 * oindexesLv1.size() - 4;
   for(size_t offset = 0; pos + offset < nAtoms - 4; offset++) {
      size_t index = pos + offset;
      if(offset % 3 == 0) {
         if(!defineds[index][static_cast<size_t>(PDBField::atomname)]) 
            setAtomname(index, "OH2");
         setResid(index, offset / 3 + 1);
      } else if(offset % 3 == 1) {
         if(!defineds[index][static_cast<size_t>(PDBField::atomname)]) 
            setAtomname(index, "H1");
         setResid(index, offset / 3 + 1);
      } else if(offset % 3 == 2) {
         if(!defineds[index][static_cast<size_t>(PDBField::atomname)]) 
            setAtomname(index, "H2");
         setResid(index, offset / 3 + 1);
      }
   }
   if(!defineds[nAtoms - 4][static_cast<size_t>(PDBField::atomname)]) 
      setAtomname(nAtoms - 4, "OH2");
   if(!defineds[nAtoms - 3][static_cast<size_t>(PDBField::atomname)]) 
      setAtomname(nAtoms - 3, "H1");
   if(!defineds[nAtoms - 2][static_cast<size_t>(PDBField::atomname)]) 
      setAtomname(nAtoms - 2, "H2");
   if(!defineds[nAtoms - 1][static_cast<size_t>(PDBField::atomname)]) 
      setAtomname(nAtoms - 1, "H3");
   for(size_t index = nAtoms - 4; index < nAtoms; ++index) {
      setResid(index, 1);
      setResname(index, "H3O");
   }
   //auto ito = oindexesLv1.begin();
   //size_t opos = nAtoms - 4, hpos = nAtoms - 1;
   //moveToWithIndexes(*ito, opos, oindexes, hindexes); ito ++; opos -= 3;
   //moveToWithIndexes(*ith, hpos, oindexes, hindexes); ith ++; hpos --;
   //moveToWithIndexes(*ith, hpos, oindexes, hindexes); ith ++; hpos --;
   //moveToWithIndexes(*ith, hpos, oindexes, hindexes); ith ++; hpos -= 2;
   //for(;ito != oindexesLv1.end(); ++ito)  {
   //   moveToWithIndexes(*ito, opos, oindexes, hindexes); ito++; opos -= 3;
   //}
   //for(;ith != hindexesLv1.end(); ith += 2)  {
   //   moveToWithIndexes(*ith, hpos, oindexes, hindexes); ith++; hpos--;
   //   moveToWithIndexes(*(ith + 1), hpos, oindexes, hindexes); 
   //   ith++; hpos--; hpos-;
   //}
   return true;
}

bool PDB::moveToWithIndexes(const size_t i1, const size_t i2,
      vector<size_t>& list1, vector<size_t>& list2)
{
   if(!moveTo(i1, i2)) return false;
   if(i1 == i2) return true;
   //find element right larger/smaller than i1 in list1 and list2
   auto start1 = list1.end() - 1;
   auto start2 = list2.end() - 1;
   if(i1 < i2) {
      for(auto it = list1.begin(); it != list1.end(); ++it) 
         if(*it > i1) { start1 = it; break;}
      for(auto it = list2.begin(); it != list2.end(); ++it) 
         if(*it > i1) { start2 = it; break;}
   } else {
      start1 = list1.begin();
      start2 = list2.begin();
      for(auto it = list1.rbegin(); it != list1.rend(); ++it) 
         if(*it < i1) { start1 = it.base() - 1; break;}
      for(auto it = list2.rbegin(); it != list2.rend(); ++it) 
         if(*it < i1) { start2 = it.base() - 1; break;}
   }
   //find element right smaller/larger than i2 in list1 and list2
   auto end1 = list1.begin();
   auto end2 = list2.begin();
   if(i1 < i2) {
      for(auto it = list1.rbegin(); it != list1.rend(); ++it) 
         if(*it <= i2) { end1 = it.base() - 1; break;}
      for(auto it = list2.rbegin(); it != list2.rend(); ++it) 
         if(*it <= i2) { end2 = it.base() - 1; break;}
   } else {
      end1 = list1.end() - 1;
      end2 = list2.end() - 1;
      for(auto it = list1.begin(); it != list1.end(); ++it) 
         if(*it >= i2) { end1 = it; break;}
      for(auto it = list2.begin(); it != list2.end(); ++it) 
         if(*it >= i2) { end2 = it; break;}
   }
   int step;
   if(i1 < i2) step = 1;
   else        step = -1;
   for(auto it = start1; it != end1; it += step)  (*it) -= step;
   (*end1) -= step;
   for(auto it = start2; it != end2; it += step)  (*it) -= step;
   (*end2) -= step;
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


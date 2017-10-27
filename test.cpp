#include "pdb.h"
#include <iostream>
#include <array>
#include <chrono>
using namespace std::chrono;

using namespace std;
using namespace PDB_NS;

int main(int argc, char **argv) {
   auto t1 = high_resolution_clock::now();
   PDB pdb(argv[1]);
   auto t2 = high_resolution_clock::now();
   cout << "reading costs "<<duration_cast<duration<double>>(t2-t1).count()<<endl;
   pdb.write2file("write.pdb");
   t1 = high_resolution_clock::now();
   cout << "writing write.pdb costs "<<duration_cast<duration<double>>(t1-t2).count()<<endl;
   auto undefineds = pdb.checkUndefined();
   t2 = high_resolution_clock::now();
   cout << "check undefined costs "<<duration_cast<duration<double>>(t2-t1).count()<<endl;
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << transField(iter->second) << " of atom " << iter->first+1
         << "\n";
   }
   cout << "************\n";

   pdb.guessAllChainids();
   undefineds = pdb.checkUndefined();
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << transField(iter->second) << " of atom " << iter->first+1
         << "\n";
   }

   cout << "************\n";
   t1 = high_resolution_clock::now();
   pdb.guessAllSegnames();
   pdb.guessAllAtomtypes();
   t2 = high_resolution_clock::now();
   cout << "geussing segnames atomtypes costs "<<duration_cast<duration<double>>(t2-t1).count()<<endl;
   undefineds = pdb.checkUndefined();
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << transField(iter->second) << " of atom " << iter->first+1
         << "\n";
   }
   pdb.write2file("write2.pdb");

   double x = 3.1415926;
   vector<PDBField> fields(1,PDBField::x);
   pdb.swapFields(214082,214083,fields);
   //pdb.swapCoordinates(214082,214083);
   pdb.write2file("write3.pdb");
   
   //pdbdef
   cout << "************\n";
   PDBDef pdbdef("def1");
   pdbdef.pushBack(PDBField::x, -1.0f);
   pdbdef.pushBack(PDBField::resid, 908);
   pdbdef.pushBack(PDBField::atomname, "OH2");
   auto defstr = pdbdef.getDefstr();
   for(auto iter=defstr.begin();iter!=defstr.end();++iter) {
      cout << "select "<< transField(iter->first)<< " as "<<iter->second<<'\n';
   }
   auto defflt = pdbdef.getDefflt();
   for(auto iter=defflt.begin();iter!=defflt.end();++iter) {
      cout << "select "<< transField(iter->first)<< " as "<<iter->second<<'\n';
   }
   auto defint = pdbdef.getDefint();
   for(auto iter=defint.begin();iter!=defint.end();++iter) {
      cout << "select "<< transField(iter->first)<< " as "<<iter->second<<'\n';
   }
   auto defchr = pdbdef.getDefchr();
   for(auto iter=defchr.begin();iter!=defchr.end();++iter) {
      cout << "select "<< transField(iter->first)<< " as "<<iter->second<<'\n';
   }

   //pbc
   cout << "************\n";
   array<float,3> x1 = {1,2,3}, x2 = {1000, 50, -1000}, x3 = {90, -90, 5000};
   auto d1 = pdb.pbcDistance(x1, x2);
   auto d2 = pdb.pbcDistance(x1, x3);
   cout << d1[0] << ' ' << d1[1] << ' ' << d1[2] << '\n';
   cout << d2[0] << ' ' << d2[1] << ' ' << d2[2] << '\n';

   //reorder
   pdb.setSegname(214073,"W17");
   //pdb.setChainid(214073,'W');
   t1 = high_resolution_clock::now();
   PDBDef defo, defh, defhyd;
   defo.pushBack(PDBField::chainid, 'W');defo.pushBack(PDBField::atomtype, "O");
   defh.pushBack(PDBField::chainid, 'W');defh.pushBack(PDBField::atomtype, "H");
   defhyd = defo;
   defhyd.pushBack(PDBField::resname, "H3O");
   pdb.reorderWater(true, true, true, defo, defh, defhyd);

   t2 = high_resolution_clock::now();
   cout << "reordering costs "<<duration_cast<duration<double>>(t2-t1).count()<<endl;





   return 0;
}

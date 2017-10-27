#include "pdb.h"
#include <iostream>

using namespace std;
using namespace PDB_NS;

int main(int argc, char **argv) {
   PDB pdb(argv[1]);
   pdb.write2file("write.pdb");
   auto undefineds = pdb.checkUndefined();
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << iter->second << " of atom " << iter->first+1
         << "\n";
   }
   cout << "************\n";

   pdb.guessAllChainids();
   undefineds = pdb.checkUndefined();
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << iter->second << " of atom " << iter->first+1
         << "\n";
   }

   cout << "************\n";
   pdb.guessAllSegnames();
   pdb.guessAllAtomtypes();
   undefineds = pdb.checkUndefined();
   for(auto iter=undefineds.begin();iter!=undefineds.end();iter++) {
      cout << "undefined " << iter->second << " of atom " << iter->first+1
         << "\n";
   }
   pdb.write2file("write2.pdb");

   double x = 3.1415926;
   vector<PDBField> fields(1,PDBField::x);
   pdb.swapFields(214082,214083,fields);
   //pdb.swapCoordinates(214082,214083);
   pdb.write2file("write3.pdb");

   return 0;
}

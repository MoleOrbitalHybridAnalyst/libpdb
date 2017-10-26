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
   return 0;
}

#include "pdb.h"

using namespace PDB_NS;

int main(int argc, char **argv) {
   PDB pdb(argv[1]);
   pdb.write2file("write.pdb");
   return 0;
}

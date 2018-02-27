#include "utili.cpp"
#include "pdbdef.cpp"
#include "pdb.cpp"

using namespace PDB_NS;

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pypdb)
{
   class_<PDBDef>("PDBDef", init<const std::string&>());
   class_<PDB>("PDB", init<const std::string&>())
      .add_property("write2file", &PDB::write2file);
   //class_<PDBDef>("PDBDef", init<std::string>());
}


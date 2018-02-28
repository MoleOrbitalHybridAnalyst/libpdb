#include "utili.cpp"
#include "pdbdef.cpp"
#include "pdb.cpp"

using namespace PDB_NS;

#include <boost/python.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pypdb)
{
   enum_<PDBField>("pdb_field")
      .value("atomname", PDBField::atomname)
      .value("resname", PDBField::resname)
      .value("segname", PDBField::segname)
      .value("atomtype", PDBField::atomtype)
      .value("chainid", PDBField::chainid)
      .value("resid", PDBField::resid)
      .value("x", PDBField::x)
      .value("y", PDBField::y)
      .value("z", PDBField::z)
      .value("occ", PDBField::occ)
      .value("tempf", PDBField::tempf)
      .value("unknown", PDBField::unknown);

   void (PDBDef::*pushBack0) (PDBField, const std::string&) = &PDBDef::pushBack;
   void (PDBDef::*pushBack1) (PDBField, const float) = &PDBDef::pushBack;
   void (PDBDef::*pushBack2) (PDBField, const int) = &PDBDef::pushBack;
   void (PDBDef::*pushBack3) (PDBField, const char) = &PDBDef::pushBack;
   class_<PDBDef>("pdb_def")
      .def(init<const std::string&>())
      .def("push_back", pushBack0)
      .def("push_back", pushBack1)
      .def("push_back", pushBack2)
      .def("push_back", pushBack3);

   size_t (PDB::*reorderWater1)
      (const PDBDef&, const PDBDef&, const PDBDef&) = &PDB::reorderWater;
   class_<PDB>("pdb_obj")
      .def(init<const std::string&>())
      .def("write2file", &PDB::write2file)
      //.def("getAtomname", &PDB::getAtomname);
      .add_property("x", &PDB::getX, &PDB::setX)
      .def_readonly("natom", &PDB::getNatoms)
      .def("set_atomname", &PDB::setAtomname)
      .def("reorder_water", reorderWater1);
}


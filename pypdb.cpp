#include "utili.cpp"
#include "pdbdef.cpp"
#include "pdb.cpp"

using namespace PDB_NS;

#include <boost/python.hpp>
using namespace boost::python;

template<class T>
struct VecToList
{
   static PyObject* convert(const std::vector<T>& vec)
   {
       list* l = new boost::python::list();
       for(size_t i = 0; i < vec.size(); i++) {
           l->append(vec[i]);
       }

       return l->ptr();
   }
};

BOOST_PYTHON_MODULE(pypdb)
{
   to_python_converter<std::vector<size_t>, VecToList<size_t>>();

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
      .def("push_back", pushBack3)
      .def("print", &PDBDef::print);

   size_t (PDB::*reorderWater0) (bool, bool, bool, 
         const PDBDef&, const PDBDef&, const PDBDef&) = &PDB::reorderWater;
   size_t (PDB::*reorderWater1)
      (const PDBDef&, const PDBDef&, const PDBDef&) = &PDB::reorderWater;
   class_<PDB>("pdb_obj")
      .def(init<const std::string&>())
      .def("write2file", &PDB::write2file)
      //.def("getAtomname", &PDB::getAtomname);
      .add_property("x", &PDB::getX, &PDB::setX)
      .def_readonly("natom", &PDB::getNatoms)
      .def("set_atomname", &PDB::setAtomname)
      .def("reorder_water", reorderWater0)
      .def("reorder_water", reorderWater1)
      .def("select_atoms", &PDB::selectAtoms);
}


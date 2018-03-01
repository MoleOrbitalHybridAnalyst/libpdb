#include "utili.h"
#include "pdbdef.h"
#include "pdb.h"

using namespace PDB_NS;

#include "pypdb.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

/// convert my Vector class to python list
//struct Vector2List
//{
//   static PyObject* convert(const Vector& v)
//   {
//      list* l = new list();
//     for(int i = 0; i < 3; ++i) 
//        l->append(v[i]);
//
//     return l->ptr();
//   }
//};

BOOST_PYTHON_MODULE(pypdb)
{
   //obsolete after I use vector_indexing_suite:
   //to_python_converter<std::vector<size_t>, VecToList<size_t>>();
   //to_python_converter<std::vector<std::string>, VecToList<std::string>>();
   to_python_converter<Vector, Vector2List<Vector>();

/// wrap vectors 
/// http://www.boost.org/doc/libs/1_51_0/libs/python/doc/v2/indexing.html
   class_<std::vector<std::string>>("std_vector_string")
      .def(vector_indexing_suite<std::vector<std::string>>());
   class_<std::vector<float>>("std_vector_float")
      .def(vector_indexing_suite<std::vector<float>>());
   class_<std::vector<size_t>>("std_vector_size_t")
      .def(vector_indexing_suite<std::vector<size_t>>());
   class_<std::vector<int>>("std_vector_int")
      .def(vector_indexing_suite<std::vector<int>>());

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
   void (PDBDef::*pushBack3) (const std::string&) = &PDBDef::pushBack;
   //void (PDBDef::*pushBack3) (PDBField, const char) = &PDBDef::pushBack;
   class_<PDBDef>("pdb_def")
      .def(init<const std::string&>())
      .def("__copy__", &generic__copy__<PDBDef>)
      .def("__deepcopy__", &generic__deepcopy__<PDBDef>)
      .def("show", &PDBDef::print)
      .def("push_back", pushBack0)
      .def("push_back", pushBack1)
      .def("push_back", pushBack2)
      .def("push_back", pushBack3)
      .def("pop_back", &PDBDef::popBack)
   ;

   size_t (PDB::*reorderWater0) (bool, bool, bool, 
         const PDBDef&, const PDBDef&, const PDBDef&) = &PDB::reorderWater;
   size_t (PDB::*reorderWater1)
      (const PDBDef&, const PDBDef&, const PDBDef&) = &PDB::reorderWater;
   void (PDB::*swapFields0)
      (const size_t, const size_t, const PDBField&) = &PDB::swapFields;
   void (PDB::*swapFields1)
      (const size_t, const size_t) = &PDB::swapFields;
   void (PDB::*assembleWater0) (bool, bool,
         const PDBDef&, const PDBDef&) = &PDB::assembleWater;
   void (PDB::*assembleWater1) 
      (const PDBDef&, const PDBDef&) = &PDB::assembleWater;
   class_<PDB>("pdb_obj")
      .def(init<const std::string&>())
      .def("__copy__", &generic__copy__<PDB>)
      .def("__deepcopy__", &generic__deepcopy__<PDB>)
      .def_readonly("natom", &PDB::getNatoms)
      .def("write2file", &PDB::write2file)
      .def("print_atom", static_cast<
            void (PDB::*)(size_t) const>(&PDB::printOneAtom))
      .def("print_atom", static_cast<
            void (PDB::*)(const PDBDef&) const>(&PDB::printAtoms))
      .add_property("atomnames", make_function(
               &PDB::getAtomnames, return_internal_reference<>()))
      .add_property("resnames", make_function(
               &PDB::getResnames, return_internal_reference<>()))
      .add_property("segnames", make_function(
               &PDB::getSegnames, return_internal_reference<>()))
      .add_property("atomtypes", make_function(
               &PDB::getAtomtypes, return_internal_reference<>()))
      .add_property("chainids", make_function(
               &PDB::getChainids, return_internal_reference<>()))
      .add_property("resids", make_function(
               &PDB::getResids, return_internal_reference<>()))
      .add_property("xs", make_function(
               &PDB::getXs, return_internal_reference<>()))
      .add_property("ys", make_function(
               &PDB::getYs, return_internal_reference<>()))
      .add_property("zs", make_function(
               &PDB::getZs, return_internal_reference<>()))
      .add_property("occs", make_function(
               &PDB::getOccs, return_internal_reference<>()))
      .add_property("tempfs", make_function(
               &PDB::getTempfs, return_internal_reference<>()))
      .def("select_atoms", &PDB::selectAtoms)
      .def("reorder_water", reorderWater0)
      .def("reorder_water", reorderWater1)
      .def("guess_chainids", &PDB::guessAllChainids)
      .def("guess_segnames", &PDB::guessAllSegnames)
      .def("guess_atomtypes", &PDB::guessAllAtomtypes)
      .def("swap", swapFields0)
      .def("swap", swapFields1)
      .def("assemble_water", assembleWater0)
      .def("assemble_water", assembleWater1)
      .def("geo_center", static_cast<
            Vector (PDB::*)(const PDBDef&) const >(&PDB::geoCenter))
   ;
}


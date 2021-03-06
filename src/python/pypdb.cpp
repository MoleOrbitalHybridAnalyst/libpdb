#include "utili.h"
#include "pdbdef.h"
#include "pdb.h"

using namespace PDB_NS;

#include "pypdb.h"
#include <boost/python/suite/indexing/vector_indexing_suite.hpp>
using namespace boost::python;

BOOST_PYTHON_MODULE(pypdb_core)
{
   //obsolete after I use vector_indexing_suite:
   //to_python_converter<std::vector<size_t>, VecToList<size_t>>();
   //to_python_converter<std::vector<std::string>, VecToList<std::string>>();
   //to_python_converter<Vector, SubscrToList<Vector>>();

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
   class_<Vector>("Vector")
      .def(vector_indexing_suite<Vector>());

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

   enum_<Direction>("Direction")
      .value("forward", Direction::forward)
      .value("backward", Direction::backward)
      .value("both", Direction::both)
      .value("either", Direction::either);

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
      .def("clear", static_cast<void (PDBDef::*)(PDBField) >(&PDBDef::clear))
      .def("clear", static_cast<
            void (PDBDef::*)(const std::string&) >(&PDBDef::clear))
   ;

   void (PDB::*write2file0) (const std::string&) const = &PDB::write2file;
   void (PDB::*write2file1) 
      (const std::string&, const PDBDef&) const = &PDB::write2file;
   void (PDB::*write2file2) 
      (const PDBDef&, const std::string&) const = &PDB::write2file;
   void (PDB::*write2file3) 
      (const std::string&, const std::vector<size_t>&) const = &PDB::write2file;
   void (PDB::*write2file4) 
      (const std::vector<size_t>&, const std::string&) const = &PDB::write2file;
   void (PDB::*show0) () const = &PDB::show;
   void (PDB::*show1) (const PDBDef&) const = &PDB::show;
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
      .def("write2file", write2file0)
      .def("write2file", write2file1)
      .def("write2file", write2file2)
      .def("write2file", write2file3)
      .def("write2file", write2file4)
      .def("show", show0)
      .def("show", show1)
      .def("print_atom", static_cast<
            void (PDB::*)(size_t) const>(&PDB::printOneAtom))
      .def("print_atom", static_cast<
            void (PDB::*)(const PDBDef&) const>(&PDB::printAtoms))
      .def("print_atom", static_cast<
            void (PDB::*)(const std::vector<size_t>&) const>(&PDB::printAtoms))
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
      .add_property("residues", make_function(
               &PDB::getResidues, return_internal_reference<>()))
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
      .def("reorder_water_fast", &PDB::reorderWaterFast)
      .def("guess_chainids", &PDB::guessAllChainids)
      .def("guess_segnames", &PDB::guessAllSegnames)
      .def("guess_atomtypes", &PDB::guessAllAtomtypes)
      .def("swap", swapFields0)
      .def("swap", swapFields1)
      .def("assemble_water", assembleWater0)
      .def("assemble_water", assembleWater1)
      .def("geo_center", static_cast<
            Vector (PDB::*)(const PDBDef&) const >(&PDB::geoCenter))
      .def("geo_center", static_cast<
            Vector (PDB::*)(const std::vector<size_t>&) const>(&PDB::geoCenter))
      .def("shift2middle", static_cast<
            Vector (PDB::*)(const std::vector<size_t>&)> (&PDB::shiftToMiddle))
      .def("shift2middle", static_cast<
            Vector (PDB::*)(const PDBDef&)> (&PDB::shiftToMiddle))
      .def("get_hb_acceptors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const std::vector<size_t>&, size_t) 
            >(&PDB::getHBAcceptors))
      .def("get_hb_acceptors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const PDBDef&, const PDBDef&) 
            >(&PDB::getHBAcceptors))
      .def("get_hb_donors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const std::vector<size_t>&, size_t) 
            >(&PDB::getHBDonors))
      .def("get_hb_donors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const PDBDef&, const PDBDef&) 
            >(&PDB::getHBDonors))
      .def("get_hb_acceptors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const std::vector<size_t>&, size_t, size_t) 
            >(&PDB::getHBAcceptors))
      .def("get_hb_acceptors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const PDBDef&, const PDBDef&, size_t) 
            >(&PDB::getHBAcceptors))
      .def("get_hb_donors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const std::vector<size_t>&, size_t, size_t) 
            >(&PDB::getHBDonors))
      .def("get_hb_donors", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, const PDBDef&, const PDBDef&, size_t) 
            >(&PDB::getHBDonors))
      .def("get_hb_network", static_cast<
            std::vector<size_t> (PDB::*)(
               int, float, float, const PDBDef&, 
               const PDBDef&, const PDBDef&, const PDBDef&, Direction, bool)
            >(&PDB::getHBNetwork))
      .def("atoms_within", static_cast<
            std::vector<size_t> (PDB::*)(
              const list&, float) const> (&PDB::atomsWithin) )
      .def("write2xyz", static_cast<
            void (PDB::*) (const std::string&, const std::string&) const
            > (&PDB::writeXYZ) )
      .def("pbc_distance", static_cast<
            Vector (PDB::*) (size_t, size_t) const
            > (&PDB::pbcDistance) )
      .def("pbc_distance", static_cast<
            Vector (PDB::*) (const Vector&, const Vector&) const
            > (&PDB::pbcDistance) )
      .def("pbc_distance2", static_cast<
            float (PDB::*) (size_t, size_t) const
            > (&PDB::pbcDistance2) )
      .def("pbc_distance2", static_cast<
            float (PDB::*) (const Vector&, const Vector&) const
            > (&PDB::pbcDistance2) )
      .def("write2ndx", static_cast<
            void (PDB::*) (const std::string&, const std::string&) const
            > (&PDB::writeIndexFile) )
      .def("moveto", &PDB::moveTo)
      .def("get_coordinates", &PDB::getCoordinates)
   ;
}


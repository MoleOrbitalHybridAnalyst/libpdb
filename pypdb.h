#include <vector>
#include <string>

#include <boost/python.hpp>

namespace bp = boost::python;

/// convert vector into python list (obsolete after I use vector_indexing_suite)
//template<class T>
//struct VecToList
//{
//   static PyObject* convert(const std::vector<T>& vec)
//   {
//      bp::list* l = new bp::list();
//       for(size_t i = 0; i < vec.size(); i++) {
//           l->append(vec[i]);
//       }
//
//       return l->ptr();
//   }
//};

/// convert c++ subscribable into python list
template<class T>
struct SubscrToList
{
   static PyObject* convert(const T& subscr)
   {
       Py_intptr_t shape[1] = { subscr.size() };
       bn::ndarray result = bn::zeros(1, shape, bn::dtype::get_builtin<double>());
       std::copy(v.begin(), v.end(), reinterpret_cast<double*>(result.get_data()));
       return result.ptr();

       //bp::list* l = new bp::list();
       //for(size_t i = 0; i < vec.size(); i++) {
       //    l->append(vec[i]);
       //}

       //return l->ptr();
   }
};

/// provide python __copy__ and __deepcopy__ for bp object
/// https://mail.python.org/pipermail/cplusplus-sig/2009-May/014505.html
template<class T>
inline PyObject * managingPyObject(T *p)
{
    return typename bp::manage_new_object::apply<T *>::type()(p);
}

template<class Copyable>
bp::object
generic__copy__(bp::object copyable)
{
    Copyable *newCopyable(new Copyable(bp::extract<const Copyable
&>(copyable)));
    bp::object
result(bp::detail::new_reference(managingPyObject(newCopyable)));

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        copyable.attr("__dict__"));

    return result;
}

template<class Copyable>
bp::object
generic__deepcopy__(bp::object copyable, bp::dict memo)
{
    bp::object copyMod = bp::import("copy");
    bp::object deepcopy = copyMod.attr("deepcopy");

    Copyable *newCopyable(new Copyable(bp::extract<const Copyable
&>(copyable)));
    bp::object
result(bp::detail::new_reference(managingPyObject(newCopyable)));

    long copyableId = (long)(copyable.ptr());
    memo[copyableId] = result;

    bp::extract<bp::dict>(result.attr("__dict__"))().update(
        deepcopy(bp::extract<bp::dict>(copyable.attr("__dict__"))(),
memo));

    return result;
}

/// wrap [] operator for vector (obsolete afer I found vector_indexing_suite)
/// https://mail.python.org/pipermail/cplusplus-sig/2003-June/004024.html
//template <class T>
//void vector_setitem(std::vector<T>& v, size_t index, T value)
//{
//    if (index < v.size()) {
//        v[index] = value;
//    }
//    else {
//        PyErr_SetString(PyExc_IndexError, "index out of range");
//        bp::throw_error_already_set();
//    }
//}
//
//template <class T>
//T vector_getitem(std::vector<T> &v, size_t index)
//{
//    if (index < v.size()) {
//        return v[index];
//    }
//    else {
//        PyErr_SetString(PyExc_IndexError, "index out of range");
//        bp::throw_error_already_set();
//    }
//}


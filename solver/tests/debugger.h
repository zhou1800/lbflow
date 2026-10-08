/*! \file debugger.h
 * \author Yu Zhou
 */
#ifndef DEBUGGER_H
#define DEBUGGER_H

#include "defs.h"

#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace DRL {

namespace TEST {

using std::ostream;
using std::map;
using std::string;
using std::gslice;
using std::vector;
using std::endl;

/*!
 * \brief The Debugger class is for generating the debug information.
 *
 * Useage:
 *  - Construct with target ostream for output.
 *  - Add() links to variables needs to be debug. Accept almost all kinds of
 *    variables.
 *  - Show() human readable debug information to ostream.
 */
class Debugger
{
private:
    /// \name Alias.
    ///@{
    using va_string    = valarray<string>;      //!< Array of string.
    using va_gslice    = valarray<gslice>;      //!< Array of gslice.
    using ve_idx_t    = vector<size_t>;        //!< Vector of size.
    using ve_va_idx_t = vector<va_idx_t>;     //!< Vector of array of size.
    using ve_double    = vector<double>;        //!< Vector of double.
    using ve_string    = vector<string>;        //!< Vector of string.
    using ve_gslice    = vector<gslice>;        //!< Vector of gslice.

    ///@}

    static const size_t MAX_LEN = 10000;  //!< Max length for showing a valarray.
    ostream& out;                         //!< The ostream for output.

    /// \name Linked variables to be debugged
    ///@{
    string m_title;                                    //!< Title in output.

    map<string, const size_t*>         m_size_t;       //!< size_t.
    map<string, const va_idx_t*>     m_va_idx_t;    //!< va_idx_t.
    map<string, const vava_idx_t*>  m_va_va_idx_t; //!< vava_idx_t.
    map<string, const ve_idx_t*>      m_ve_idx_t;    //!< vector of size_t.
    map<string, const ve_va_idx_t*>   m_ve_va_idx_t; //!< vector of va_idx_t.

    map<string, const va_int_t*>      m_va_int;       //!< arr_int.
    map<string, const vava_int_t*>  m_va_va_int;    //!< arr_arr_int.

    map<string, const double*>         m_double;       //!< double.
    map<string, const va_val_t*>      m_va_double;    //!< arr of double.
    map<string, const vava_val_t*>  m_va_va_val;    //!< arr_arr_double.
    map<string, const ve_double*>      m_ve_double;    //!< vector of double.

    map<string, const string*>         m_string;       //!< string.
    map<string, const va_string*>      m_va_string;    //!< arr of string.
    map<string, const ve_string*>      m_ve_string;    //!< vector of string.

    map<string, const gslice*>         m_gslice;       //!< gslice type.
    map<string, const va_gslice*>      m_va_gslice;    //!< arr of gslice.
    map<string, const ve_gslice*>      m_ve_gslice;    //!< vector of gslice.

    ///@}

public:
    ///@{
    /*!
     * \brief Constructor.
     *
     * Link an ostream for output.
     *
     * \param os - ostream for output.
     */
    Debugger(ostream& os);

    /// \name Add variables to Debugger
    ///@{
    /*!
     * \brief Add title.
     * \param title - title to be add to output.
     */
    void Add(string title) { m_title = title; }

    /*!
     * \brief Add a size_t variable.
     * \param name - name of the variable.
     * \param v - value.
     */
    void Add(string name, const size_t& v) { m_size_t[name] = &v; }

    /*!
     * \brief Add a valarray of size_t variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const va_idx_t& v) { m_va_idx_t[name] = &v; }

    /*!
     * \brief Add a valarray of valarray of size_t variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const vava_idx_t& v) { m_va_va_idx_t[name] = &v; }

    /*!
     * \brief Add a vector of size_t variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const ve_idx_t& v) { m_ve_idx_t[name] = &v; }

    /*!
     * \brief Add a vector of array of size_t variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const ve_va_idx_t& v) { m_ve_va_idx_t[name] = &v; }

    /*!
     * \brief Add a valarray of int variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const va_int_t& v) { m_va_int[name] = &v; }

    /*!
     * \brief Add a valarray of valarray of int variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const vava_int_t& v) { m_va_va_int[name] = &v; }

    /*!
     * \brief Add a double percision float variable.
     * \param name - name of the variable.
     * \param v - value.
     */
    void Add(string name, const double& v) { m_double[name] = &v; }

    /*!
     * \brief Add a valarray of double precision float variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const va_val_t& v) { m_va_double[name] = &v; }

    /*!
     * \brief Add a valarray of valarray of double variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const vava_val_t& v) { m_va_va_val[name] = &v; }

    /*!
     * \brief Add a vector of double precision float variable.
     * \param name - name of the array.
     * \param v - value of the array.
     */
    void Add(string name, const ve_double& v) { m_ve_double[name] = &v; }

    /*!
     * \brief Add a string.
     * \param name - name of the string.
     * \param v - string value.
     */
    void Add(string name, const string& v) { m_string[name] = &v; }

    /*!
     * \brief Add a valarray of string.
     * \param name - name of the array.
     * \param v - value of the string array.
     */
    void Add(string name, const va_string& v) { m_va_string[name] = &v; }

    /*!
     * \brief Add a vector of string.
     * \param name - name of the array.
     * \param v - value of the string array.
     */
    void Add(string name, const ve_string& v) { m_ve_string[name] = &v; }

    /*!
     * \brief Add a gslice.
     * \param name - name of the gslice.
     * \param v - value of the gslice.
     */
    void Add(string name, const gslice& v) { m_gslice[name] = &v; }

    /*!
     * \brief Add a valarray of gslices.
     * \param name - name of the array.
     * \param v - value of the gslice array.
     */
    void Add(string name, const va_gslice& v) { m_va_gslice[name] = &v; }

    /*!
     * \brief Add a vector of gslices.
     * \param name - name of the vector.
     * \param v - value of the gslice vector.
     */
    void Add(string name, const ve_gslice& v) { m_ve_gslice[name] = &v; }
    ///@}

    /// \name Output operations
    ///@{
    /*!
     * \brief Show debug information.
     *
     * Generate debug information to ostream.
     */
    void Show() const;

    /*!
     * \brief Clear all variables.
     */
    void Clear();
    ///@}

private:
    /// \name Helper functions for Show()
    ///@{
    /*!
     * \brief Show a value.
     *
     * Generate human readable information to output ostream. The value must be
     *  supported by ostream.operator<<().
     *
     * \param a value to be shown.
     */
    template <typename T>
    void showValue(const T& a) const;

    /*!
     * \brief Show an array of value.
     *
     * Generate human readable information to output ostream. The values in the
     *  array must be supported by ostream.operator<<().
     *
     * \param a - the array of value to be shown.
     */
    template <typename T>
    void showArray(const T& a) const;

    /*!
     * \brief Show an array of arrays of value.
     *
     * Generate human readable information to output ostream. The value must be
     *  supported by ostream.operator<<().
     *
     * \param a value to be shown.
     */
    template <typename T>
    void showArrArr(const T& a) const;

    /*!
     * \brief Show a gslice value.
     *
     * Generate human readable information of gslice to output ostream.
     *
     * \param gs - the array of value to be shown.
     */
    void showGslice(const gslice& gs) const;

    template <typename T>
    /*!
     * \brief Show an array of gslice.
     *
     * Generate human readable information of an array of gslice to output
     *  ostream.
     *
     * \param a - array of gslice.
     */
    void showGsArr(const T& a) const;
    ///@}
};

template <typename T> inline
void Debugger::showValue(const T &a) const
{
    out << a.first << " = " << *a.second << ", ";
}

template <typename T> inline
void Debugger::showArray(const T &a) const
{
    out << a.first << "={ ";
    if ((*a.second).size() > MAX_LEN)
        out << "Warnning! Size exceeded (size = "
            << (*a.second).size() << ")! }";
    else {
        for (const auto& s : *a.second)
            out << s << " ";
    }
    out << "}";
}

template <typename T> inline
void Debugger::showArrArr(const T &a) const
{
    out << a.first << "= ";
    if ((*a.second).size() > MAX_LEN)
        out << "{Warnning! Size exceeded (size = "
            << (*a.second).size() << ")! }";
    else {
        out <<endl;
        for (size_t i = 0; i < (*a.second).size(); ++i) {
            out << a.first << "-" << i << "={";
            if ((*a.second)[i].size() > MAX_LEN)
                out << "Warnning! Size exceeded (size = "
                    << (*a.second)[i].size() << ")! }";
            else {
                for (const auto& s : (*a.second)[i])
                    out << s << " ";
            }
            out << "}" << endl;
        }
    }
}

inline
void Debugger::showGslice(const gslice &gs) const
{
    out << "{";
    out << gs.start() << ", [ ";
    for (const auto& v : gs.size())
        out << v << " ";
    out << "], [ ";
    for (const auto& v : gs.stride())
        out << v << " ";
    out << "]}";
}

template <typename T> inline
void Debugger::showGsArr(const T &a) const
{
    out << a.first << "=" << endl;
    if ((*a.second).size() > MAX_LEN)
        out << "Warnning! Size exceeded (size = "
            << (*a.second).size() << ")!" << endl;
    else {
        for (const auto& gs : *a.second) {
            showGslice(gs);
            out << endl;
        }
    }

}

} // namespace TEST

} // namespace DRL

#endif // DEBUGGER_H

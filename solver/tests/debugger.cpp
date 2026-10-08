/*! \file debugger.cpp
 * \author Yu Zhou
 */
#include "debugger.h"
#include "thread.h"

namespace DRL {

namespace TEST {

using std::pair;

Debugger::Debugger(std::ostream &os)
    : out(os)
{
}

void Debugger::Show() const
{
    // show title and mpi info:
    out << "[processor " << mpi().getId() << "of" << mpi().getCnt();
    if (m_title != "")
        out << "-" << m_title;
    out << "]" << endl;

    // show m_size_t
    if (m_size_t.size() != 0) {
        for(const auto& a : m_size_t)
            showValue(a);
        out << endl;
    }

    // show m_va_idx_t
    if (m_va_idx_t.size() != 0)
        for(const auto& a : m_va_idx_t) {
            showArray(a);
            out << endl;
        }

    // show m_va_va_idx_t
    if (m_va_va_idx_t.size() != 0)
        for(const auto& a : m_va_va_idx_t) {
            showArrArr(a);
            out << endl;
        }

    // show m_ve_idx_t
    if (m_ve_idx_t.size() != 0)
        for(const auto& a : m_ve_idx_t) {
            showArray(a);
            out << endl;
        }

    // show m_ve_va_idx_t
    if (m_ve_va_idx_t.size() != 0)
        for(const auto& a : m_ve_va_idx_t) {
            showArrArr(a);
            out << endl;
        }

    // show m_va_int
    if (m_va_int.size() != 0)
        for(const auto& a : m_va_int) {
            showArray(a);
            out << endl;
        }

    // show m_va_va_int
    if (m_va_va_int.size() != 0)
        for(const auto& a : m_va_va_int) {
            showArrArr(a);
            out << endl;
        }

    // show m_double
    if (m_double.size() != 0) {
        for(const auto& a : m_double)
            showValue(a);
        out << endl;
    }

    // show m_va_double
    if (m_va_double.size() != 0)
        for(const auto& a : m_va_double) {
            showArray(a);
            out << endl;
        }

    // show m_va_va_double
    if (m_va_va_val.size() != 0)
        for(const auto& a : m_va_va_val) {
            showArrArr(a);
            out << endl;
        }

    // show m_ve_double
    if (m_ve_double.size() != 0)
        for(const auto& a : m_ve_double) {
            showArray(a);
            out << endl;
        }

    // show m_string
    if (m_string.size() != 0) {
        for(const auto& a : m_string)
            showValue(a);
        out << endl;
    }

    // show m_va_string
    if (m_va_string.size() != 0)
        for(const auto& a : m_va_string) {
            showArray(a);
            out << endl;
        }

    // show m_ve_string
    if (m_ve_string.size() != 0)
        for(const auto& a : m_ve_string) {
            showArray(a);
            out << endl;
        }

    // show m_gslice
    if ( m_gslice.size() != 0) {
        for(const auto& a : m_gslice) {
            out << a.first << "=";
            showGslice(*a.second);
            out << endl;
        }
    }

    // show m_va_gslice
    if ( m_va_gslice.size() != 0) {
        for(const auto& a : m_va_gslice) {
            showGsArr(a);
        }
    }

    // show m_ve_gslice
    if ( m_ve_gslice.size() != 0) {
        for(const auto& a : m_ve_gslice) {
            showGsArr(a);
        }
    }
}

void Debugger::Clear()
{
    m_size_t.clear();
    m_va_idx_t.clear();
    m_va_va_idx_t.clear();
    m_ve_idx_t.clear();
    m_va_int.clear();
    m_va_va_int.clear();
    m_double.clear();
    m_va_double.clear();
    m_ve_double.clear();
    m_string.clear();
    m_va_string.clear();
    m_ve_string.clear();
    m_gslice.clear();
    m_va_gslice.clear();
}

} // namespace TEST

} // namespace DRL

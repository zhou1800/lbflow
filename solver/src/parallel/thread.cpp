/*! \file thread.cpp
 * \author Yu Zhou
 */
#include "thread.h"
#include "sstream"

namespace DRL {

using std::runtime_error;
using std::stringstream;

Threads::Threads()
{
    MPI_Init(0, 0);
    MPI_Comm_rank(MPI_COMM_WORLD, &id);
    MPI_Comm_size(MPI_COMM_WORLD, &sz);
}

Threads::~Threads()
{
    MPI_Finalize();
}

void *Threads::Pack(MPI_Datatype *dataType, char **pStrs, int *strLen,
                    unsigned strCnt)
{
    MPI_Datatype*   arrType = new MPI_Datatype[strCnt];
    MPI_Aint*       arrDisplace = new MPI_Aint[strCnt];

    for (size_t i = 0; i < strCnt; ++i) {
        arrType[i] = MPI_CHAR;
        if (MPI_Get_address(pStrs[i], &arrDisplace[i]) != MPI_SUCCESS){
            stringstream errMsg;
            errMsg << "MPI_Get_address() failed to get address from string \'"
                   << pStrs[i] << "\' in Threads::Pack(string).";
            throw runtime_error(errMsg.str());
        }
    }
    for (auto i = strCnt; i > 0; --i)
        arrDisplace[i-1] -= arrDisplace[0];

    if (MPI_Type_create_struct(strCnt, strLen, arrDisplace, arrType,
                               dataType) != MPI_SUCCESS)
        throw runtime_error("MPI_Type_create_struct() failed in "
                            "Threads::pack(string).");
    if (MPI_Type_commit(dataType) != MPI_SUCCESS)
        throw runtime_error("MPI_Type_commit() failed in "
                            "Threads::pack(string).");

    delete [] arrType;
    delete [] arrDisplace;
    return pStrs[0];
}

void* Threads::Pack(MPI_Datatype* dataType, GsliceOpen *gs) const
{
    const int       count = 3;
    int size = static_cast<int>(gs->size.size());
    int             arrBlockLength[] = {1, size, size};
    MPI_Datatype    arrType[] = {MPI_UNSIGNED_LONG, MPI_UNSIGNED_LONG,
                                 MPI_UNSIGNED_LONG};
    MPI_Aint        arrDisplace[count];

    if (MPI_Get_address(&(gs->start), &arrDisplace[0]) != MPI_SUCCESS)
        throw runtime_error("MPI_Get_address() failed when getting address for "
                            "GsliceOpen gs::start in "
                            "Threads::pack(GsliceOpen).");
    if (MPI_Get_address(&(gs->size[0]), &arrDisplace[1]) != MPI_SUCCESS)
        throw runtime_error("MPI_Get_address() failed when getting address for "
                            "GsliceOpen gs::size in "
                            "Threads::pack(GsliceOpen).");
    if (MPI_Get_address(&(gs->stride[0]), &arrDisplace[2]) != MPI_SUCCESS)
        throw runtime_error("MPI_Get_address() failed when getting address for "
                            "GsliceOpen gs::stride in "
                            "Threads::pack(GsliceOpen).");

    for (auto i = count; i > 0; --i)
        arrDisplace[i-1] -= arrDisplace[0];

    if (MPI_Type_create_struct(count, arrBlockLength, arrDisplace, arrType,
                               dataType) != MPI_SUCCESS)
        throw runtime_error("MPI_Type_create_struct() failed in "
                            "Threads::pack(GsliceOpen).");
    if (MPI_Type_commit(dataType) != MPI_SUCCESS)
        throw runtime_error("MPI_Type_commit() failed in "
                            "Threads::pack(GsliceOpen).");

    return &(gs->start);
}

} // namespace DRL

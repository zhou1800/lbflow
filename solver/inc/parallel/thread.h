/*! \file thread.h
 *
 * \author Yu Zhou (zhouyu@live.com)
 */
#ifndef THREAD_H
#define THREAD_H

#include "defs.h"
#include <stdexcept>
#include <iostream>
#include "mpi.h"

namespace DRL {

class GsliceOpen;

/*!
 * \brief The Thread class wrap up mpi functions.
 */
class Threads
{
private:
    int id;                     //!< Local rank.
    int sz;                     //!< Communication size.
    mutable MPI_Status status;  //!< Mpi status after receive.

public:
    /// \name Constructors and destructor
    ///@{
    /*!
     * \brief Construct threads.
     *
     * Initialize mpi with MPI_Init(), MPI_Comm_rank(), and MPI_Comm_size().
     */
    Threads();

    /*!
     * \brief Destructor.
     *
     * Release mpi resource with MPI_Finalize().
     */
    ~Threads();
    ///@}

    /// \name MPI transfer functions
    ///@{
    /*!
     * \brief Implements MPI_Send().
     *
     * Prefer Ssend() for safety.
     *
     * \sa Ssend(), Recv()
     * \exception std::runtime_error on MPI operator failed.
     * \param buf - buffer to be sent.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param tag - tag attached to the message.
     * \param dest - destination rank(id).
     */
    void Send(void *buf, size_t size, MPI_Datatype type, int tag, int dest)
    const;

    /*!
     * \brief Implements MPI_Ssend().
     * \sa Send(), Recv()
     * \exception std::runtime_error on MPI operator failed.
     * \param buf - buffer to be sent.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param tag - tag attached to the message.
     * \param dest - destination rank(id).
     */
    void Ssend(void *buf, size_t size, MPI_Datatype type, int tag, int dest)
    const;

    /*!
     * \brief Implements MPI_Recv().
     *
     * MPI status is updated after receive.
     *
     * \sa Send(), Ssend()
     * \exception std::runtime_error on MPI operator failed.
     * \param buf - buffer to store the received data.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param tag - tag attached to the message.
     * \param src - source rank(id).
     */
    void Recv(void* buf, size_t size, MPI_Datatype type, int tag = MPI_ANY_TAG,
              int src = MPI_ANY_SOURCE) const;

    /*!
     * \brief Implements MPI_Bcast().
     * \exception std::runtime_error on MPI operator failed.
     * \param buf - buffer to send and receive data.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param src - source rank(id).
     */
    void Bcast(void* buf, size_t size, MPI_Datatype type, int src) const;

    /*!
     * \brief Implements MPI_Scatter().
     *
     * Share the same size and type for both sendbuf and recvBuf.
     *
     * \sa Gather()
     * \exception std::runtime_error on MPI operator failed.
     * \param sendBuf - buffer to be sent.
     * \param recvBuf - buffer to store the received data.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param src - source rank(id).
     */
    void Scatter(void *sendBuf, void* recvBuf, size_t size, MPI_Datatype type,
                 int src) const;

    /*!
     * \brief Implements MPI_Gather().
     *
     * Share the same size and type for both sendbuf and recvBuf.
     *
     * \sa Scatter()
     * \exception std::runtime_error on MPI operator failed.
     * \param sendBuf - buffer to be sent.
     * \param recvBuf - buffer to store the received data.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param dest - destination rank(id).
     */
    void Gather(void* sendBuf, void* recvBuf, size_t size, MPI_Datatype type,
                int dest) const;

    /*!
     * \brief Implements MPI_Sendrecv().
     *
     * Send sendBuf to dest while receive from src.
     *
     * \exception std::runtime_error on MPI operator failed.
     * \param sendBuf - buffer to be sent.
     * \param recvBuf - buffer to store the received data.
     * \param size - size of buffer.
     * \param type - buffer type.
     * \param tag - tag attached to the message.
     * \param src - source rank(id).
     * \param dest - destination rank(id).
     */
    void SendRecv(void *sendBuf, void* recvBuf, size_t size, MPI_Datatype type,
                  int tag, int src, int dest) const;

    /*!
     * \brief Implements MPI_Barrier().
     *
     * Synchronizes all processes.
     */
    void Barrier() const;

    /*!
     * \brief Implements MPI_Reduce().
     *
     * The destination process is always the head processor (rank 0).
     *
     * \exception std::runtime_error on MPI operator failed.
     * \param inBuf - input buffer.
     * \param outBuf - output buffer.
     * \param size - array size.
     * \param type - data type.
     * \param op - MPI reduce operator.
     */
    void Reduce(void *inBuf, void* outBuf, size_t size, MPI_Datatype type,
                MPI_Op op) const;

    /*!
     * \brief Implements MPI_Allreduce().
     * \exception std::runtime_error on MPI operator failed.
     * \param inBuf - input buffer.
     * \param outBuf - output buffer.
     * \param size - array size.
     * \param type - data type.
     * \param op - MPI reduce operator.
     */
    void AllReduce(void *inBuf, void* outBuf, size_t size, MPI_Datatype type,
                MPI_Op op) const;
    ///@}

    /// \name Get MPI status
    ///@{
    /*!
     * \brief Get processor rank(Id).
     * \return current processor rank(Id).
     */
    int getId() const { return id; }

    /*!
     * \brief Get communication size.
     * \return communication size.
     */
    int getCnt() const { return sz; }

    /*!
     * \brief Get sender (mpi source) id.
     * \return MPI_SOURCE id.
     */
    int getSrc() const { return status.MPI_SOURCE; }

    /*!
     * \brief Get sent tag of pervious message.
     * \return MPI_TAG.
     */
    int getTag() const { return status.MPI_TAG; }

    /*!
     * \brief Get received element count.
     * \param type - received buffer type.
     * \return number of received elements.
     */
    int getCount(MPI_Datatype type) const;

    /*!
     * \brief getTime implements MPI_Wtime().
     * \return current time.
     */
    double getTime() const { return MPI_Wtime(); }
    ///@}

    /// \name MPI pack functions
    ///@{
    /*!
     * \brief Pack an array of strings into MPI_Datatype.
     *
     * Calculate array type and array displacement, then create dataType use
     *  MPI_Type_create_struct, and then commit the dataType with
     *  MPI_Type_commit.
     *
     * \sa use TypeFree() to free the dataType.
     * \exception std::runtime_error on MPI pack operation failure.
     * \param dataType - MPI_Datatype.
     * \param pStrs - pointer to array of strings.
     * \param strLen - array of string length.
     * \param strCnt - number of strings.
     * \return start pointer to the packed data.
     */
    void* Pack(MPI_Datatype* dataType, char** pStrs, int *strLen,
               unsigned strCnt);

    /*!
     * \brief Pack a open gslice (GsliceOpen) into MPI_Datatype.
     *
     * Create a MPI_Datatype for GsliceOpen, then commit the dataType with
     *  MPI_Type_commit.
     *
     * \sa use TypeFree() to free the dataType.
     * \sa the GsliceOpen class for defination of open gslice.
     * \exception std::runtime_error on MPI pack operation failure.
     * \param dataType - MPI_Datatype.
     * \param gs - GsliceOpen subject to be packed.
     * \return start pointer to the packed data.
     */
    void* Pack(MPI_Datatype* dataType, GsliceOpen* gs) const;

    /*!
     * \brief Implements MPI_Type_free().
     *
     * Free user created MPI_Datatype.
     *
     * \sa Pack()
     * \exception std::runtime_error on MPI operator failed.
     * \param dataType - user created MPI_Datatype.
     */
    void TypeFree(MPI_Datatype* dataType);
    ///@}
};

/*!
 * \brief The GsliceOpen class is a replica of std::gslice with public start,
 *  size, and stride.
 *
 * Provides methods to convert between std::gslice and GsliceOpen.
 *
 * This class is intend to be used with Thread::Pack().
 */
class GsliceOpen
{
public:
    /// \name replication of gslice variables.
    ///@{
    size_t start;         //!< Compatible with gslice::start().
    va_idx_t size;        //!< Compatible with gslice::size().
    va_idx_t stride;      //!< Compatible with gslice::stride().

    ///@}

public:
    /// \name Constructors and destructor
    ///@{
    /*!
     * \brief Construct an empty GsliceOpen.
     */
    GsliceOpen() { }

    /*!
     * \brief Construct with a copy from std::gslice.
     * \param gs - std::gslice to be copied.
     */
    GsliceOpen(const gslice& gs) { setGs(gs); }

    /*!
     * \brief Construct by start, size and stride.
     * \param start - gslice start.
     * \param size - gslice size.
     * \param stride - gslice stride.
     */
    GsliceOpen(size_t start, va_idx_t size, va_idx_t stride);
    ///@}

    /// \name Get and set functions
    ///@{
    /*!
     * \brief Set GsliceOpen by std::gslice.
     * \param gs - std::gslice to be copied.
     */
    void setGs(const gslice& gs);

    /*!
     * \brief Set GsliceOpen by start, size, and stride.
     * \param start - gslice start.
     * \param size - gslice size.
     * \param stride - gslice stride.
     */
    void setGs(size_t start, va_idx_t size, va_idx_t stride);

    /*!
     * \brief Get std::gslice from GsliceOpen.
     * \return std::gslice.
     */
    gslice getGs();
    ///@}
};

inline
void Threads::Send(void *buf, size_t size, MPI_Datatype type, int tag, int dest)
const
{
    if (tag == MPI_ANY_TAG) tag = 0;
    if (MPI_Send(buf, static_cast<int>(size), type,
                 dest, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Send() failed in Threads::Send().");
}

inline
void Threads::Ssend(void *buf, size_t size, MPI_Datatype type, int tag,
                    int dest) const
{
    if (tag == MPI_ANY_TAG) tag = 0;
    if (MPI_Ssend(buf, static_cast<int>(size), type,
                  dest, tag, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Ssend() failed in Threads::SSend().");
}

inline
void Threads::Recv(void *buf, size_t size, MPI_Datatype type, int tag, int src)
const
{
    if (MPI_Recv(buf, static_cast<int>(size), type,
                 src, tag, MPI_COMM_WORLD, &status) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Recv() failed in Threads::Recv().");
}

inline
void Threads::Bcast(void *buf, size_t size, MPI_Datatype type, int src) const
{
    if (MPI_Bcast(buf, static_cast<int>(size), type,
                  src, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Bcast() failed in Threads::Bcast().");
}

inline
void Threads::Scatter(void *sendBuf, void *recvBuf, size_t size,
                      MPI_Datatype type, int src) const
{
    if (MPI_Scatter(sendBuf, static_cast<int>(size), type,
                    recvBuf, static_cast<int>(size), type,
                    src, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Scatter() failed in Threads::Scatter().");
}

inline
void Threads::Gather(void *sendBuf, void *recvBuf, size_t size,
                     MPI_Datatype type, int dest) const
{
    if (MPI_Gather(sendBuf, static_cast<int>(size), type,
                   recvBuf, static_cast<int>(size), type,
                   dest, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Gather() failed in Threads::Gather().");
}

inline
void Threads::SendRecv(void *sendBuf, void *recvBuf, size_t size,
                       MPI_Datatype type, int tag, int src, int dest) const
{
    if (MPI_Sendrecv(sendBuf, static_cast<int>(size), type, dest, tag,
                     recvBuf, static_cast<int>(size), type, src, tag,
                     MPI_COMM_WORLD, &status) != MPI_SUCCESS)
        throw std::runtime_error("MPI_Sendrecv() failed in "
                                 "Threads::SendRecv().");
}

inline
void Threads::Barrier() const
{
    MPI_Barrier(MPI_COMM_WORLD);
}

inline
void Threads::Reduce(void *inBuf, void *outBuf, size_t size, MPI_Datatype type,
                     MPI_Op op) const
{
    if (MPI_Reduce(inBuf, outBuf, static_cast<int>(size), type,
                   op, 0, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error ("MPI_Reduce() failed in Threads::Reduce().");
}

inline
void Threads::AllReduce(void *inBuf, void *outBuf, size_t size,
                        MPI_Datatype type, MPI_Op op) const
{
    if (MPI_Allreduce(inBuf, outBuf, static_cast<int>(size), type,
                      op, MPI_COMM_WORLD) != MPI_SUCCESS)
        throw std::runtime_error ("MPI_Allreduce() failed in "
                                  "Threads::AllReduce().");
}

inline
int Threads::getCount(MPI_Datatype type) const
{
    int count;
    if (MPI_Get_count(&status, type, &count) != MPI_SUCCESS)
        throw std::runtime_error ("MPI_Get_count() failed in "
                                  "MpiThread::getCount().");
    return count;
}

inline
void Threads::TypeFree(MPI_Datatype *dataType)
{
    if (MPI_Type_free(dataType) != MPI_SUCCESS)
        throw std::runtime_error ("MPI_Type_free() failed in "
                                  "MpiThread::TypeFree().");
}

inline
GsliceOpen::GsliceOpen(size_t start, va_idx_t size, va_idx_t stride)
{
    setGs(start, size, stride);
}

inline
void GsliceOpen::setGs(const gslice &gs)
{
    start = gs.start();
    size = gs.size();
    stride = gs.stride();
}

inline
void GsliceOpen::setGs(size_t start, va_idx_t size, va_idx_t stride)
{
    this->start = start;
    this->size = size;
    this->stride = stride;
}

inline
gslice GsliceOpen::getGs()
{
    return gslice(start, size, stride);
}

/*!
 * \brief Create and get thread instance.
 *
 * Create the only mpi instance if this file is included.
 * End the programm if failed to create the mpi instance.
 *
 * \return the mpi instance of Threads.
 */
inline
Threads& mpi()
{
    try {
        static Threads instance;
        return instance;
    }
    catch (std::exception& e){
        std::cout << "Failed to start mpi: " << e.what() << std::endl;
        exit(EXIT_FAILURE);
    }
}

} // namespace DRL

#endif // THREADS_H

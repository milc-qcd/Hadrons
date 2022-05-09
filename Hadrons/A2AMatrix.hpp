/*
 * A2AMatrix.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Fionn O hOgain <fionn.o.hogain@ed.ac.uk>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
 * Author: Michael Lynch <michaellynch628@gmail.com>
 *
 * Hadrons is free software: you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation, either version 2 of the License, or
 * (at your option) any later version.
 *
 * Hadrons is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with Hadrons.  If not, see <http://www.gnu.org/licenses/>.
 *
 * See the full license in the file "LICENSE" in the top level distribution 
 * directory.
 */

/*  END LEGAL */
#ifndef A2A_Matrix_hpp_
#define A2A_Matrix_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/TimerArray.hpp>
#include <Grid/Eigen/unsupported/CXX11/Tensor>
#ifdef USE_MKL
#include "mkl.h"
#include "mkl_cblas.h"
#endif

#ifndef HADRONS_A2AM_NAME 
#define HADRONS_A2AM_NAME "a2aMatrix"
#endif

#ifndef HADRONS_A2AM_IO_TYPE
#define HADRONS_A2AM_IO_TYPE ComplexF
#endif

#define HADRONS_A2AM_PARALLEL_IO

BEGIN_HADRONS_NAMESPACE

// general A2A matrix set based on Eigen tensors and Grid-allocated memory
// Dimensions:
//   0 - ext - external field (momentum, EM field, ...)
//   1 - str - spin-color structure
//   2 - t   - timeslice
//   3 - i   - left  A2A mode index
//   4 - j   - right A2A mode index
template <typename T>
using A2AMatrixSet = Eigen::TensorMap<Eigen::Tensor<T, 5, Eigen::RowMajor>>;

template <typename T>
using A2AMatrix = Eigen::Matrix<T, -1, -1, Eigen::RowMajor>;

template <typename T>
using A2AMatrixTr = Eigen::Matrix<T, -1, -1, Eigen::ColMajor>;

/******************************************************************************
 *                      Abstract class for A2A kernels                        *
 ******************************************************************************/
template <typename T, typename Field>
class A2AKernel
{
public:
    A2AKernel(void) = default;
    virtual ~A2AKernel(void) = default;
    virtual void operator()(A2AMatrixSet<T> &m, const Field *left, const Field *right,
                          const unsigned int orthogDim, double *time = nullptr, double *tg = nullptr) = 0;
    virtual double flops(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
    virtual double bytes(const unsigned int blockSizei, const unsigned int blockSizej) = 0;
};

/******************************************************************************
 *                  Class to handle A2A matrix block HDF5 I/O                 *
 ******************************************************************************/
template <typename T>
class A2AMatrixIo
{
public:
    // constructors
    A2AMatrixIo(void) = default;
    A2AMatrixIo(std::string filename, std::string dataname, 
                const unsigned int nt, const unsigned int ni = 0,
                const unsigned int nj = 0);
    // destructor
    ~A2AMatrixIo(void) = default;
    // access
    unsigned int getNi(void) const;
    unsigned int getNj(void) const;
    unsigned int getNt(void) const;
    size_t       getSize(void) const;
    // file allocation
    template <typename MetadataType>
    void initFile(const MetadataType &d, const unsigned int chunkSize);
    // block I/O
    void saveBlock(const T *data, const unsigned int i, const unsigned int j,
                   const unsigned int blockSizei, const unsigned int blockSizej);
    void saveBlock(const A2AMatrixSet<T> &m, const unsigned int ext, const unsigned int str,
                   const unsigned int i, const unsigned int j);
    template <template <class> class Vec, typename VecT>
    void load(Vec<VecT> &v, double *tRead = nullptr, GridBase *grid = nullptr);
private:
    std::string  filename_{""}, dataname_{""};
    unsigned int nt_{0}, ni_{0}, nj_{0};
};

/******************************************************************************
 *                  Wrapper for A2A matrix block computation                  *
 ******************************************************************************/
template <typename T, typename Field, typename MetadataType, typename TIo = T>
class A2AMatrixBlockComputation
{
private:
    struct IoHelper
    {
        A2AMatrixIo<TIo> io;
        MetadataType     md;
        unsigned int     e, s, i, j;
    };
    typedef std::function<std::string(const unsigned int, const unsigned int)>  FilenameFn;
    typedef std::function<MetadataType(const unsigned int, const unsigned int)> MetadataFn;
    typedef std::function<void(int)> SwapFn;
public:
    // constructor
    A2AMatrixBlockComputation(GridBase *grid,
                              const unsigned int orthogDim,
                              const unsigned int next,
                              const unsigned int nstr,
                              const unsigned int blockSize,
                              const unsigned int cacheBlockSize,
                              TimerArray *tArray = nullptr);
    // execution
    void execute(const std::vector<Field> &left, 
                 const std::vector<Field> &right,
                 A2AKernel<T, Field> &kernel,
                 const FilenameFn &ionameFn,
                 const FilenameFn &filenameFn,
                 const MetadataFn &metadataFn,
                 std::vector<Field> *evecs = nullptr,
                 const std::vector<ComplexD> &evals = {},
                 const SwapFn *swapEvecCheckerFn = nullptr);
private:
    // I/O handler
    void saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h);
private:
    TimerArray            *tArray_;
    GridBase              *grid_;
    unsigned int          orthogDim_, nt_, next_, nstr_, blockSize_, cacheBlockSize_;
    Vector<T>             mCache_;
    Vector<Field>         lowBufi_, lowBufj_;
    Vector<TIo>           mBuf_;
    std::vector<IoHelper> nodeIo_;
};

/******************************************************************************
 *                       A2A matrix contraction kernels                       *
 ******************************************************************************/
class A2AContraction
{
public:
    // accTrMul(acc, a, b): acc += tr(a*b)
    template <typename C, typename MatLeft, typename MatRight>
    static inline void accTrMul(C &acc, const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;
        if ((MatLeft::Options  == RowMajor) and
            (MatRight::Options == ColMajor))
        {
      thread_for(r,a.rows(),
            {
                C tmp;
#ifdef USE_MKL
                dotuRow(tmp, r, a, b);
#else
                tmp = a.row(r).conjugate().dot(b.col(r));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
        else
      {
            thread_for(c,a.cols(),
            {
                C tmp;
#ifdef USE_MKL 
                dotuCol(tmp, c, a, b);
#else
                tmp = a.col(c).conjugate().dot(b.row(c));
#endif
                thread_critical
                {
                    acc += tmp;
                }
            });
        }
//         const int RowMajor = Eigen::RowMajor;
//         const int ColMajor = Eigen::ColMajor;
// #ifdef USE_MKL
//         if ((MatLeft::Options  == RowMajor) and
//             (MatRight::Options == ColMajor)) {
//       	  thread_for(r,a.rows(),
//                 {
//                     C tmp;
//                     dotuRow(tmp, r, a, b);
//                     thread_critical
//                     {
//                         acc += tmp;
//                     }
//                 });
//         } else {
//             thread_for(c,a.cols(),
//                 {
//                     C tmp;
//                     dotuCol(tmp, c, a, b);
//                     thread_critical
//                     {
//                         acc += tmp;
//                     }
//                 });
//         }
// #else
//         int aRows = a.rows();
//         int aCols = a.cols();
//         int bRows = b.rows();
//         int bCols = b.cols();
//         const C *aPt = a.data();
//         const C *bPt = b.data();

//         Vector<C> res(aRows,0.0);
//         C *resPt = &res[0];

//         if ((MatLeft::Options  == RowMajor) and
//             (MatRight::Options == ColMajor)) {
//             accelerator_for(r,aRows, 1,{
//                 C tmp = 0.0;
//                 for (int c = 0;c < aCols;c++) {
//                     tmp += innerProduct(*(aPt + r*aCols + c),*(bPt+r*bRows + c));
//                 }
//                 *(resPt+r) = tmp;
//             });

//         } else {
//             accelerator_for(c,aCols, 1,{
//                 C tmp = 0.0;
//                 for (int r = 0;r < aRows;r++) {
//                     tmp += innerProduct(*(aPt + c*aRows + r),*(bPt+c*bCols + r));
//                 }
//                 *(resPt+c) = tmp;
//             });


//         }
//         for (int r=0;r<aRows;r++) {
//             acc += res[r];
//         }
// #endif
    }

    template <typename MatLeft, typename MatRight>
    static inline double accTrMulFlops(const MatLeft &a, const MatRight &b)
    {
        double n = a.rows()*a.cols();

        return 8.*n;
    }

    // mul(res, a, b): res = a*b
#ifdef USE_MKL
    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexD, Opts...> &res, 
                           const Mat<ComplexD, Opts...> &a, 
                           const Mat<ComplexD, Opts...> &b)
    {
        static const ComplexD one(1., 0.), zero(0., 0.);
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexD, Opts...>::Options == RowMajor)
        {
            cblas_zgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexD, Opts...>::Options == ColMajor)
        {
            cblas_zgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }

    template <template <class, int...> class Mat, int... Opts>
    static inline void mul(Mat<ComplexF, Opts...> &res, 
                           const Mat<ComplexF, Opts...> &a, 
                           const Mat<ComplexF, Opts...> &b)
    {
        static const ComplexF one(1., 0.), zero(0., 0.);
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if ((res.rows() != a.rows()) or (res.cols() != b.cols()))
        {
            res.resize(a.rows(), b.cols());
        }
        if (Mat<ComplexF, Opts...>::Options == RowMajor)
        {
            cblas_cgemm(CblasRowMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.cols(), b.data(), b.cols(), &zero,
                        res.data(), res.cols());
        }
        else if (Mat<ComplexF, Opts...>::Options == ColMajor)
        {
            cblas_cgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, a.rows(), b.cols(),
                        a.cols(), &one, a.data(), a.rows(), b.data(), b.rows(), &zero,
                        res.data(), res.rows());
        }
    }
#else
    template <typename Mat>
    static inline void mul(Mat &res, const Mat &a, const Mat &b)
    {
        res = a*b;
    }
#endif
    template <typename Mat>
    static inline double mulFlops(const Mat &a, const Mat &b)
    {
        double nr = a.rows(), nc = a.cols();

        return nr*nr*(6.*nc + 2.*(nc - 1.));
    }
private:
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotRowPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aRow, 
                                    const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;

        if (MatLeft::Options == RowMajor)
        {
            aPt  = a.data() + aRow*a.cols();
            aInc = 1;
        }
        else if (MatLeft::Options == ColMajor)
        {
            aPt  = a.data() + aRow;
            aInc = a.rows();
        }
        if (MatRight::Options == RowMajor)
        {
            bPt  = b.data() + aRow;
            bInc = b.cols();
        }
        else if (MatRight::Options == ColMajor)
        {
            bPt  = b.data() + aRow*b.rows();
            bInc = 1;
        }
    }

#ifdef USE_MKL
    template <typename C, typename MatLeft, typename MatRight>
    static inline void makeDotColPt(C * &aPt, unsigned int &aInc, C * &bPt, 
                                    unsigned int &bInc, const unsigned int aCol, 
                                    const MatLeft &a, const MatRight &b)
    {
        const int RowMajor = Eigen::RowMajor;
        const int ColMajor = Eigen::ColMajor;
        if (MatLeft::Options == RowMajor)
        {
            aPt  = a.data() + aCol;
            aInc = a.cols();
        }
        else if (MatLeft::Options == ColMajor)
        {
            aPt  = a.data() + aCol*a.rows();
            aInc = 1;
        }
        if (MatRight::Options == RowMajor)
        {
            bPt  = b.data() + aCol*b.cols();
            bInc = 1;
        }
        else if (MatRight::Options == ColMajor)
        {
            bPt  = b.data() + aCol;
            bInc = b.rows();
        }
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexF &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_cdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexF &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexF *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_cdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuRow(ComplexD &res, const unsigned int aRow,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotRowPt(aPt, aInc, bPt, bInc, aRow, a, b);
        cblas_zdotu_sub(a.cols(), aPt, aInc, bPt, bInc, &res);
    }

    template <typename MatLeft, typename MatRight>
    static inline void dotuCol(ComplexD &res, const unsigned int aCol,
                               const MatLeft &a, const MatRight &b)
    {
        const ComplexD *aPt, *bPt;
        unsigned int   aInc, bInc;

        makeDotColPt(aPt, aInc, bPt, bInc, aCol, a, b);
        cblas_zdotu_sub(a.rows(), aPt, aInc, bPt, bInc, &res);
    }
#endif
};

/******************************************************************************
 *                     A2AMatrixIo template implementation                    *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T>
A2AMatrixIo<T>::A2AMatrixIo(std::string filename, std::string dataname, 
                            const unsigned int nt, const unsigned int ni,
                            const unsigned int nj)
: filename_(filename), dataname_(dataname)
, nt_(nt), ni_(ni), nj_(nj)
{}

// access //////////////////////////////////////////////////////////////////////
template <typename T>
unsigned int A2AMatrixIo<T>::getNt(void) const
{
    return nt_;
}

template <typename T>
unsigned int A2AMatrixIo<T>::getNi(void) const
{
    return ni_;
}

template <typename T>
unsigned int A2AMatrixIo<T>::getNj(void) const
{
    return nj_;
}

template <typename T>
size_t A2AMatrixIo<T>::getSize(void) const
{
    return nt_*ni_*nj_*sizeof(T);
}

// file allocation /////////////////////////////////////////////////////////////
template <typename T>
template <typename MetadataType>
void A2AMatrixIo<T>::initFile(const MetadataType &d, const unsigned int chunkSize)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t>    dim = {static_cast<hsize_t>(nt_), 
                                   static_cast<hsize_t>(ni_), 
                                   static_cast<hsize_t>(nj_)},
                            chunk = {static_cast<hsize_t>(nt_), 
                                     static_cast<hsize_t>(chunkSize), 
                                     static_cast<hsize_t>(chunkSize)};
    H5NS::DataSpace         dataspace(dim.size(), dim.data());
    H5NS::DataSet           dataset;
    H5NS::DSetCreatPropList plist;
    
    // create empty file just with metadata
    {
        Hdf5Writer writer(filename_);
        write(writer, dataname_, d);
    }

    // create the dataset
    Hdf5Reader reader(filename_, false);

    push(reader, dataname_);
    auto &group = reader.getGroup();
    plist.setChunk(chunk.size(), chunk.data());
    plist.setFletcher32();
    dataset = group.createDataSet(HADRONS_A2AM_NAME, Hdf5Type<T>::type(), dataspace, plist);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

// block I/O ///////////////////////////////////////////////////////////////////
template <typename T>
void A2AMatrixIo<T>::saveBlock(const T *data, 
                               const unsigned int i, 
                               const unsigned int j,
                               const unsigned int blockSizei,
                               const unsigned int blockSizej)
{
#ifdef HAVE_HDF5
    Hdf5Reader           reader(filename_, false);
    std::vector<hsize_t> count = {nt_, blockSizei, blockSizej},
                         offset = {0, static_cast<hsize_t>(i),
                                   static_cast<hsize_t>(j)},
                         stride = {1, 1, 1},
                         block  = {1, 1, 1}; 
    H5NS::DataSpace      memspace(count.size(), count.data()), dataspace;
    H5NS::DataSet        dataset;
    //    size_t               shift;

    push(reader, dataname_);
    auto &group = reader.getGroup();
    dataset     = group.openDataSet(HADRONS_A2AM_NAME);
    dataspace   = dataset.getSpace();
    dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                              stride.data(), block.data());
    dataset.write(data, Hdf5Type<T>::type(), memspace, dataspace);
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

template <typename T>
void A2AMatrixIo<T>::saveBlock(const A2AMatrixSet<T> &m,
                               const unsigned int ext, const unsigned int str,
                               const unsigned int i, const unsigned int j)
{
    unsigned int blockSizei = m.dimension(3);
    unsigned int blockSizej = m.dimension(4);
    unsigned int nstr       = m.dimension(1);
    size_t       offset     = (ext*nstr + str)*nt_*blockSizei*blockSizej;

    saveBlock(m.data() + offset, i, j, blockSizei, blockSizej);
}

template <typename T>
template <template <class> class Vec, typename VecT>
void A2AMatrixIo<T>::load(Vec<VecT> &v, double *tRead, GridBase *grid)
{
#ifdef HAVE_HDF5
    std::vector<hsize_t> hdim;
    H5NS::DataSet        dataset;
    H5NS::DataSpace      dataspace;
    H5NS::CompType       datatype;

    if (!(grid) || grid->IsBoss())
    {
        Hdf5Reader reader(filename_);
        push(reader, dataname_);
        auto &group = reader.getGroup();
        dataset = group.openDataSet(HADRONS_A2AM_NAME);
        datatype = dataset.getCompType();
        dataspace = dataset.getSpace();
        hdim.resize(dataspace.getSimpleExtentNdims());
        dataspace.getSimpleExtentDims(hdim.data());
        if ((nt_ * ni_ * nj_ != 0) and
            ((hdim[0] < nt_) or (hdim[1] < ni_) or (hdim[2] < nj_)))
        {
            HADRONS_ERROR(Size, "all-to-all matrix size mismatch (got "
                + std::to_string(hdim[0]) + "x" + std::to_string(hdim[1]) + "x"
                + std::to_string(hdim[2]) + ", expected "
                + std::to_string(nt_) + "x" + std::to_string(ni_) + "x"
                + std::to_string(nj_));
        }
        else if (ni_*nj_ == 0)
        {
            if (hdim[0] != nt_)
            {
                HADRONS_ERROR(Size, "all-to-all time size mismatch (got "
                    + std::to_string(hdim[0]) + ", expected "
                    + std::to_string(nt_) + ")");
            }
            ni_ = hdim[1];
            nj_ = hdim[2];
        }
    }
    if (grid)
    {
        grid->Broadcast(grid->BossRank(), &ni_, sizeof(unsigned int));
        grid->Broadcast(grid->BossRank(), &nj_, sizeof(unsigned int));
    }

    A2AMatrix<T>         buf(ni_, nj_);
    int broadcastSize =  sizeof(T) * buf.size();
    std::vector<hsize_t> count    = {1, static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)},
                         stride   = {1, 1, 1},
                         block    = {1, 1, 1},
                         memCount = {static_cast<hsize_t>(ni_),
                                     static_cast<hsize_t>(nj_)};
    H5NS::DataSpace      memspace(memCount.size(), memCount.data());

    std::cout << "Loading timeslice";
    std::cout.flush();
    *tRead = 0.;
    if (grid) {
        grid->Barrier();
        unsigned int myRank = grid->ThisRank(), nRank  = grid->RankCount();

        for (unsigned int tp1 = nt_-myRank; tp1 > 0; tp1-=nRank) {

            unsigned int         t      = tp1 - 1;
            std::vector<hsize_t> offset = {static_cast<hsize_t>(t), 0, 0};

            if (t % 10 == 0)
            {
                std::cout << " " << t;
                std::cout.flush();
            }
            if (grid->IsBoss())
            {
                dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                          stride.data(), block.data());
            }
            if (tRead) *tRead -= usecond();
            dataset.read(buf.data(), datatype, memspace, dataspace);

            v[t] = buf.template cast<VecT>();
            if (!grid->IsBoss()) {
                grid->SendToRecvFrom(buf.data(),myRank,(void *)&v[t],grid->BossRank(), broadcastSize);
            }

            if (tRead) *tRead += usecond();
        }
        grid->Barrier();
    } else {

        for (unsigned int tp1 = nt_; tp1 > 0; --tp1)
        {
            unsigned int         t      = tp1 - 1;
            std::vector<hsize_t> offset = {static_cast<hsize_t>(t), 0, 0};
            
            if (t % 10 == 0)
            {
                std::cout << " " << t;
                std::cout.flush();
            }
            dataspace.selectHyperslab(H5S_SELECT_SET, count.data(), offset.data(),
                                      stride.data(), block.data());
            if (tRead) *tRead -= usecond();
            dataset.read(buf.data(), datatype, memspace, dataspace);
            if (tRead) *tRead += usecond();
            v[t] = buf.template cast<VecT>();
        }
    }

    std::cout << std::endl;
#else
    HADRONS_ERROR(Implementation, "all-to-all matrix I/O needs HDF5 library");
#endif
}

/******************************************************************************
 *               A2AMatrixBlockComputation template implementation            *
 ******************************************************************************/
// constructor /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::A2AMatrixBlockComputation(GridBase *grid,
                            const unsigned int orthogDim,
                            const unsigned int next, 
                            const unsigned int nstr,
                            const unsigned int blockSize, 
                            const unsigned int cacheBlockSize,
                            TimerArray *tArray)
: grid_(grid), nt_(grid->GlobalDimensions()[orthogDim]), orthogDim_(orthogDim)
, next_(next), nstr_(nstr), blockSize_(blockSize), cacheBlockSize_(cacheBlockSize)
, tArray_(tArray)
{
    mCache_.resize(nt_*next_*nstr_*cacheBlockSize_*cacheBlockSize_);
    mBuf_.resize(nt_*next_*nstr_*blockSize_*blockSize_);
}

#define START_TIMER(name) if (tArray_) tArray_->startTimer(name)
#define STOP_TIMER(name)  if (tArray_) tArray_->stopTimer(name)
#define GET_TIMER(name)   ((tArray_ != nullptr) ? tArray_->getDTimer(name) : 0.)

// execution ///////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::execute(const std::vector<Field> &left, const std::vector<Field> &right,
          A2AKernel<T, Field> &kernel, const FilenameFn &ionameFn,
          const FilenameFn &filenameFn, const MetadataFn &metadataFn,
          std::vector<Field> *evecs,const std::vector<ComplexD> &evals, const SwapFn *swapEvecCheckerFn)
{
    //////////////////////////////////////////////////////////////////////////
    // i,j   is first  loop over blockSize_ factors
    // ii,jj is second loop over cacheBlockSize_ factors for high perf contractions
    // iii,jjj are loops within cacheBlock
    // Total index is sum of these  i+ii+iii etc...
    //////////////////////////////////////////////////////////////////////////

    bool checkerboarded_low = (swapEvecCheckerFn != nullptr);
    int Ncb = checkerboarded_low?2:1; // Ncb == 2 if the low modes are checkerboarded

    RealD norm = 1.0;

    int N_low = 0;
    if (evecs != nullptr) {
        norm = 1.0/::sqrt(norm2(evecs->at(0))); //Calculate norm of eigenvectors
        N_low = Ncb*evecs->size(); // N_low is the number of evecs for M + evecs for Mdag
    }

    int    N_i = left.size()+N_low; // Total number of bra vectors to contract
    int    N_j = right.size()+N_low; // Total number of ket vectors to contract

    double flops, bytes, t_kernel, t_gsum;
    double nodes = grid_->NodeCount();
    
    if (checkerboarded_low) {
        norm = norm/::sqrt(2); // Reduce checkerboarded norm to 1/sqrt(2)

        if (blockSize_%2 != 0 || cacheBlockSize_%2 != 0) {
            HADRONS_ERROR(Implementation, "Blocksize must be divisible by 2 for checkerboarded low modes");
        }

        lowBufi_.resize(blockSize_/2,evecs->at(0).Grid()); // storage for caching checkerboards
        lowBufj_.resize(blockSize_/2,evecs->at(0).Grid());
    }

    int NBlock_i = N_i/blockSize_ + (((N_i % blockSize_) != 0) ? 1 : 0); // Round up on the number of blocks to compute
    int NBlock_j = N_j/blockSize_ + (((N_j % blockSize_) != 0) ? 1 : 0);

    bool low_i, low_ii, low_j, low_jj;
    int i, j, evec_i, evec_j, N_ii, N_jj;

    i = 0, evec_i = 0;
    while (i < N_i) { // While we still have bra vectors to contract

        low_i = i < N_low;

        if (low_i) {
            N_ii = MIN(N_low-i,blockSize_);

            if (checkerboarded_low) {
                for (int idxi=evec_i;idxi<(MIN(N_low,i+N_ii)/2);idxi++) {
                    lowBufi_[idxi-evec_i] = evecs->at(idxi); // Cache original evecs to avoid excessive Meooe ops.
                    (*swapEvecCheckerFn)(idxi); // Swap original evec checkerboard to complementary checkerboard.
                }
            }
        } else {
            N_ii = MIN(N_i-i,blockSize_);
        }

        j = 0, evec_j = 0;
        while (j < N_j) {

            low_j = j < N_low;

            if (low_j) {
                N_jj = MIN(N_low-j,blockSize_);

                if (checkerboarded_low && i != j) { // Only cache and swap kets if it hasn't already been done for bras
                    for (int idxj=evec_j;idxj<(MIN(N_low,j+N_jj)/2);idxj++) {
                        lowBufj_[idxj-evec_j] = evecs->at(idxj);
                        (*swapEvecCheckerFn)(idxj);
                    }
                }
            } else {
                N_jj = MIN(N_j-j,blockSize_);
            }

            // Get the W and V vectors for this block^2 set of terms
            A2AMatrixSet<TIo> mBlock(mBuf_.data(), next_, nstr_, nt_, N_ii, N_jj);

            LOG(Message) << "All-to-all matrix block " 
                         << j/blockSize_ + NBlock_j*i/blockSize_ + 1 
                         << "/" << NBlock_i*NBlock_j << " [" << i <<" .. " 
                         << i+N_ii-1 << ", " << j <<" .. " << j+N_jj-1 << "]" 
                         << std::endl;
            // Series of cache blocked chunks of the contractions within this block
            flops    = 0.0;
            bytes    = 0.0;
            t_kernel = 0.0;

            double t;
            int ii, jj, evec_ii, evec_jj, N_iii, N_jjj;

            for (int cbi=0;cbi<Ncb;cbi++) { // For checkered low modes, loop through cacheBlock twice.
                
                // If we are only left with high modes (no checkerboard), continue.
                if (cbi && !low_i && !low_j)
                    continue;

                ii = 0, evec_ii = 0;
                while (ii < N_ii) {

                    low_ii = (i+ii) < N_low;

                    const Field *l_temp;
                    // If there are still low modes to process
                    if (low_ii) {
                        // Pick the min of how many low modes are left vs. cacheBlockSize_
                        N_iii = MIN(N_low-(i+ii),cacheBlockSize_);
                        if (cbi)
                            l_temp = &lowBufi_[0];
                        else
                            l_temp = &evecs->at(evec_i+evec_ii);
                    } else {
                        // Pick the min of how many high modes are left vs. cacheBlockSize_
                        N_iii = MIN(N_ii-ii,cacheBlockSize_);
                        l_temp = &left[i+ii-N_low];
                    }

                    jj = 0, evec_jj = 0;
                    while (jj < N_jj) {

                        low_jj = (j+jj) < N_low;

                        const Field *r_temp;
                        // If there are still low modes to process
                        if (low_jj) {
                            // Pick the min of how many low modes are left vs. cacheBlockSize_
                            N_jjj = MIN(N_low-(j+jj),cacheBlockSize_);
                            if (cbi) {
                                if (i == j)
                                    r_temp = &lowBufi_[0];
                                else
                                    r_temp = &lowBufj_[0];
                            } else {
                                r_temp = &evecs->at(evec_j+evec_jj);
                            }
                        } else {
                            // Pick the min of how many high modes are left vs. cacheBlockSize_
                            N_jjj = MIN(N_jj-jj,cacheBlockSize_);
                            r_temp = &right[j+jj-N_low];
                        }

                        A2AMatrixSet<T> mCacheBlock(mCache_.data(), next_, nstr_, nt_, N_iii, N_jjj);

                        START_TIMER("kernel");
                        kernel(mCacheBlock, l_temp, r_temp, orthogDim_, &t);
                        STOP_TIMER("kernel");
                        t_kernel += t;
                        flops    += kernel.flops(N_iii, N_jjj);
                        bytes    += kernel.bytes(N_iii, N_jjj);

                        START_TIMER("cache copy");
                        ComplexD coeff;
                        int evec_jjj=evec_j+evec_jj;
                        for(int jjj=0;jjj< N_jjj;jjj++) {

                            // If the ket vectors (corresponding to the solves) are low modes, multiply by the eigenvals
                            if (low_ii || low_jj) {
                                coeff = ComplexD(norm); // Normalize low modes appropriately
                                if (low_ii && low_jj) {
                                    coeff *= coeff;
                                }
                            } else {
                                coeff = ComplexD(1.0);
                            }

                            if (low_jj) {
                                if (jjj & 0x1) {
                                    coeff = coeff/conjugate(evals[evec_jjj]); // Mdaginv evals
                                    evec_jjj+=1;
                                } else {
                                    coeff = coeff/evals[evec_jjj]; // Minv evals
                                }
                            }

                            thread_for_collapse(4,e,next_,{
                              for(int iii = 0;iii < N_iii;iii++)
                              for(int s = 0;s < nstr_;s++)
                              for(int t = 0;t < nt_  ;t++) {
                                if (cbi)
                                    mBlock(e,s,t,ii+iii,jj+jjj) += coeff*mCacheBlock(e,s,t,iii,jjj);
                                else
                                    mBlock(e,s,t,ii+iii,jj+jjj)  = coeff*mCacheBlock(e,s,t,iii,jjj);
                            }});
                        }
                        STOP_TIMER("cache copy");

                        jj += N_jjj;
                        if (low_jj)
                            evec_jj += N_jjj/Ncb;
                    }

                    ii += N_iii;
                    if (low_ii)
                        evec_ii += N_iii/Ncb;
                }
            } // End for(cbi) loop

            if (checkerboarded_low) {
                t_gsum = -usecond();
                grid_->GlobalSumVector(&mBlock(0,0,0,0,0),mBlock.size());
                t_gsum += usecond();
            }
            // perf
            LOG(Message) << "Kernel perf " << flops/t_kernel/1.0e3/nodes 
                         << " Gflop/s/node " << std::endl;
            LOG(Message) << "Kernel Time: " << t_kernel << " us." << std::endl;
            LOG(Message) << "Global Sum Time: " << t_gsum << " us." << std::endl;

            // IO
            double       blockSize, ioTime;
            unsigned int myRank = grid_->ThisRank(), nRank  = grid_->RankCount();
        
            LOG(Message) << "Writing block to disk" << std::endl;
            ioTime = -GET_TIMER("IO: write block");
            START_TIMER("IO: total");
            makeFileDir(filenameFn(0, 0), grid_);

#ifdef HADRONS_A2AM_PARALLEL_IO
            grid_->Barrier();
            // make task list for current node
            nodeIo_.clear();
            for(int f = myRank; f < next_*nstr_; f += nRank)
            {
                IoHelper h;

                h.i  = i;
                h.j  = j;
                h.e  = f/nstr_;
                h.s  = f % nstr_;
                h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                        ionameFn(h.e, h.s), nt_, N_i, N_j);
                h.md = metadataFn(h.e, h.s);
                nodeIo_.push_back(h);
            }
            // parallel IO
            for (auto &h: nodeIo_)
            {
                saveBlock(mBlock, h);
            }
            grid_->Barrier();
#else
            // serial IO, for testing purposes only
            for(int e = 0; e < next_; e++)
            for(int s = 0; s < nstr_; s++)
            {
                IoHelper h;

                h.i  = i;
                h.j  = j;
                h.e  = e;
                h.s  = s;
                h.io = A2AMatrixIo<TIo>(filenameFn(h.e, h.s), 
                                        ionameFn(h.e, h.s), nt_, N_i, N_j);
                h.md = metadataFn(h.e, h.s);
                saveBlock(mfBlock, h);
            }
#endif
            STOP_TIMER("IO: total");
            blockSize  = static_cast<double>(next_*nstr_*nt_*N_ii*N_jj*sizeof(TIo));
            ioTime    += GET_TIMER("IO: write block");
            LOG(Message) << "HDF5 IO done " << sizeString(blockSize) << " in "
                         << ioTime  << " us (" 
                         << blockSize/ioTime*1.0e6/1024/1024
                         << " MB/s)" << std::endl;

            if (checkerboarded_low && low_j && i != j) {
                for (int idxj=evec_j;idxj<(MIN(N_low,j+N_jj)/2);idxj++)
                    evecs->at(idxj) = lowBufj_[idxj-evec_j];
            }
            j+=N_jj;
            evec_j+=(N_jj/Ncb);
        } // End while (j < N_j) Loop

        if (checkerboarded_low && low_i) {
            for (int idxi=evec_i;idxi<(MIN(N_low,i+N_ii)/2);idxi++)
                evecs->at(idxi) = lowBufi_[idxi-evec_i];
        }
        i+=N_ii;
        evec_i+=(N_ii/Ncb);
    } // End while (i < N_i) Loop
}

// I/O handler /////////////////////////////////////////////////////////////////
template <typename T, typename Field, typename MetadataType, typename TIo>
void A2AMatrixBlockComputation<T, Field, MetadataType, TIo>
::saveBlock(const A2AMatrixSet<TIo> &m, IoHelper &h)
{
    if ((h.i == 0) and (h.j == 0))
    {
        START_TIMER("IO: file creation");
        h.io.initFile(h.md, blockSize_);
        STOP_TIMER("IO: file creation");
    }
    START_TIMER("IO: write block");
    h.io.saveBlock(m, h.e, h.s, h.i, h.j);
    STOP_TIMER("IO: write block");
}

#undef START_TIMER
#undef STOP_TIMER
#undef GET_TIMER

END_HADRONS_NAMESPACE

#endif // A2A_Matrix_hpp_

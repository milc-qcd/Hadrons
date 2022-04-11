/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
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
#ifndef A2A_Vectors_hpp_
#define A2A_Vectors_hpp_

#include <Hadrons/Global.hpp>
#include <Hadrons/Environment.hpp>
#include <Hadrons/Solver.hpp>

BEGIN_HADRONS_NAMESPACE

/******************************************************************************
 *                 Class to generate V & W all-to-all vectors                 *
 ******************************************************************************/
template <typename FImpl>
class A2AVectorsSchur
{
public:
    FERM_TYPE_ALIASES(FImpl,);
    SOLVER_TYPE_ALIASES(FImpl,);
    HADRONS_DEFINE_SCHUR_OP(SchurOp,FImpl);
public:
    A2AVectorsSchur(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchur(void) = default;
    void makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d, const FermionField &evec, const Real &eval);

    void removeLowModeProj(std::vector<FermionField> &wout, const std::vector<FermionField> &evecs, const std::vector<ComplexD> evals);
    inline void makeLowModeCBeooe(FermionField &out, const FermionField &evec, const Complex eval);
    void makeLowModePairs(typename std::vector<FermionField>::iterator vecOut, 
                          const typename std::vector<FermionField>::iterator evec, const Complex eval);
    void makeLowModePairs(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator wout, 
                          const typename std::vector<FermionField>::iterator evec, const Complex eval);
    void makeLowModePairs(typename std::vector<FermionField>::iterator vecOut, typename std::vector<Complex>::iterator evalOut, 
                          const typename std::vector<FermionField>::iterator evec, const Complex eval);
    void makeLowModePairs5D(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator vout5,
                                                typename std::vector<FermionField>::iterator wout, typename std::vector<FermionField>::iterator wout5,
                                                const typename std::vector<FermionField>::iterator evec, const Complex eval);

    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d, 
                         const FermionField &noise_5d);
    void makeHighModeW(FermionField &wout, const FermionField &noise);
    void makeHighModeW(FermionField &wout, const FermionField &noise, std::vector<FermionField> &evecs, int size);
    void makeHighModeW5D(FermionField &vout_5d, FermionField &wout_5d, 
                         const FermionField &noise_5d);
protected:
    FMat                                     &action_;
    Solver                                   &solver_;
    GridBase                                 *frbGrid_, *gGrid_, *fGrid_;
    FermionField                             src_rb_, sol_rb1_, sol_rb2_, temp_, temp5_;

    SchurOp<FMat,FermionField> op_;
};

/******************************************************************************
 *                  Methods for V & W all-to-all vectors I/O                  *
 ******************************************************************************/
class A2AVectorsIo
{
public:
    struct Record: Serializable
    {
        GRID_SERIALIZABLE_CLASS_MEMBERS(Record,
                                        unsigned int, index);
        Record(void): index(0) {}
    };
public:
    template <typename Field>
    static void write(const std::string fileStem, std::vector<Field> &vec, 
                      const bool multiFile, const int trajectory = -1);
    template <typename Field>
    static void read(std::vector<Field> &vec, const std::string fileStem,
                     const bool multiFile, const int trajectory = -1);
private:
    static inline std::string vecFilename(const std::string stem, const int traj, 
                                          const bool multiFile)
    {
        std::string t = (traj < 0) ? "" : ("." + std::to_string(traj));

        if (multiFile)
        {
            return stem + t;
        }
        else
        {
            return stem + t + ".bin";
        }
    }
};

/******************************************************************************
 *            A2AVectorsSchur template implementation                  *
 ******************************************************************************/
template <typename FImpl>
A2AVectorsSchur<FImpl>::A2AVectorsSchur(FMat &action, Solver &solver)
: action_(action)
, solver_(solver)
, fGrid_(action_.FermionGrid())
, frbGrid_(action_.FermionRedBlackGrid())
, gGrid_(action_.GaugeGrid())
, src_rb_(frbGrid_)
, sol_rb1_(frbGrid_)
, sol_rb2_(frbGrid_)
, temp_(frbGrid_)
, temp5_(fGrid_)
, op_(action_)
{}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval)
{
    src_rb_ = evec;
    src_rb_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_rb1_, vout);
    pickCheckerboard(Odd, sol_rb2_, vout);

    /////////////////////////////////////////////////////
    // v_ie = -(1/eval_i) * MeeInv Meo MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_rb_, temp_);
    assert(temp_.Checkerboard() == Odd);
    action_.Meooe(temp_, sol_rb1_);
    assert(sol_rb1_.Checkerboard() == Even);
    action_.MooeeInv(sol_rb1_, temp_);
    assert(temp_.Checkerboard() == Even);
    sol_rb1_ = (-1.0 / eval) * temp_;
    assert(sol_rb1_.Checkerboard() == Even);

    /////////////////////////////////////////////////////
    // v_io = (1/eval_i) * MooInv evec_i
    /////////////////////////////////////////////////////
    action_.MooeeInv(src_rb_, temp_);
    assert(temp_.Checkerboard() == Odd);
    sol_rb2_ = (1.0 / eval) * temp_;
    assert(sol_rb2_.Checkerboard() == Odd);
    setCheckerboard(vout, sol_rb1_);
    assert(sol_rb1_.Checkerboard() == Even);
    setCheckerboard(vout, sol_rb2_);
    assert(sol_rb2_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, 
                                                    const FermionField &evec, const Real &eval)
{
    makeLowModeV(vout_5d, evec, eval);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval)
{
    src_rb_ = evec;
    src_rb_.Checkerboard() = Odd;
    pickCheckerboard(Even, sol_rb1_, wout);
    pickCheckerboard(Odd, sol_rb2_, wout);

    /////////////////////////////////////////////////////
    // w_ie = - MeeInvDag MoeDag Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_rb_, temp_);
    assert(temp_.Checkerboard() == Odd);
    action_.MeooeDag(temp_, sol_rb1_);
    assert(sol_rb1_.Checkerboard() == Even);
    action_.MooeeInvDag(sol_rb1_, temp_);
    assert(temp_.Checkerboard() == Even);
    sol_rb1_ = (-1.0) * temp_;

    /////////////////////////////////////////////////////
    // w_io = Doo evec_i
    /////////////////////////////////////////////////////
    op_.Mpc(src_rb_, sol_rb2_);
    assert(sol_rb2_.Checkerboard() == Odd);
    setCheckerboard(wout, sol_rb1_);
    assert(sol_rb1_.Checkerboard() == Even);
    setCheckerboard(wout, sol_rb2_);
    assert(sol_rb2_.Checkerboard() == Odd);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModeW5D(FermionField &wout_4d, 
                                                   FermionField &wout_5d, 
                                                   const FermionField &evec, 
                                                   const Real &eval)
{
    makeLowModeW(temp5_, evec, eval);
    action_.DminusDag(temp5_, wout_5d);
    action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
}

template <typename FImpl>
inline void A2AVectorsSchur<FImpl>::makeLowModeCBeooe(FermionField &out, const FermionField &evec, const Complex eval)
{
    int cb = evec.Checkerboard();
    int cbNeg = (cb==Even) ? Odd : Even;

    temp_ = Zero();
    out.Checkerboard() = cbNeg;
    temp_.Checkerboard() = cbNeg;
    action_.Meooe(evec, temp_);
    out = (1.0/eval) * temp_;
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::removeLowModeProj(std::vector<FermionField> &wout, const std::vector<FermionField> &evecs, const std::vector<ComplexD> evals)
{
    // Alternate form
    int cb = evecs[0].Checkerboard();
    int cbNeg = (cb==Even) ? Odd : Even;
    
    FermionField rbw(frbGrid_), rbwNeg(frbGrid_);

    rbw = Zero();
    rbwNeg = Zero();
    
    for (auto &w:wout) {
        rbw.Checkerboard() = cb;
        rbwNeg.Checkerboard() = cbNeg;

        pickCheckerboard(cb,rbw,w);
        pickCheckerboard(cbNeg,rbwNeg,w);

        // Add up W vector projection onto provided evec checkerboard
        temp_ = Zero();
        temp_.Checkerboard() = cb;
        for (int i=0;i<evecs.size();i++) {
          const FermionField& e = evecs[i];
          axpy(temp_,TensorRemove(innerProduct(e,rbw)),e,temp_);
        }
        // Subtract projected component from original. (factor of 2 compensates for normalization of checkerboard to 1/2)
        axpy(rbw,-2.0,temp_,rbw);
        setCheckerboard(w,rbw);

        
        action_.Meooe(rbwNeg, rbw); // Move cbNeg component of W to cb

        // Add up cbNeg checkerboard of W vector projection
        temp_ = Zero();
        temp_.Checkerboard() = cb;
        for (int i=0;i<evecs.size();i++) {
            RealD eval_Dinv = -1.0/pow(evals[i].imag(),2); // using Meooe twice brings two factors of 1/eval_D
            const FermionField& e = evecs[i];
            axpy(temp_,eval_Dinv*TensorRemove(innerProduct(e,rbw)),e,temp_);
        }
        rbw.Checkerboard() = cbNeg;
        action_.Meooe(temp_, rbw); // Move projection back to cbNeg checkerboard
        axpy(rbwNeg,-2.0,rbw,rbwNeg); // Subtract projected component from original. 
        setCheckerboard(w,rbwNeg);
    }

    // int cb = evecs[0].Checkerboard();
    // int cbNeg = (cb==Even) ? Odd : Even;
    
    // FermionField Mevec(fGrid_), Mdagevec(fGrid_), evecNeg(frbGrid_);

    // evecNeg.Checkerboard() = cbNeg;

    // for (int i=0;i<evecs.size();i++) {
    //     ComplexD eval_D = ComplexD(0.0,evals[i].imag());
    //     makeLowModeCBeooe(evecNeg,evecs[i],eval_D);

    //     setCheckerboard(Mevec,evecNeg);
    //     setCheckerboard(Mevec,evecs[i]);

    //     if (cb == Even) {
    //         temp_.Checkerboard() = cbNeg;
    //         temp_ = -evecNeg;
    //         setCheckerboard(Mdagevec,temp_);
    //         setCheckerboard(Mdagevec,evecs[i]);
    //     } else {
    //         temp_.Checkerboard() = cb;
    //         temp_ = -evecs[i];
    //         setCheckerboard(Mdagevec,temp_);
    //         setCheckerboard(Mdagevec,evecNeg);
    //     }

    //     for (auto &w:wout) {
    //         auto ip = innerProduct(Mevec,w);
    //         w = w - ip*Mevec;
    //         ip = innerProduct(Mdagevec,w);
    //         w = w - ip*Mdagevec;
    //     }
    // }    
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModePairs(typename std::vector<FermionField>::iterator vecOut, 
                                              const typename std::vector<FermionField>::iterator evec, const Complex eval)
{
    double cbEven = (*evec).Checkerboard() == Even;
    int cbParity = cbEven ? Even : Odd;
    int cbParityNeg = !cbEven ? Even : Odd;


    //Expects eigenvalues of M
    ComplexD eval_D = ComplexD(0.0,eval.imag());

    makeLowModeCBeooe(sol_rb1_,evec,eval_D);

    setCheckerboard(*vecOut, sol_rb1_);
    setCheckerboard(*vecOut, *evec);

    if (cbEven){
        pickCheckerboard(cbParityNeg, temp_, *(vecOut+1));
        temp_ = -1.0 * sol_rb1_;

        setCheckerboard(*(vecOut+1), temp_);
        setCheckerboard(*(vecOut+1), *evec);
    } else {
        pickCheckerboard(cbParity, temp_, *(vecOut+1));
        temp_ = -1.0 * (*evec);

        setCheckerboard(*(vecOut+1), temp_);
        setCheckerboard(*(vecOut+1), sol_rb1_);
    }
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModePairs(typename std::vector<FermionField>::iterator vecOut, typename std::vector<Complex>::iterator evalOut, 
                                              const typename std::vector<FermionField>::iterator evec, const Complex eval)
{
    makeLowModePairs(vecOut,evec,eval);

    *evalOut = 1.0/eval;
    *(evalOut+1) = 1.0/conjugate(eval);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModePairs(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator wout, 
                                              const typename std::vector<FermionField>::iterator evec, const Complex eval)
{
    std::vector<ComplexD> evals(2);

    makeLowModePairs(wout,evals.begin(),evec,eval);

    *vout = evals[0]*(*wout);
    *(vout+1) = evals[1]*(*(wout+1));

}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModePairs5D(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator vout5,
                                                typename std::vector<FermionField>::iterator wout, typename std::vector<FermionField>::iterator wout5,
                                                const typename std::vector<FermionField>::iterator evec, const Complex eval)
{
    makeLowModePairs(vout5,wout5, evec, eval);
    action_.ExportPhysicalFermionSolution(*vout5, *vout);
    action_.ExportPhysicalFermionSolution(*(vout5+1), *(vout+1));
    action_.DminusDag(temp5_, *wout5);
    action_.ExportPhysicalFermionSolution(temp5_, *wout);
    action_.DminusDag(temp5_, *(wout5+1));
    action_.ExportPhysicalFermionSolution(temp5_, *(wout+1));
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeHighModeV(FermionField &vout, 
                                                  const FermionField &noise)
{
    solver_(vout, noise);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeHighModeV5D(FermionField &vout_4d, 
                                                    FermionField &vout_5d, 
                                                    const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportPhysicalFermionSource(noise, temp5_);
    }
    else
    {
        temp5_ = noise;
    }
    makeHighModeV(vout_5d, temp5_);
    action_.ExportPhysicalFermionSolution(vout_5d, vout_4d);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeHighModeW(FermionField &wout, 
                                                  const FermionField &noise)
{
    wout = noise;
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeHighModeW(FermionField &wout, const FermionField &noise,
                                            std::vector<FermionField> &evecs, int size)
{
    wout = noise;
    basisOrthogonalize(evecs, wout, size);
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeHighModeW5D(FermionField &wout_4d, 
                                                    FermionField &wout_5d, 
                                                    const FermionField &noise)
{
    if (noise.Grid()->Dimensions() == fGrid_->Dimensions() - 1)
    {
        action_.ImportUnphysicalFermion(noise, wout_5d);
        wout_4d = noise;
    }
    else
    {
        wout_5d = noise;
        action_.ExportPhysicalFermionSource(wout_5d, wout_4d);
    }
}

/******************************************************************************
 *               all-to-all vectors I/O template implementation               *
 ******************************************************************************/
template <typename Field>
void A2AVectorsIo::write(const std::string fileStem, std::vector<Field> &vec, 
                         const bool multiFile, const int trajectory)
{
    Record       record;
    GridBase     *grid = vec[0].Grid();
    ScidacWriter binWriter(grid->IsBoss());
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Writing vector " << i << std::endl;
            makeFileDir(fullFilename, grid);
            binWriter.open(fullFilename);
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
            binWriter.close();
        }
    }
    else
    {
        makeFileDir(filename, grid);
        binWriter.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Writing vector " << i << std::endl;
            record.index = i;
            binWriter.writeScidacFieldRecord(vec[i], record);
        }
        binWriter.close();
    }
}

template <typename Field>
void A2AVectorsIo::read(std::vector<Field> &vec, const std::string fileStem, 
                        const bool multiFile, const int trajectory)
{
    Record       record;
    ScidacReader binReader;
    std::string  filename = vecFilename(fileStem, trajectory, multiFile);

    if (multiFile)
    {
        std::string fullFilename;

        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            fullFilename = filename + "/elem" + std::to_string(i) + ".bin";

            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.open(fullFilename);
            binReader.readScidacFieldRecord(vec[i], record);
            binReader.close();
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
    }
    else
    {
        binReader.open(filename);
        for (unsigned int i = 0; i < vec.size(); ++i)
        {
            LOG(Message) << "Reading vector " << i << std::endl;
            binReader.readScidacFieldRecord(vec[i], record);
            if (record.index != i)
            {
                HADRONS_ERROR(Io, "vector index mismatch");
            }
        }
        binReader.close();
    }
}

END_HADRONS_NAMESPACE

#endif // A2A_Vectors_hpp_

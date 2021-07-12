/*
 * A2AVectors.hpp, part of Hadrons (https://github.com/aportelli/Hadrons)
 *
 * Copyright (C) 2015 - 2020
 *
 * Author: Antonin Portelli <antonin.portelli@me.com>
 * Author: Peter Boyle <paboyle@ph.ed.ac.uk>
 * Author: fionnoh <fionnoh@gmail.com>
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
    // using SchurOp = typename std::conditional<HADRONS_IS_STAGGERED_IMPLEMENTATION(FImpl), HADRONS_DEFAULT_SCHUR_OP_STAGGERED<FMat,FermionField>, HADRONS_DEFAULT_SCHUR_OP<FMat,FermionField> >::type;
public:
    A2AVectorsSchur(FMat &action, Solver &solver);
    virtual ~A2AVectorsSchur(void) = default;
    void makeLowModeV(FermionField &vout, const FermionField &evec, const Real &eval);
    void makeLowModeV5D(FermionField &vout_4d, FermionField &vout_5d, const FermionField &evec, const Real &eval);
    void makeLowModeW(FermionField &wout, const FermionField &evec, const Real &eval);
    void makeLowModeW5D(FermionField &wout_4d, FermionField &wout_5d, const FermionField &evec, const Real &eval);

    void makeLowModePairs(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator wout, 
                          const typename std::vector<FermionField>::iterator evec, const Real mass, const Real eval, bool cbEven = true);
    void makeLowModePairs5D(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator vout5,
                                                typename std::vector<FermionField>::iterator wout, typename std::vector<FermionField>::iterator wout5,
                                                const typename std::vector<FermionField>::iterator evec, const Real mass, const Real eval, bool cbEven = true);

    void makeHighModeV(FermionField &vout, const FermionField &noise);
    void makeHighModeV5D(FermionField &vout_4d, FermionField &vout_5d, 
                         const FermionField &noise_5d);
    void makeHighModeW(FermionField &wout, const FermionField &noise);
    void makeHighModeW5D(FermionField &vout_5d, FermionField &wout_5d, 
                         const FermionField &noise_5d);
public:
    template <typename T = FImpl>
    typename std::enable_if<HADRONS_IS_STAGGERED_IMPLEMENTATION(T),bool>::type isStaggered(){ return true; }
    template <typename T = FImpl>
    typename std::enable_if<!HADRONS_IS_STAGGERED_IMPLEMENTATION(T),bool>::type isStaggered(){ return false; }
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
void A2AVectorsSchur<FImpl>::makeLowModePairs(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator wout, 
                                              const typename std::vector<FermionField>::iterator evec, const Real mass, const Real eval, bool cbEven)
{
    int cbParity = cbEven ? Even : Odd;
    int cbParityNeg = !cbEven ? Even : Odd;

    //Expects eigenvalues of Dslash squarred
    ComplexD eval_D = ComplexD(0,sqrt(eval));

    src_rb_ = *evec;
    src_rb_.Checkerboard() = cbParity;
    pickCheckerboard(cbParityNeg, sol_rb1_, *wout);
    
    action_.Meooe(src_rb_, temp_);
    sol_rb1_ = (1.0/eval_D) * temp_;

    setCheckerboard(*wout, sol_rb1_);
    setCheckerboard(*wout, src_rb_);

    if (cbEven){
        pickCheckerboard(cbParityNeg, temp_, *(wout+1));
        temp_ = -1.0 * sol_rb1_;

        setCheckerboard(*(wout+1), temp_);
        setCheckerboard(*(wout+1), src_rb_);
    } else {
        pickCheckerboard(cbParity, temp_, *(wout+1));
        temp_ = -1.0 * src_rb_;

        setCheckerboard(*(wout+1), temp_);
        setCheckerboard(*(wout+1), sol_rb1_);
    }
    *vout = (1.0/(mass+eval_D))*(*wout);
    *(vout+1) = (1.0/(mass-eval_D))*(*(wout+1));
}

template <typename FImpl>
void A2AVectorsSchur<FImpl>::makeLowModePairs5D(typename std::vector<FermionField>::iterator vout, typename std::vector<FermionField>::iterator vout5,
                                                typename std::vector<FermionField>::iterator wout, typename std::vector<FermionField>::iterator wout5,
                                                const typename std::vector<FermionField>::iterator evec, const Real mass, const Real eval, bool cbEven)
{
    makeLowModePairs(vout5,wout5, evec, mass, eval, cbEven);
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

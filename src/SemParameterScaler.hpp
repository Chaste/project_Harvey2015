/*

Copyright (c) 2005-2015, University of Oxford.
All rights reserved.

University of Oxford means the Chancellor, Masters and Scholars of the
University of Oxford, having an administrative office at Wellington
Square, Oxford OX1 2JD, UK.

This file is part of Chaste.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:
 * Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.
 * Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.
 * Neither the name of the University of Oxford nor the names of its
   contributors may be used to endorse or promote products derived from this
   software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.

*/

#ifndef SEMPARAMETERSCALER_HPP_
#define SEMPARAMETERSCALER_HPP_

#include "ChasteSerialization.hpp"
#include "SerializableSingleton.hpp"

#include <map>

/**
 * @class SemParameterScaler
 *
 * A class for uniformly scaling parameters for a cell-centre simulation.
 * This class is instantiated using the singleton pattern to ensure
 * that parameters are scaled consistently throughout a simulation.
 */
template<unsigned DIM>
class SemParameterScaler
{
private:

    /**A pointer to the singleton instance of this class. */
    static SemParameterScaler<DIM>* mpInstance;

    /** Whether we have set the number of elements. */
    bool mNumElementsAreSet;

    /** The number of elements. */
    unsigned mNumSemElements;

    /** Map to the number of elements. */
    std::map<unsigned, unsigned> mNumScaledSemElements;

    /** Scaling factors for the interaction force. */
    double mSpringConstantScalingFactor;
    double mNewmanScalingFactor;

    /** The equilibrium distance of nodes in the mesh. */
    double mRestLength;

    /** Test class is a friend to enable better testing. */
    friend class TestSemParameterScaler;

protected:

    /** 
     * Protected default constructor. See "Singleton Pattern",
     * Design Patterns, Gamma et al.
     */
    SemParameterScaler();

public:

    /**
     * Access the single instance of this class.
     */
    static SemParameterScaler<DIM>* Instance();

    /**
     * Teardown the single instance of this class.
     */
    static void Destroy();

    /**
     * Set the private members of this class.
     */
    void SetNumElementsPerCell(unsigned numElements);
    void SetNumElementsForScaledCell(unsigned semIndex, unsigned numElements);
    void SetSpringConstantScalingFactor(double scalingFactor);
    void SetNewmanScalingFactor(double scalingFactor);

    /**
     * Get the public members of this class.
     */
    unsigned GetNumElements(unsigned semIndex);
    double GetRestLength();
    double GetNewmanScalingFactor();

    /**
     * Scale a cell damping constant.
     */
    double ScaleDampingConstant(unsigned semIndex, double unscaledDampingConstant);
    double ScaleSpringConstant(unsigned semIndex, double unscaledSpringConstant);
};

#endif /*SEMPARAMETERSCALER_HPP_*/

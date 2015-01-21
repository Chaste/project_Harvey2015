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

#include "SemParameterScaler.hpp"
#include "Exception.hpp"
#include <assert.h>
#include <cmath>

template<unsigned DIM>
SemParameterScaler<DIM>* SemParameterScaler<DIM>::mpInstance = NULL;

template<unsigned DIM>
SemParameterScaler<DIM>* SemParameterScaler<DIM>::Instance()
{
    if (mpInstance == NULL)
    {
        mpInstance = new SemParameterScaler<DIM>;
        std::atexit(Destroy);
    }
    return mpInstance;
}

template<unsigned DIM>
void SemParameterScaler<DIM>::Destroy()
{
    if (mpInstance)
    {
        delete mpInstance;
        mpInstance = NULL;
    }
}

template<unsigned DIM>
SemParameterScaler<DIM>::SemParameterScaler()
    :   mNumElementsAreSet(false),
        mSpringConstantScalingFactor(0.5),
        mNewmanScalingFactor(2.0)
{
    assert(mpInstance == NULL);
}

template<unsigned DIM>
void SemParameterScaler<DIM>::SetNumElementsPerCell(unsigned numElements)
{
    if (mNumElementsAreSet)
    {
        EXCEPTION("Num SEM elements per cell can only be set once");
    }

    assert(numElements > 0);

    mNumSemElements = numElements;

    switch (DIM)
    {
        case 1:
        {
            mRestLength = 1.0 / ((double)(mNumSemElements) - 1.0);
            break;
        }
        case 2:
        {
            double packing_density_2d = M_PI / sqrt(12.0);
            mRestLength = sqrt(packing_density_2d / (double)mNumSemElements);
            break;
        }
        case 3:
        {
            double packing_density_3d = M_PI / sqrt(18.0);
            mRestLength = pow(packing_density_3d / (double)mNumSemElements, 1.0/3.0);
            break;
        }
    }
    mNumElementsAreSet = true;
}

template<unsigned DIM>
void SemParameterScaler<DIM>::SetNumElementsForScaledCell(unsigned semIndex, unsigned numElements)
{
    mNumScaledSemElements[semIndex] = numElements;
}

template<unsigned DIM>
void SemParameterScaler<DIM>::SetSpringConstantScalingFactor(double scalingFactor)
{
    assert(scalingFactor > 0.0 && scalingFactor < 1.0);

    mSpringConstantScalingFactor = scalingFactor;
}

template<unsigned DIM>
void SemParameterScaler<DIM>::SetNewmanScalingFactor(double scalingFactor)
{
    assert(scalingFactor > 0.0);

    mNewmanScalingFactor = scalingFactor;
}

template<unsigned DIM>
double SemParameterScaler<DIM>::GetRestLength()
{
    assert(mNumElementsAreSet);
    return mRestLength;
}

template<unsigned DIM>
unsigned SemParameterScaler<DIM>::GetNumElements(unsigned semIndex)
{
    unsigned num_elements;

    if (mNumScaledSemElements.find(semIndex) == mNumScaledSemElements.end())
    {
        num_elements = mNumSemElements;
    }
    else
    {
        num_elements = mNumScaledSemElements[semIndex];
    }

    return num_elements;
}

template<unsigned DIM>
double SemParameterScaler<DIM>::GetNewmanScalingFactor()
{
    assert(mNumElementsAreSet);
    return mNewmanScalingFactor;
}

template<unsigned DIM>
double SemParameterScaler<DIM>::ScaleDampingConstant(unsigned semIndex, double unscaledDampingConstant)
{
    assert(mNumElementsAreSet);
    return unscaledDampingConstant / (double)(GetNumElements(semIndex));
}

template<unsigned DIM>
double SemParameterScaler<DIM>::ScaleSpringConstant(unsigned semIndex, double unscaledSpringConstant)
{
    assert(mNumElementsAreSet);

    double scaling_factor = pow( GetRestLength() / mNewmanScalingFactor , 2.0) / 8.0;

    switch (DIM)
    {
        case 1:
        {
            scaling_factor *= (double) (GetNumElements(semIndex));
            break;
        }
        case 2:
        {
            scaling_factor *= (1.0 - mSpringConstantScalingFactor / sqrt((double)(GetNumElements(semIndex))));
            break;
        }
        case 3:
        {
            scaling_factor *= pow((double)(GetNumElements(semIndex)), -1.0 / 3.0) * (1.0 - mSpringConstantScalingFactor * pow((double)(GetNumElements(semIndex)), -1.0 / 3.0));
            break;
        }
    }

    return scaling_factor * unscaledSpringConstant;
}

// Explicit instantiation
template class SemParameterScaler<1>;
template class SemParameterScaler<2>;
template class SemParameterScaler<3>;

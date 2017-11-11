//=================================================================================================
/*!
//  \file blazetest/utiltest/AlignedResource.h
//  \brief Header file for the AlignedResource class template
//
//  Copyright (C) 2012-2017 Klaus Iglberger - All Rights Reserved
//
//  This file is part of the Blaze library. You can redistribute it and/or modify it under
//  the terms of the New (Revised) BSD License. Redistribution and use in source and binary
//  forms, with or without modification, are permitted provided that the following conditions
//  are met:
//
//  1. Redistributions of source code must retain the above copyright notice, this list of
//     conditions and the following disclaimer.
//  2. Redistributions in binary form must reproduce the above copyright notice, this list
//     of conditions and the following disclaimer in the documentation and/or other materials
//     provided with the distribution.
//  3. Neither the names of the Blaze development group nor the names of its contributors
//     may be used to endorse or promote products derived from this software without specific
//     prior written permission.
//
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY
//  EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
//  OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT
//  SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//  TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR
//  BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
//  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN
//  ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
//  DAMAGE.
*/
//=================================================================================================

#ifndef _BLAZETEST_UTILTEST_ALIGNEDRESOURCE_H_
#define _BLAZETEST_UTILTEST_ALIGNEDRESOURCE_H_


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <blaze/util/StaticAssert.h>
#include <blaze/util/typetraits/AlignmentOf.h>
#include <blazetest/utiltest/InstanceCounter.h>


namespace blazetest {

namespace utiltest {

//=================================================================================================
//
//  CLASS DEFINITION
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Implementation of an instance counted, aligned resource.
//
// The AlignedResource class represents an important resource for testing purposes. It is instance
// counted via the InstanceCounter class and guaranteed to be 64-bit aligned. Additionally, it
// provides a data member that is guaranteed to be initialized to 7 via the default constructor.
*/
class alignas( 64UL ) AlignedResource : public InstanceCounter<AlignedResource>
{
 public:
   //**Constructors********************************************************************************
   /*!\name Constructors */
   //@{
   inline AlignedResource();
   //@}
   //**********************************************************************************************

   //**Utility functions***************************************************************************
   /*!\name Utility functions */
   //@{
   inline unsigned int getValue();
   //@}
   //**********************************************************************************************

 private:
   //**Member variables****************************************************************************
   /*!\name Member variables */
   //@{
   unsigned int value_;  // The value of the aligned resource.
   //@}
   //**********************************************************************************************
};
//*************************************************************************************************




//=================================================================================================
//
//  CONSTRUCTOR
//
//=================================================================================================

//*************************************************************************************************
/*!\brief The constructor of AlignedResource.
*/
inline AlignedResource::AlignedResource()
   : InstanceCounter<AlignedResource>()
   , value_( 7 )
{
   BLAZE_STATIC_ASSERT( blaze::AlignmentOf<AlignedResource>::value == 64UL );
}
//*************************************************************************************************




//=================================================================================================
//
//  UTILITY FUNCTIONS
//
//=================================================================================================

//*************************************************************************************************
/*!\brief Returns the current value of the resource.
//
// \return The current value of the resource.
*/
inline unsigned int AlignedResource::getValue()
{
   return value_;
}
//*************************************************************************************************

} // namespace utiltest

} // namespace blazetest

#endif

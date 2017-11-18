//=================================================================================================
/*!
//  \file src/mathtest/svecdvecouter/VCaVDb.cpp
//  \brief Source file for the VCaVDb sparse vector/dense vector outer product math test
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


//*************************************************************************************************
// Includes
//*************************************************************************************************

#include <cstdlib>
#include <iostream>
#include <blaze/math/CompressedVector.h>
#include <blaze/math/DynamicVector.h>
#include <blazetest/mathtest/Creator.h>
#include <blazetest/mathtest/svecdvecouter/OperationTest.h>
#include <blazetest/system/MathTest.h>


//=================================================================================================
//
//  MAIN FUNCTION
//
//=================================================================================================

//*************************************************************************************************
int main()
{
   std::cout << "   Running 'VCaVDb'..." << std::endl;

   using blazetest::mathtest::TypeA;
   using blazetest::mathtest::TypeB;

   try
   {
      // Vector type definitions
      typedef blaze::CompressedVector<TypeA>  VCa;
      typedef blaze::DynamicVector<TypeB>     VDb;

      // Creator type definitions
      typedef blazetest::Creator<VCa>  CVCa;
      typedef blazetest::Creator<VDb>  CVDb;

      // Running tests with small vectors
      for( size_t i=0UL; i<=8UL; ++i ) {
         for( size_t j=0UL; j<=8UL; ++j ) {
            for( size_t k=0UL; k<=i; ++k ) {
               RUN_SVECDVECOUTER_OPERATION_TEST( CVCa( i, k ), CVDb( j ) );
            }
         }
      }

      // Running tests with large vectors
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa(  67UL,  7UL ), CVDb(  67UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa(  67UL, 13UL ), CVDb( 127UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa( 127UL,  7UL ), CVDb(  67UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa( 127UL, 13UL ), CVDb( 127UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa(  64UL,  8UL ), CVDb(  64UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa(  64UL, 16UL ), CVDb( 128UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa( 128UL,  8UL ), CVDb(  64UL ) );
      RUN_SVECDVECOUTER_OPERATION_TEST( CVCa( 128UL, 16UL ), CVDb( 128UL ) );
   }
   catch( std::exception& ex ) {
      std::cerr << "\n\n ERROR DETECTED during sparse vector/dense vector outer product:\n"
                << ex.what() << "\n";
      return EXIT_FAILURE;
   }

   return EXIT_SUCCESS;
}
//*************************************************************************************************
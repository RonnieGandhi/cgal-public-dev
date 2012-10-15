// TODO: Add licence
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL:$
// $Id: $
// 
//
// Author(s)     :  
//
// ============================================================================


/*! \file Extended_descartes.C
 This is the test file for the class Extended_descartes.h
*/

#include <CGAL/basic.h>

#ifndef CGAL_ISOLATOR_USES_INPUT_RANGE
#define CGAL_ISOLATOR_USES_INPUT_RANGE 1
#endif

#ifndef CGAL_DESCARTES_VERBOSE
#define CGAL_DESCARTES_VERBOSE 1
#endif

#ifndef CGAL_DESCARTES_EXTENDED_VERBOSE
#define CGAL_DESCARTES_EXTENDED_VERBOSE 1
#endif

// include these traits here by 'hand', since not in release 3.3
#include <CGAL/Algebraic_extension_traits.h>
#include <CGAL/Scalar_factor_traits.h>

#include <CGAL/Polynomial_type_generator.h>
#include <_test_real_root_isolator.h>

#include <CGAL/Arithmetic_kernel.h>
#include <CGAL/Algebraic_kernel_d/Descartes_star.h>

template <class AT>
void test_descartes_star(){
  typedef typename AT::Integer Integer;
  typedef typename AT::Rational Rational;
  { 
    // general test of concept RealRootIsolator for Integer coefficients
    typedef ::CGAL::internal::Descartes_star< Integer, Rational > Isolator;
    CGAL::internal::test_real_root_isolator<Isolator>();
  }{
    // general test of concept RealRootIsolator for Rational coefficients
    typedef ::CGAL::internal::Descartes_star< Rational, Rational > Isolator;
    CGAL::internal::test_real_root_isolator<Isolator>();
  }    
}

template < class ArithmeticKernel >
void simple_test_descartes_star(){
  
  typedef ArithmeticKernel Arithmetic_kernel;
  typedef typename Arithmetic_kernel::Rational Rational;
  typedef typename Arithmetic_kernel::Integer Integer;
/*
  std::cout << "Rational polynomial:" << std::endl;

  std::list<Rational> coeffs;
  
//  coeffs.push_back(-100);
//  coeffs.push_back(0);
//  coeffs.push_back(80);
//  coeffs.push_back(0);
//  coeffs.push_back(-17);
//  coeffs.push_back(0);
//  coeffs.push_back(1);

//  coeffs.push_back(0);
//  coeffs.push_back(Rational(1)/Rational(1024));
//  coeffs.push_back(0);
//  coeffs.push_back(Rational(-5)/Rational(64));
//  coeffs.push_back(0);
//  coeffs.push_back(1);

//  coeffs.push_back(-1);
//  coeffs.push_back(0);
//  coeffs.push_back(25);

	coeffs.push_back(Rational(-1)/Rational(36));
  coeffs.push_back(0);
  coeffs.push_back(1);
  
  typedef ::CGAL::internal::Descartes_star< Rational, Rational > Isolator;
  
  Isolator isolator(coeffs.begin(), coeffs.end());
  std::cout << "#Roots: " <<  isolator.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator.left_bound(0) << "," 
            << isolator.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator.left_bound(1) << "," 
            << isolator.right_bound(1) << "]" << std::endl;
//  std::cout << "3rd root: [" << isolator.left_bound(2) << "," 
//            << isolator.right_bound(2) << "]" << std::endl;
//  std::cout << "4th root: [" << isolator.left_bound(3) << "," 
//            << isolator.right_bound(3) << "]" << std::endl;
//  std::cout << "5th root: [" << isolator.left_bound(4) << "," 
//            << isolator.right_bound(4) << "]" << std::endl;  
//  std::cout << "6th root: [" << isolator.left_bound(5) << "," 
//            << isolator.right_bound(5) << "]" << std::endl;  
  */

  std::cout << "Einfaches Integer-Polynom:" << std::endl;
  std::list<Integer> coeffs_integer;
  coeffs_integer.push_back(2000);
  coeffs_integer.push_back(-3000);
  coeffs_integer.push_back(1000);

  typedef ::CGAL::internal::Descartes_star< Integer, Rational > Isolator_integer;
  Isolator_integer isolator_integer(coeffs_integer.begin(), coeffs_integer.end());

  std::cout << "#Roots: " <<  isolator_integer.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator_integer.left_bound(0) << "," 
            << isolator_integer.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator_integer.left_bound(1) << "," 
            << isolator_integer.right_bound(1) << "]" << std::endl;
  
/*

  std::cout << "Mignotte-like polynomial:" << std::endl;

  std::list<Integer> coeff_mig;
  coeff_mig.push_back(-2);
  coeff_mig.push_back(800);
  coeff_mig.push_back(-80000);
  for (int i = 0; i < 57; i++) {
    coeff_mig.push_back(0);
  }
  coeff_mig.push_back(1048576);
  
  typedef ::CGAL::internal::Descartes_star< Integer, Rational > Isolator_integer;
  
  Isolator_integer isolator_mignotte(coeff_mig.begin(), coeff_mig.end());
  std::cout << "#Roots: " <<  isolator_mignotte.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator_mignotte.left_bound(0) << "," 
            << isolator_mignotte.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator_mignotte.left_bound(1) << "," 
            << isolator_mignotte.right_bound(1) << "]" << std::endl;
  std::cout << "3nd root: [" << isolator_mignotte.left_bound(2) << "," 
            << isolator_mignotte.right_bound(2) << "]" << std::endl;
  std::cout << "4nd root: [" << isolator_mignotte.left_bound(3) << "," 
            << isolator_mignotte.right_bound(3) << "]" << std::endl;
*/
/*
  std::cout << "Algebraic polynomial:" << std::endl;

  typedef typename CGAL::Algebraic_kernel_d_1_generator< Integer, Rational >::
    Algebraic_kernel_with_qir_and_descartes_1  Algebraic_kernel_1;
  typedef typename Algebraic_kernel_1::Algebraic_real_1 Algebraic_real_1;
  typedef std::list< std::pair< Algebraic_real_1, unsigned int > > Roots;

  std::list<Algebraic_real_1> coeff_real;
  Roots roots;
  typename Algebraic_kernel_1::Solve_1 solve;

  typedef typename Algebraic_kernel_1::Polynomial_1 Polynomial_1;
  typedef CGAL::Polynomial_traits_d<Polynomial_1> PT_1;

  Polynomial_1 x = typename PT_1::Shift()(Polynomial_1(1),1);
  Polynomial_1 polynomial = x*x - 2;
  solve(polynomial, std::back_inserter(roots));

  Algebraic_real_1 a = roots.back().first;    
  coeff_real.push_back(a);
  coeff_real.push_back(Algebraic_real_1(1));
  coeff_real.push_back(Algebraic_real_1(-2));

  typedef ::CGAL::internal::Descartes_star< Algebraic_real_1, Rational >
    Isolator_real;
  Isolator_real isolator_real(coeff_real.begin(), coeff_real.end());
  std::cout << "#Roots: " <<  isolator_real.number_of_real_roots() << std::endl;
  std::cout << "1st root: [" << isolator_real.left_bound(0) << "," 
            << isolator_real.right_bound(0) << "]" << std::endl;
  std::cout << "2nd root: [" << isolator_real.left_bound(1) << "," 
            << isolator_real.right_bound(1) << "]" << std::endl;
*/
}


int main(){
  
  CGAL::set_pretty_mode(std::cout);
  CGAL::set_pretty_mode(std::cerr);
  CGAL::set_pretty_mode(std::clog);

#ifdef CGAL_HAS_GMP_ARITHMETIC_KERNEL
  simple_test_descartes_star< CGAL::GMP_arithmetic_kernel >();
//  test_descartes_star<CGAL::GMP_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_LEDA_ARITHMETIC_KERNEL
//  simple_test_descartes_star< CGAL::LEDA_arithmetic_kernel >();
//  test_descartes_star<CGAL::LEDA_arithmetic_kernel>();
#endif
#ifdef CGAL_HAS_CORE_ARITHMETIC_KERNEL
//  simple_test_descartes_star< CGAL::CORE_arithmetic_kernel >();
//  test_descartes_star<CGAL::CORE_arithmetic_kernel>();
#endif

  return EXIT_SUCCESS;
}
// EOF

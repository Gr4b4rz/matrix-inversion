#define BOOST_TEST_MAIN
#include <iostream>
#include <boost/test/unit_test.hpp>
#include "../Matrix.h"

BOOST_AUTO_TEST_CASE(two_on_two_matrix_test)
{
  matrix m = matrix {
    {4, 7},
    {2, 6}
  };

  matrix expected = matrix {
    {0.6,	-0.7},
    {-0.2, 0.4}
  };
  matrix inversed = inverse_matrix(m);
  // TODO add compare method to vector

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_EQUAL(inversed[i][j], expected[i][j]);

}

BOOST_AUTO_TEST_CASE(three_on_three_matrix_test)
{
  std::cout<<"WHAAT";
  matrix m = matrix {
    {0.0, 0.0, 1.0},
    {0.0, 1.0, 5.0},
    {5.0, 6.0, 0.0}
  };

  matrix expected = matrix {
    {-6.0, 3.6, 1.4},
    {5.0, -3.0, -1.0},
    {-1.0, 0.8, -0.2}
  };
  matrix inversed = inverse_matrix(m);

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_EQUAL(inversed[i][j], expected[i][j]);

}

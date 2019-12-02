#define BOOST_TEST_MAIN
#include <iostream>
#include <boost/test/unit_test.hpp>
#include <boost/test/included/unit_test.hpp>
#include <cmath>
#include <fstream>
#include "../Matrix.h"

namespace tt = boost::test_tools;

matrix read_matrix_from_file(std::string path){
  matrix m;
  std::ifstream inputFile(path);
  std::string line;
  std::cout<<path;
  double val;

  while (getline(inputFile, line)){
    std::vector<double> v;
    std::istringstream IS(line);
    while( IS >> val)
       v.push_back(val);

    m.push_back(v);
  }
  inputFile.close();

return m;
}
//TODO change CHECK_CLOSE to BOOST_EQUA

BOOST_AUTO_TEST_CASE(two_on_two_matrix_test)
{
  matrix m = matrix {
    {4.0, 7.0},
    {2.0, 6.0}
  };

  matrix expected = matrix {
    {0.6,	-0.7},
    {-0.2, 0.4}
  };
  matrix inversed = inverse_matrix_iterative(m);

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_CLOSE(inversed[i][j], expected[i][j], EPSILON);

}

BOOST_AUTO_TEST_CASE(three_on_three_matrix_test)
{
  matrix m = matrix {
    {0.0, 0.0, 1.0},
    {0.0, 1.0, 5.0},
    {5.0, 6.0, 0.0}
  };

  matrix expected = matrix {
    {6.0, -1.2, 0.2},
    {-5.0, 1, 0.0},
    {1.0, 0.0, 0.0}
  };
  matrix inversed = inverse_matrix_iterative(m);

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_CLOSE(inversed[i][j], expected[i][j], EPSILON);

}
BOOST_AUTO_TEST_CASE(four_on_four_matrix_test)
{
 matrix m = matrix {
   {1, 1, 1, 11},
   {1, 15, 1, 1},
   {21, 1, 14, 13},
   {1, 1, 1, 5}
 };

 matrix expected {
   {1.4455782312925170066, 0.1326530612244897959, 0.14285714285714285714, -3.5782312925170068025},
   {0.047619047619047619039, 0.071428571428571428571, 0, -0.11904761904761904761},
   {-2.3265306122448979588, -0.20408163265306122448, -0.14285714285714285714, 5.5306122448979591832},
   {0.16666666666666666666, 0, 0, -0.16666666666666666666}
 };

 matrix inversed = inverse_matrix_iterative(m);
 for (int i =0; i < expected.size(); ++i)
     for (int j = 0; j < expected.size(); ++j)
         BOOST_CHECK_CLOSE(inversed[i][j], expected[i][j], EPSILON);
}

BOOST_AUTO_TEST_CASE(twenty_on_twent_matrix_test){
  matrix m = read_matrix_from_file("../inMatrix20.txt");
  matrix expected = read_matrix_from_file("../outMatrix20.txt");
  matrix inversed = inverse_matrix_iterative(m);
  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_CLOSE(inversed[i][j], expected[i][j], EPSILON);

}

BOOST_AUTO_TEST_CASE(transpose_matrix_test)
{
  matrix m = matrix {
    {1, 2, 3, 4},
    {1, 2, 3, 4},
    {1, 2, 3, 4},
    {1, 2, 3, 4}

  };
  matrix expected = matrix{
    {1,1,1,1},
    {2,2,2,2},
    {3,3,3,3},
    {4,4,4,4}
  };

  matrix transposed = transpose_matrix(m);

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_EQUAL(transposed[i][j], expected[i][j]);

}

BOOST_AUTO_TEST_CASE(matrix_trace_test)
{
  matrix m1 = matrix {
    {100, 1, 1, 1},
    {1, 200, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1}
  };

  double trace = calculate_matrix_trace(m1);
  BOOST_CHECK_EQUAL(trace, 302);

}


BOOST_AUTO_TEST_CASE(multiply_matrices_test)
{
  matrix m1 = matrix {
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1}
  };

  matrix m2 = matrix {
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1},
    {1, 1, 1, 1}
  };

  matrix expected = matrix {
    {4, 4, 4, 4},
    {4, 4, 4, 4},
    {4, 4, 4, 4},
    {4, 4, 4, 4}
  };

  matrix multiplied = multiply_matrices(m1, m2);
  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_EQUAL(multiplied[i][j], expected[i][j]);

}

BOOST_AUTO_TEST_CASE(check_R_matrix_when_should_return_false_test)
{
  matrix m = matrix {
    {EPSILON, EPSILON, EPSILON},
    {EPSILON, EPSILON, EPSILON},
    {EPSILON, EPSILON, EPSILON + 0.00001}
  };

  BOOST_CHECK_EQUAL(check_R_matrix(m), false);
}

BOOST_AUTO_TEST_CASE(check_R_matrix_when_should_return_false_true)
{
  matrix m = matrix {
    {EPSILON, EPSILON, EPSILON},
    {EPSILON, EPSILON, EPSILON},
    {EPSILON, EPSILON, EPSILON}
  };

  BOOST_CHECK_EQUAL(check_R_matrix(m), true);
}

BOOST_AUTO_TEST_CASE(multiply_by_scalar_test)
{
  matrix m = matrix {
    {1, 1, 1},
    {1, 1, 1},
    {1, 1, 1}
  };

  matrix expected = matrix {
    {3, 3, 3},
    {3, 3, 3},
    {3, 3, 3}
  };

  matrix multiplied = multiply_by_scalar(m, 3);

  for (int i =0; i < expected.size(); ++i)
      for (int j = 0; j < expected.size(); ++j)
          BOOST_CHECK_EQUAL(multiplied[i][j], expected[i][j]);
}

BOOST_AUTO_TEST_CASE(compare_different_of_inversions_50_X_50)
{
  matrix m = create_random_matrix(50);

  matrix m1 = inverse_matrix_iterative(m);
  matrix m2 = inverse_matrix(m);

  for (int i =0; i < m.size(); ++i)
      for (int j = 0; j < m.size(); ++j)
          BOOST_CHECK_CLOSE(m1[i][j], m2[i][j], 0.00001);

}

BOOST_AUTO_TEST_CASE(compare_different_of_inversions_100_X_100)
{
  matrix m = create_random_matrix(100);

  matrix m1 = inverse_matrix_iterative(m);
  matrix m2 = inverse_matrix(m);

  for (int i =0; i < m.size(); ++i)
      for (int j = 0; j < m.size(); ++j)
          BOOST_CHECK_CLOSE(m1[i][j], m2[i][j], 0.00001);

}

BOOST_AUTO_TEST_CASE(compare_different_of_inversions_200_X_200)
{
  matrix m = create_random_matrix(200);

  matrix m1 = inverse_matrix_iterative(m);
  matrix m2 = inverse_matrix(m);

  for (int i =0; i < m.size(); ++i)
      for (int j = 0; j < m.size(); ++j)
          BOOST_CHECK_CLOSE(m1[i][j], m2[i][j], 0.0001);

}

/******************************************************************************
 * Copyright (c): (2019-2020) GEIRINA
 * All rights reserved.
 * Project: Power system transient simulation
 *
 * - Author:      Peng Wei, peng.wei@geirina.net
 * - Created on:  Jan. 29, 2020
 * - Last Update: Jan. 29, 2020
 *
 * - This library contains general-purpose functions.
 *
*******************************************************************************/

#include "GraphDyn_util.hpp"

namespace transient_analysis {

template <typename T>
void sortrows(vector<vector<T>>& matrix) {
  std::sort(matrix.begin(), matrix.end(),
            [](vector<T> &lhs, vector<T> &rhs)
            {return (lhs[0] != rhs[0]) ? lhs[0] < rhs[0] : lhs[1] < rhs[1];});
}

/* print matrix for debug */
template <typename T>
void print_matrix(const vector<vector<T>>& mtx, const string& info) {
  printf("\n%s\n", info.c_str());
  for (int j = 0; j < mtx.size(); ++j) {
    for (int i = 0; i < mtx[j].size() - 1; ++i)
      printf("%4.6f \t", mtx[j][i]);
    printf("%4.6f\n", mtx[j][mtx[j].size() - 1]);
  }
}

/* print vector for debug */
template <typename T>
void print_vector(const vector<T>& vec, const string& info) {
  printf("\n%s\n", info.c_str());
  for (int i = 0; i < vec.size(); ++i)
    printf("%4.8f \t", 1. * i);
  printf("\n");
  for (int i = 0; i < vec.size(); ++i)
    printf("%4.8f", vec[i]);
  printf("\n\n");
}

/* print array as a complex value for debug */
template <typename T>
void print_array(T* vec, const int size, const string& info) {
  printf("\n%s\n", info.c_str());
  for (int i = 0; i < size; ++i)
    printf("%+2.6f%+2.6fi\n", vec[i], vec[i + size]);
  printf("\n");
}

void printLogo() {
  std::vector<std::string> logo{
      R"(                                               )",
      R"( _____                 _    ______             )",
      R"(|  __ \               | |   |  _  \            )",
      R"(| |  \/_ __ __ _ _ __ | |__ | | | |_   _ _ __  )",
      R"(| | __| '__/ _` | '_ \| '_ \| | | | | | | '_ \ )",
      R"(| |_\ \ | | (_| | |_) | | | | |/ /| |_| | | | |)",
      R"( \____/_|  \__,_| .__/|_| |_|___/  \__, |_| |_|)",
      R"(                | |                 __/ |      )",
      R"(                |_|                |___/       )",
      R"(                                               )"};

  for (auto s : logo)
    std::cout << s << "\n";
}

void print_line() {
  cout << (string(60, '-') + "\n");
}

template void sortrows<real__t>(vector<vector<real__t>>& mtx);
template void print_vector<real__t>(const vector<real__t>& vec, const string& info);
template void print_array<real__t>(real__t* vec, const int size, const string& info);
template void print_matrix<real__t>(const vector<vector<real__t>>& mtx, const string& info);

}  // namespace transient_analysis

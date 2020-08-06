/******************************************************************************
 * Copyright (c) 2016, GraphSQL Inc.
 * All rights reserved.
 * Project: GraphSQL Query Language
 *
 * - This library is for defining struct and helper functions that will be used
 *   in the user-defined functions in "ExprFunctions.hpp". Note that functions
 *   defined in this file cannot be directly called from GraphSQL Query scripts.
 *   Please put such functions into "ExprFunctions.hpp" under the same directory
 *   where this file is located.
 *
 * - Please don't remove necessary codes in this file
 *
 * - A backup of this file can be retrieved at
 *     <graphsql_root_path>/dev_<backup_time>/gdk/gsql/src/QueryUdf/ExprUtil.hpp
 *   after upgrading the system.
 *
 ******************************************************************************/

#ifndef EXPRUTIL_HPP_
#define EXPRUTIL_HPP_

#include <stdio.h>
#include <stdlib.h>
#include <string>

#include <cassert>
#include <mutex>
#include <set>
#include <unordered_map>
#include <unordered_set>
#include <vector>
#include "graphlu.h"
#include "graphlu_util.h"

#include <assert.h>

// for transient simulation
#include <ctime>
#include <fstream>
#include <iostream>
#include <utility>
#include "boost/numeric/odeint.hpp"

#define PI 3.141592653589793
using namespace std;
using namespace boost::numeric::odeint;

typedef boost::numeric::ublas::matrix<real__t> matrix_type;
typedef vector<real__t> state_type;

typedef std::string string;  //XXX DON'T REMOVE

//JIn: Add rpi, cp in loadMatrix and loadAuxiliaryData; Add _rpi, _cp, getRpi, getCp
class CSRMatrix {
 private:
  uint__t _n;
  uint__t _nnz;
  vector<uint__t> _bp;
  //vector<uint__t> _bpp;
  vector<size_t> _p;
  vector<uint__t> _i;
  vector<real__t> _x;
  vector<uint__t> _rp;
  vector<uint__t> _cpi;
  vector<real__t> _rows;
  vector<real__t> _cols;
  vector<uint__t> _rpi;
  vector<uint__t> _cp;
  vector<real__t> _Pn;
  vector<real__t> _Qn;
  vector<real__t> _Va;
  vector<real__t> _Vm;
  vector<uint> _btype;

 public:
  CSRMatrix() : _n(0), _nnz(0) {}

  //JIn: Add rpi, cp in loadMatrix and loadAuxiliaryData
  void loadMatrix(uint__t n, uint__t nnz, size_t *p, uint__t *i, real__t *x, uint__t *rp, uint__t *cpi, real__t *rows, real__t *cols, uint__t *rpi, uint__t *cp) {
    loadMatrix(n, nnz, p, i, x);
    loadAuxiliaryData(n, rp, cpi, rows, cols, rpi, cp);
  }

  void loadMatrix(uint__t n, uint__t nnz, size_t *p, uint__t *i, real__t *x, uint__t *rp, uint__t *cpi, real__t *rows, real__t *cols) {
    loadMatrix(n, nnz, p, i, x);
    loadAuxiliaryData(n, rp, cpi, rows, cols);
  }

  void loadMatrix(uint__t n, uint__t nnz, size_t *p, uint__t *i, real__t *x) {
    assert(n >= 1);
    _n = n;
    _nnz = nnz;
    _p.resize(_n + 1);
    _i.resize(p[_n]);
    _x.resize(p[_n]);

    std::copy(p, p + n + 1, _p.begin());
    std::copy(i, i + p[_n], _i.begin());
    std::copy(x, x + p[_n], _x.begin());
  }

  void loadAuxiliaryData(uint__t n, uint__t *rp, uint__t *cpi, real__t *rows, real__t *cols) {
    assert(n >= 1);
    _n = n;
    _rp.resize(_n);
    _cpi.resize(_n);
    _rows.resize(_n);
    _cols.resize(_n);

    std::copy(rp, rp + n, _rp.begin());
    std::copy(cpi, cpi + n, _cpi.begin());
    std::copy(rows, rows + n, _rows.begin());
    std::copy(cols, cols + n, _cols.begin());
  }

  void loadAuxiliaryData(uint__t n, uint__t *rp, uint__t *cpi, real__t *rows, real__t *cols, uint__t *rpi, uint__t *cp) {
    assert(n >= 1);
    _n = n;
    _rp.resize(_n);
    _cpi.resize(_n);
    _rows.resize(_n);
    _cols.resize(_n);
    _rpi.resize(_n);
    _cp.resize(_n);

    std::copy(rp, rp + n, _rp.begin());
    std::copy(cpi, cpi + n, _cpi.begin());
    std::copy(rows, rows + n, _rows.begin());
    std::copy(cols, cols + n, _cols.begin());
    std::copy(rpi, rpi + n, _rpi.begin());
    std::copy(cp, cp + n, _cp.begin());
  }

  void loadMatrix(uint__t n, uint__t nnz, uint__t *p, uint__t *i, real__t *x, uint__t *btype, real__t *Va, real__t *Vm, real__t *Pn) {
    loadMatrix(n, nnz, p, i, x);
    loadBasecaseStatus(n, btype, Va, Vm, Pn);
  }

  void loadMatrix(uint__t n, uint__t nnz, uint__t *p, uint__t *i, real__t *x) {
    assert(n >= 1);
    _n = n;
    _nnz = nnz;
    _bp.resize(_n + 1);
    _i.resize(p[_n]);
    _x.resize(p[_n]);

    std::copy(p, p + n + 1, _bp.begin());
    std::copy(i, i + p[_n], _i.begin());
    std::copy(x, x + p[_n], _x.begin());
  }

  void loadBasecaseStatus(uint__t n, uint__t *btype, real__t *Va, real__t *Vm, real__t *Pn) {
    assert(n >= 1);
    _n = n;
    _btype.resize(_n);
    _Va.resize(_n);
    _Vm.resize(_n);
    _Pn.resize(_n);

    std::copy(btype, btype + n, _btype.begin());
    std::copy(Va, Va + n, _Va.begin());
    std::copy(Vm, Vm + n, _Vm.begin());
    std::copy(Pn, Pn + n, _Pn.begin());
  }

  bool isValid() const { return (_n >= 1) && (_p.size() == _n + 1) && (_p[_n] >= _n) && (_i.size() == _p[_n]) && (_x.size() == _p[_n]); }

  bool hasAuxiliaryData() const { return !_rp.empty(); }

  const uint__t *getBtype() const { return &_btype[0]; }
  const real__t *getVa() const { return &_Va[0]; }
  const real__t *getVm() const { return &_Vm[0]; }
  const real__t *getPn() const { return &_Pn[0]; }
  const uint__t *getBp() const { return &_bp[0]; }

  const size_t *getP() const { return &_p[0]; }
  const uint__t *getI() const { return &_i[0]; }
  const real__t *getX() const { return &_x[0]; }
  const uint__t *getRp() const { return &_rp[0]; }
  const uint__t *getCpi() const { return &_cpi[0]; }
  const real__t *getRows() const { return &_rows[0]; }
  const real__t *getCols() const { return &_cols[0]; }
  const uint__t *getRpi() const { return &_rpi[0]; }
  const uint__t *getCp() const { return &_cp[0]; }
  const uint__t getn() const { return _n; }
  const uint__t getnnz() const { return _nnz; }
};

class SingletonMatrixInterface {
 private:
  static CSRMatrix *getOrDeleteMatrixeG(bool isDelete) {
    static CSRMatrix *eG = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (eG != NULL) {
        delete eG;
        eG = NULL;
      }
      assert(!eG);
    } else {
      if (eG == NULL) {
        eG = new CSRMatrix();
      }
      assert(eG);
    }
    mtx.unlock();
    return eG;
  }

  static CSRMatrix *getOrDeleteMatrixeB(bool isDelete) {
    static CSRMatrix *eB = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (eB != NULL) {
        delete eB;
        eB = NULL;
      }
      assert(!eB);
    } else {
      if (eB == NULL) {
        eB = new CSRMatrix();
      }
      assert(eB);
    }
    mtx.unlock();
    return eB;
  }

  static CSRMatrix *getOrDeleteMatrixLineQ1(bool isDelete) {
    static CSRMatrix *LineQ1 = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (LineQ1 != NULL) {
        delete LineQ1;
        LineQ1 = NULL;
      }
      assert(!LineQ1);
    } else {
      if (LineQ1 == NULL) {
        LineQ1 = new CSRMatrix();
      }
      assert(LineQ1);
    }
    mtx.unlock();
    return LineQ1;
  }

  static CSRMatrix *getOrDeleteMatrixBp(bool isDelete) {
    static CSRMatrix *Bp = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Bp != NULL) {
        delete Bp;
        Bp = NULL;
      }
      assert(!Bp);
    } else {
      if (Bp == NULL) {
        Bp = new CSRMatrix();
      }
      assert(Bp);
    }
    mtx.unlock();
    return Bp;
  }

  static CSRMatrix *getOrDeleteMatrixBpp(bool isDelete) {
    static CSRMatrix *Bpp = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Bpp != NULL) {
        delete Bpp;
        Bpp = NULL;
      }
      assert(!Bpp);
    } else {
      if (Bpp == NULL) {
        Bpp = new CSRMatrix();
      }
      assert(Bpp);
    }
    mtx.unlock();
    return Bpp;
  }

  static CSRMatrix *getOrDeleteMatrixLp(bool isDelete) {
    static CSRMatrix *Lp = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Lp != NULL) {
        delete Lp;
        Lp = NULL;
      }
      assert(!Lp);
    } else {
      if (Lp == NULL) {
        Lp = new CSRMatrix();
      }
      assert(Lp);
    }
    mtx.unlock();
    return Lp;
  }

  static CSRMatrix *getOrDeleteMatrixUp(bool isDelete) {
    static CSRMatrix *Up = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Up != NULL) {
        delete Up;
        Up = NULL;
      }
      assert(!Up);
    } else {
      if (Up == NULL) {
        Up = new CSRMatrix();
      }
      assert(Up);
    }
    mtx.unlock();
    return Up;
  }

  static CSRMatrix *getOrDeleteMatrixLpp(bool isDelete) {
    static CSRMatrix *Lpp = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Lpp != NULL) {
        delete Lpp;
        Lpp = NULL;
      }
      assert(!Lpp);
    } else {
      if (Lpp == NULL) {
        Lpp = new CSRMatrix();
      }
      assert(Lpp);
    }
    mtx.unlock();
    return Lpp;
  }

  static CSRMatrix *getOrDeleteMatrixUpp(bool isDelete) {
    static CSRMatrix *Upp = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (Upp != NULL) {
        delete Upp;
        Upp = NULL;
      }
      assert(!Upp);
    } else {
      if (Upp == NULL) {
        Upp = new CSRMatrix();
      }
      assert(Upp);
    }
    mtx.unlock();
    return Upp;
  }

  // static CSRMatrix *getOrDeleteMatrixGainP(bool isDelete) {
  // static CSRMatrix *GainP = NULL;
  // static std::mutex mtx;
  // mtx.lock();
  // if (isDelete) {
  // if (GainP != NULL) {
  // delete GainP;
  // GainP = NULL;
  // }
  // assert(!GainP);
  // } else {
  // if (GainP == NULL) {
  // GainP = new CSRMatrix();
  // }
  // assert(GainP);
  // }
  // mtx.unlock();
  // return GainP;
  // }

  // static CSRMatrix *getOrDeleteMatrixGainQ(bool isDelete) {
  // static CSRMatrix *GainQ = NULL;
  // static std::mutex mtx;
  // mtx.lock();
  // if (isDelete) {
  // if (GainQ != NULL) {
  // delete GainQ;
  // GainQ = NULL;
  // }
  // assert(!GainQ);
  // } else {
  // if (GainQ == NULL) {
  // GainQ = new CSRMatrix();
  // }
  // assert(GainQ);
  // }
  // mtx.unlock();
  // return GainQ;
  // }

  static CSRMatrix *getOrDeleteMatrixGainP_L(bool isDelete) {
    static CSRMatrix *GainP_L = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (GainP_L != NULL) {
        delete GainP_L;
        GainP_L = NULL;
      }
      assert(!GainP_L);
    } else {
      if (GainP_L == NULL) {
        GainP_L = new CSRMatrix();
      }
      assert(GainP_L);
    }
    mtx.unlock();
    return GainP_L;
  }

  static CSRMatrix *getOrDeleteMatrixGainP_U(bool isDelete) {
    static CSRMatrix *GainP_U = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (GainP_U != NULL) {
        delete GainP_U;
        GainP_U = NULL;
      }
      assert(!GainP_U);
    } else {
      if (GainP_U == NULL) {
        GainP_U = new CSRMatrix();
      }
      assert(GainP_U);
    }
    mtx.unlock();
    return GainP_U;
  }

  static CSRMatrix *getOrDeleteMatrixGainQ_L(bool isDelete) {
    static CSRMatrix *GainQ_L = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (GainQ_L != NULL) {
        delete GainQ_L;
        GainQ_L = NULL;
      }
      assert(!GainQ_L);
    } else {
      if (GainQ_L == NULL) {
        GainQ_L = new CSRMatrix();
      }
      assert(GainQ_L);
    }
    mtx.unlock();
    return GainQ_L;
  }

  static CSRMatrix *getOrDeleteMatrixGainQ_U(bool isDelete) {
    static CSRMatrix *GainQ_U = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isDelete) {
      if (GainQ_U != NULL) {
        delete GainQ_U;
        GainQ_U = NULL;
      }
      assert(!GainQ_U);
    } else {
      if (GainQ_U == NULL) {
        GainQ_U = new CSRMatrix();
      }
      assert(GainQ_U);
    }
    mtx.unlock();
    return GainQ_U;
  }

 public:
  static CSRMatrix *getMatrixeG() { return getOrDeleteMatrixeG(false); }
  static void deleteMatrixeG() { getOrDeleteMatrixeG(true); }

  static CSRMatrix *getMatrixeB() { return getOrDeleteMatrixeB(false); }
  static void deleteMatrixeB() { getOrDeleteMatrixeB(true); }

  static CSRMatrix *getMatrixLineQ1() { return getOrDeleteMatrixLineQ1(false); }
  static void deleteMatrixLineQ1() { getOrDeleteMatrixLineQ1(true); }

  static CSRMatrix *getMatrixBp() { return getOrDeleteMatrixBp(false); }
  static void deleteMatrixBp() { getOrDeleteMatrixBp(true); }

  static CSRMatrix *getMatrixBpp() { return getOrDeleteMatrixBpp(false); }
  static void deleteMatrixBpp() { getOrDeleteMatrixBpp(true); }

  static CSRMatrix *getMatrixLp() { return getOrDeleteMatrixLp(false); }
  static void deleteMatrixLp() { getOrDeleteMatrixLp(true); }

  static CSRMatrix *getMatrixUp() { return getOrDeleteMatrixUp(false); }
  static void deleteMatrixUp() { getOrDeleteMatrixUp(true); }

  static CSRMatrix *getMatrixLpp() { return getOrDeleteMatrixLpp(false); }
  static void deleteMatrixLpp() { getOrDeleteMatrixLpp(true); }

  static CSRMatrix *getMatrixUpp() { return getOrDeleteMatrixUpp(false); }
  static void deleteMatrixUpp() { getOrDeleteMatrixUpp(true); }

  // static CSRMatrix *getMatrixGainP() {return getOrDeleteMatrixGainP(false); }
  // static void deleteMatrixGainP() { getOrDeleteMatrixGainP(true); }

  // static CSRMatrix *getMatrixGainQ() {return  getOrDeleteMatrixGainQ(false); }
  // static void deleteMatrixGainQ() { getOrDeleteMatrixGainQ(true); }

  static CSRMatrix *getMatrixGainP_L() { return getOrDeleteMatrixGainP_L(false); }
  static void deleteMatrixGainP_L() { getOrDeleteMatrixGainP_L(true); }

  static CSRMatrix *getMatrixGainP_U() { return getOrDeleteMatrixGainP_U(false); }
  static void deleteMatrixGainP_U() { getOrDeleteMatrixGainP_U(true); }

  static CSRMatrix *getMatrixGainQ_L() { return getOrDeleteMatrixGainQ_L(false); }
  static void deleteMatrixGainQ_L() { getOrDeleteMatrixGainQ_L(true); }

  static CSRMatrix *getMatrixGainQ_U() { return getOrDeleteMatrixGainQ_U(false); }
  static void deleteMatrixGainQ_U() { getOrDeleteMatrixGainQ_U(true); }
};

// **********************************************************************************
////Jingjin: 2018/08/18
////Split Symbolic and Numurical factorization, mainly include 3 functions: symbolic_factorize, numerical_factorize, and Refactorize_preprocess;
////Refactorize can be performed by storing the permuation and scaling information first.
////symbolic factorize can be skipped during refactorizatio
////symbolic factorize can be skipped during refactorization.
//// *************************************************************************************
class Graphlu_Util {
 public:
  static int symbolic_factorize(SGraphLU *matrix) {
    if (matrix == NULL) {
      return -1;
    }
    int ret = GraphLU_Analyze(matrix);
    // std::cout << "analysis time: " << matrix->stat[0] << "\n"
    //           << std::endl;
    //  std::cout << "Symbolic debug: " << ret << "\n" <<  std::endl;
    return ret;
  }

  static int numerical_factorize(SGraphLU *matrix) {
    int error;
    if (matrix == NULL) {
      return -1;
    }

    if (matrix->cfgi[9] == 1) {
      error = GraphLU_Factorize(matrix);
      // printf("Factorization time: %.8g\n", matrix->stat[1]);
      if (error < 0)  //there is an error, print the code
        std::cout << "GraphLU_Factorize error code:" << error << std::endl;
    } else {
      int ret = GraphLU_CreateScheduler(matrix);
      printf("time of creating scheduler: %.8g\n", matrix->stat[4]);
      printf("suggestion: %s.\n", ret == 0 ? "parallel" : "sequential");

      // This function creates threads for parallel computation. The second argument (thread)
      // specifies the number of threads, including the main thread. The last argument (check)
      // specifies whether to check the number of threads or not.
      if (!ret) {  // parallel factorization
        GraphLU_CreateThreads(matrix, 2, TRUE);
        printf("total cores: %d, threads created: %d\n", (int)(matrix->stat[9]), (int)(matrix->cfgi[5]));

        // This function binds threads to cores (unbind = FALSE) or unbinds threads from cores (unbind = TRUE).
        GraphLU_BindThreads(matrix, FALSE);
        // Numerical LU factorization with partial pivoting, parallel
        error = GraphLU_Factorize_MT(matrix);
        printf("factorization time: %.8g\n", matrix->stat[1]);
        if (error < 0)  //there is an error, print the code
          std::cout << "GraphLU_Factorize error code:" << error << std::endl;
      }
      // [Ting][end] comment out if to try sequential factorization
      else {
        // Sequential factorization
        error = GraphLU_Factorize(matrix);
        // printf("Factorization time: %.8g\n", matrix->stat[1]);
        if (error < 0)  //there is an error, print the code
          std::cout << "GraphLU_Factorize error code:" << error << std::endl;
      }
      //printf("Numerical debug: %d\n", error);
    }
    return error;
  }

  static int Refactorize_preprocess(uint__t n, SGraphLU *matrix, CSRMatrix *refMatrix) {
    if (refMatrix == NULL) {
      return -1;
    }
    assert(n >= 1);
    matrix->flag[1] = TRUE;
    std::copy(refMatrix->getRp(), refMatrix->getRp() + n, matrix->row_perm);
    std::copy(refMatrix->getRpi(), refMatrix->getRpi() + n, matrix->row_perm_inv);
    std::copy(refMatrix->getCp(), refMatrix->getCp() + n, matrix->col_perm);
    std::copy(refMatrix->getCpi(), refMatrix->getCpi() + n, matrix->col_perm_inv);
    std::copy(refMatrix->getRows(), refMatrix->getRows() + n, matrix->row_scale);
    std::copy(refMatrix->getCols(), refMatrix->getCols() + n, matrix->col_scale_perm);

    permute_and_scale(matrix);
    return 0;
  }
};

class MatrixGraph {
 public:
  MatrixGraph(int n, uint__t *ai, uint__t *ap)
      : _n(n) {
    _node2edges.resize(_n + _n);
    for (int i = 0; i < _n; ++i) {
      int end = ap[i + 1];
      for (int j = ap[i]; j < end; ++j) {
        if (ai[j] != i) {
          _node2edges[i].insert(ai[j]);
          _node2edges[ai[j]].insert(i);
        }
      }
    }
  }

  void findSingleDegreeNodes(vector<int> &singleDegreeNodes) const {
    singleDegreeNodes.clear();
    for (int i = 0; i < _n; ++i) {
      if (isSingleDegreeNode(i)) {
        singleDegreeNodes.push_back(i);
      }
    }
  }

  void findBridges(vector<pair<int, int>> &edges) const {
    edges.clear();
    int t = 0;
    vector<bool> visited(_n, false);
    vector<int> disc(_n);
    vector<int> low(_n);
    vector<int> parent(_n, -1);
    for (int i = 0; i < _n; ++i) {
      if (!visited[i]) {
        findBridgesInt(i, t, visited, disc, low, parent, edges);
      }
    }
  }

  void findBridgesWithoutSingleDegreeNodes(vector<pair<int, int>> &edges) const {
    edges.clear();
    vector<pair<int, int>> bridges;
    findBridges(bridges);
    for (int i = 0; i < bridges.size(); ++i) {
      if (isSingleDegreeNode(bridges[i].first) ||
          isSingleDegreeNode(bridges[i].second)) {
        continue;
      }
      edges.push_back(bridges[i]);
    }
  }

 private:
  void findBridgesInt(int u,
                      int &t,
                      vector<bool> &visited,
                      vector<int> &disc,
                      vector<int> &low,
                      vector<int> &parent,
                      vector<pair<int, int>> &edges) const {
    visited[u] = true;
    ++t;
    disc[u] = t;
    low[u] = t;
    for (set<int>::const_iterator iter = _node2edges[u].begin();
         iter != _node2edges[u].end(); ++iter) {
      int v = *iter;
      if (!visited[v]) {
        parent[v] = u;
        findBridgesInt(v, t, visited, disc, low, parent, edges);
        low[u] = min(low[u], low[v]);
        if (low[v] > disc[u]) {
          edges.push_back(make_pair(u, v));
        }
      } else if (v != parent[u]) {
        low[u] = min(low[u], disc[v]);
      }
    }
  }
  bool isSingleDegreeNode(int i) const {
    return _node2edges[i].size() == 1;
  }

 private:
  int _n;
  vector<set<int>> _node2edges;
};

// **********************************************************************************
//Jingjin: 2019/01/07
//static unordered map to store information, which can be used in other functions.
//Also implement some basic operations for unordered map: get, delete, insert, find,check, etc.
// *************************************************************************************
typedef unordered_map<string, vector<int64_t>> EdgeToIslandNodes;

class SingletonEdgeToIslandNodesMap {
 private:
  static EdgeToIslandNodes *edgeToIslandNodesMapOpsInt(bool isTableDelete,
                                                       bool isEntryAdd, bool isEntryDelete,
                                                       bool isEntryCheck, bool isEntryFind,
                                                       EdgeToIslandNodes *&table,
                                                       const string &key, const vector<int64_t> &value,
                                                       bool &success, EdgeToIslandNodes::const_iterator &iter) {
    static EdgeToIslandNodes *edgeToIslandNodesMap = NULL;
    static std::mutex mtx;
    mtx.lock();
    if (isTableDelete) {
      if (edgeToIslandNodesMap != NULL) {
        delete edgeToIslandNodesMap;
        edgeToIslandNodesMap = NULL;
      }
      assert(!edgeToIslandNodesMap);
    } else {
      if (edgeToIslandNodesMap == NULL) {
        edgeToIslandNodesMap = new EdgeToIslandNodes();
      }
      assert(edgeToIslandNodesMap);
    }
    table = edgeToIslandNodesMap;
    if (isEntryAdd) {
      success = edgeToIslandNodesMap->insert(EdgeToIslandNodes::value_type(key, value)).second;
      //   std::cout << "Insert debug: " << success << std::endl;

    } else if (isEntryDelete) {
      success = (edgeToIslandNodesMap->erase(key) == 1);
    } else if (isEntryCheck) {
      success = (edgeToIslandNodesMap->count(key) == 1);
    } else if (isEntryFind) {
      iter = edgeToIslandNodesMap->find(key);
    }
    mtx.unlock();
    /*
     //Test
     if (isEntryFind) {
     iter = edgeToIslandNodesMap->find(key);
     }
     */
    return edgeToIslandNodesMap;
  }

  static EdgeToIslandNodes *getOrDeleteEdgeToIslandNodesMap(bool isDelete) {
    bool isTableDelete = isDelete;
    bool isEntryAdd = false;
    bool isEntryDelete = false;
    bool isEntryCheck = false;
    bool isEntryFind = false;
    EdgeToIslandNodes *table = NULL;
    string key = "";
    vector<int64_t> value;
    bool success = false;
    EdgeToIslandNodes::const_iterator iter;
    edgeToIslandNodesMapOpsInt(isTableDelete,
                               isEntryAdd, isEntryDelete, isEntryCheck, isEntryFind,
                               table, key, value, success, iter);
    return table;
  }

 public:
  static bool insert(const string &key, const vector<int64_t> &value) {
    bool isTableDelete = false;
    bool isEntryAdd = true;
    bool isEntryDelete = false;
    bool isEntryCheck = false;
    bool isEntryFind = false;
    EdgeToIslandNodes *table = NULL;
    bool success = false;
    EdgeToIslandNodes::const_iterator iter;
    edgeToIslandNodesMapOpsInt(isTableDelete,
                               isEntryAdd, isEntryDelete, isEntryCheck, isEntryFind,
                               table, key, value, success, iter);
    return success;
  }

  static bool erase(const string &key) {
    bool isTableDelete = false;
    bool isEntryAdd = false;
    bool isEntryDelete = true;
    bool isEntryCheck = false;
    bool isEntryFind = false;
    EdgeToIslandNodes *table = NULL;
    vector<int64_t> value;
    bool success = false;
    EdgeToIslandNodes::const_iterator iter;
    edgeToIslandNodesMapOpsInt(isTableDelete, isEntryAdd,
                               isEntryDelete, isEntryCheck, isEntryFind,
                               table, key, value, success, iter);
    return success;
  }

  static bool contains(const string &key) {
    bool isTableDelete = false;
    bool isEntryAdd = false;
    bool isEntryDelete = false;
    bool isEntryCheck = true;
    bool isEntryFind = false;
    EdgeToIslandNodes *table = NULL;
    vector<int64_t> value;
    bool success = false;
    EdgeToIslandNodes::const_iterator iter;
    edgeToIslandNodesMapOpsInt(isTableDelete,
                               isEntryAdd, isEntryDelete, isEntryCheck, isEntryFind,
                               table, key, value, success, iter);
    return success;
  }

  static EdgeToIslandNodes::const_iterator find(const string &key) {
    bool isTableDelete = false;
    bool isEntryAdd = false;
    bool isEntryDelete = false;
    bool isEntryCheck = false;
    bool isEntryFind = true;
    EdgeToIslandNodes *table = NULL;
    vector<int64_t> value;
    bool success = false;
    EdgeToIslandNodes::const_iterator iter;
    edgeToIslandNodesMapOpsInt(isTableDelete,
                               isEntryAdd, isEntryDelete, isEntryCheck, isEntryFind,
                               table, key, value, success, iter);
    return iter;
  }

  static EdgeToIslandNodes *getEdgeToIslandNodesMap() { return getOrDeleteEdgeToIslandNodesMap(false); }
  static void deleteEdgeToIslandNodesMap() { getOrDeleteEdgeToIslandNodesMap(true); }
};

#endif /* EXPRUTIL_HPP_ */

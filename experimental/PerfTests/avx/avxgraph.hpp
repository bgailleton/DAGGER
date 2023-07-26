#pragma once

#include "avxcon.hpp"

template <class Connector_t> class graphavx {
public:
  Connector_t *connector;
  std::vector<uint32_t> stack;
  std::vector<int> argstack;
  std::vector<int> argrec;
  std::vector<std::array<int, 8>> stackarray;
  std::vector<int> leveler;

  graphavx() { ; };

  graphavx(Connector_t &con) {
    this->connector = &con;
    this->stack = std::vector<uint32_t>(this->connector->nxy);
  }

  // This is my implementation of Braun and willett 2014
  // Slightly modified:
  // - First it is based on the fastscapelib version, which bypass the need of
  // the delta vectors and all by simply applying the recursion directly feeding
  // the stack
  // - Secondly, recursion is not the best practice for most languages so I am
  // using a stack data structure instead Results are identical and performances
  // similar (marginally better, but that is linked to c++ not being a heavy
  // recursion friendly language)
  void toposort_SFD_classic() {
    // The stack container helper
    std::stack<int, std::vector<int>> stackhelper;
    // std::vector<bool> isdone(this->nnodes,false);
    // going through all the nodes
    std::array<int, 8> arr;
    int istack = 0;
    for (int i = 0; i < this->connector->nxy; ++i) {
      // if they are base level I include them in the stack
      if (this->connector->Sreceivers(i) == i) {
        stackhelper.emplace(i);
        // ++istack;
      }

      // While I still have stuff in the stack helper
      while (stackhelper.empty() == false) {
        // I get the next node and pop it from the stack helper
        int nextnode = stackhelper.top();
        stackhelper.pop();
        this->stack[istack] = static_cast<uint32_t>(nextnode);
        ++istack;
        // std::cout << istack << "/" << this->connector->nxy << std::endl;
        int nn = this->connector->SDonors(nextnode, arr);
        // as well as all its donors which will be processed next
        for (int j = 0; j < nn; ++j) {
          stackhelper.emplace(arr[j]);
        }
      }
    }
  }

  void toposort_SFD() {
    this->stack.clear();
    this->stack.reserve(this->connector->nxy);
    this->leveler.clear();
    this->leveler.reserve(this->connector->nxy);
    this->argstack.clear();
    this->argstack.resize(this->connector->nxy);
    this->argrec.clear();
    this->argrec.resize(this->connector->nxy);

    std::queue<int> QA;
    int tlevel = 0;
    std::array<int, 8> arr;

    for (int i = 0; i < this->connector->nxy; ++i) {
      if (this->connector->Sreceivers(i) == i) {
        this->stack.emplace_back(static_cast<uint32_t>(i));
        this->argstack[i] = this->stack.size() - 1;
        ++tlevel;
        int nn = this->connector->SDonors(i, arr);
        for (int j = 0; j < nn; ++j)
          QA.emplace(arr[j]);
      }
    }

    QA.emplace(-1);
    this->leveler.emplace_back(tlevel);
    tlevel = 0;

    int next = 0;
    while (QA.empty() == false) {
      next = QA.front();
      QA.pop();
      if (next == -1) {
        this->leveler.emplace_back(tlevel);
        tlevel = 0;
        if (QA.empty() == false)
          QA.emplace(-1);
        continue;
      }
      this->stack.emplace_back(static_cast<uint32_t>(next));
      this->argstack[next] = this->stack.size() - 1;
      ++tlevel;
      int nn = this->connector->SDonors(next, arr);
      for (int j = 0; j < nn; ++j)
        QA.emplace(arr[j]);
    }

    for (int i = 0; i < this->connector->nxy; ++i) {
      this->argrec[i] =
          this->argstack[this->connector->Sreceivers(this->stack[i])];
    }

    std::cout << "DEBUGLOG::Stacksize::" << this->stack.size() << " vs "
              << this->connector->nxy << std::endl;
  }

  void toposort_SFD_beta() {
    this->stack.clear();
    this->stack.reserve(this->connector->nxy);

    bool opti_space = false;
    // std::cout << "B" << std::endl;

    std::array<std::stack<int, std::vector<int>>, 8> nextnodes;
    // std::cout << "C" << std::endl;

    std::vector<int> bezlvl;
    bezlvl.reserve(this->connector->nx * 10);
    // std::cout << "D" << std::endl;

    for (int i = 0; i < this->connector->nxy; ++i) {
      if (this->connector->Sreceivers(i) == i) {
        bezlvl.emplace_back(i);
      }
    }
    // std::cout << "E" << std::endl;

    for (int i = 0; i < 8; ++i)
      nextnodes[i].emplace(bezlvl[i]);

    // std::cout << "A" << std::endl;
    int bslr = 8;
    int nempty = 1;
    std::array<int, 8> neighbours;
    int nn;

    while (nempty != 8) {
      nempty = 0;
      bool reverse = false;
      for (int i = 0; i < 8; ++i) {
        if (nextnodes[i].empty()) {
          ++nempty;
          continue;
        }

        int node = nextnodes[i].top();
        nextnodes[i].pop();

        this->stack.emplace_back(node);
        nn = this->connector->SDonors(node, neighbours);
        for (int j = 0; j < nn; ++j)
          nextnodes[i].emplace(neighbours[j]);

        if (nextnodes[i].empty()) {
          if (bslr < static_cast<int>(bezlvl.size())) {
            nextnodes[i].emplace(bezlvl[bslr]);
            ++bslr;
          } else if (opti_space) {
            for (int ii = 0; ii < 8; ++ii) {
              if (nextnodes[ii].size() > 1) {
                // reverse = true;
                std::vector<int> &cont = Container(nextnodes[ii]);
                nextnodes[i].emplace(cont.back());
                cont.pop_back();
                break;
              }
            }
          }
        }
      }

      // padding
      for (int i = 0; i < nempty; ++i) {
        this->stack.emplace_back(0);
      }

      // if(reverse)
      // {
      //   for(int i =0; i<4; ++i)
      //     std::swap(nextnodes[i],nextnodes[7-i]);
      // }
    }
  }

  void toposort_SFD_beta_1() {
    this->stack.clear();
    this->stack.reserve(this->connector->nxy);
    // this->leveler.clear();
    // this->leveler.reserve(this->connector->nxy);
    this->argstack.clear();
    this->argstack.resize(this->connector->nxy);
    this->argrec.clear();
    this->argrec.resize(this->connector->nxy);

    std::queue<int> QA;
    int tlevel = 0;
    std::array<int, 8> arr;

    for (int i = 0; i < this->connector->nxy; ++i) {
      if (this->connector->Sreceivers(i) == i) {
        this->stack.emplace_back(static_cast<uint32_t>(i));

        // this->argstack[i] = this->stack.size()-1;
        ++tlevel;
        int nn = this->connector->SDonors(i, arr);
        for (int j = 0; j < nn; ++j)
          QA.emplace(arr[j]);
      }
    }

    // QA.emplace(-1);
    std::cout << this->stack.size() - tlevel << std::endl;
    reversedRadixSort(this->stack, this->stack.size() - tlevel,
                      this->stack.size() - 1);
    // radixSort(this->stack, this->stack.size() - tlevel,this->stack.size()-1);
    // this->leveler.emplace_back(tlevel);
    tlevel = 0;

    int next = 0;
    while (QA.empty() == false) {
      next = QA.front();
      QA.pop();
      if (next == -1) {
        // this->leveler.emplace_back(tlevel);
        reversedRadixSort(this->stack, this->stack.size() - tlevel,
                          this->stack.size() - 1);
        // radixSort(this->stack, this->stack.size() -
        // tlevel,this->stack.size()-1);
        tlevel = 0;

        // this->leveler.emplace_back(tlevel);
        if (QA.empty() == false)
          QA.emplace(-1);
        continue;
      }
      this->stack.emplace_back(static_cast<uint32_t>(next));
      // this->argstack[next] = this->stack.size()-1;
      ++tlevel;
      int nn = this->connector->SDonors(next, arr);
      for (int j = 0; j < nn; ++j)
        QA.emplace(arr[j]);
    }

    // for(int i=0;i<this->connector->nxy;++i)
    // {
    // 	this->argrec[i] =
    // this->argstack[this->connector->Sreceivers(this->stack[i])];
    // }

    std::cout << "DEBUGLOG::Stacksize::" << this->stack.size() << " vs "
              << this->connector->nxy << std::endl;
  }

  std::vector<float> _compute_drainage_area() {
    std::vector<float> tA(this->connector->nxy,
                          this->connector->get_cell_area(0));
    float *ptrA = tA.data();
    int nremaining = this->connector->nxy % 8;
    __m256 DA = _mm256_set1_ps(this->connector->get_cell_area(0));
    __m256 m1 = _mm256_set1_ps(-1);
    __m256 p1 = _mm256_set1_ps(1);
    __m256i nodes;
    __m256i recs;
    __m256 ttAs;
    __m256 ttArecs;
    __m256i cmpResult;
    __m256 cmpResultf;
    __m256 tadd;

    for (int i = this->connector->nxy - 1; i > 7; i -= 8) {
      // // int node = this->stack[i];
      // // __m256i nodes = _mm256_loadu_si256(reinterpret_cast<const
      // __m256i*>(&this->stack[i+7])); nodes =
      // _mm256_setr_epi32(this->stack[i-0],this->stack[i-1],this->stack[i-2],this->stack[i-3],this->stack[i-4],this->stack[i-5],this->stack[i-6],this->stack[i-7]);
      // recs =
      // _mm256_setr_epi32(this->connector->Sreceivers(this->stack[i-0]),this->connector->Sreceivers(this->stack[i-1]),this->connector->Sreceivers(this->stack[i-2]),this->connector->Sreceivers(this->stack[i-3]),this->connector->Sreceivers(this->stack[i-4]),this->connector->Sreceivers(this->stack[i-5]),this->connector->Sreceivers(this->stack[i-6]),this->connector->Sreceivers(this->stack[i-7]));
      // // printM256i(nodes);
      // // printM256i(recs);
      // ttAs = _mm256_i32gather_ps(ptrA, nodes, sizeof(float));
      // ttArecs = _mm256_set1_ps(0.);
      // ttAs = _mm256_add_ps(ttAs, DA);

      // // if(i > this->connector->nx * 10)
      // // {
      // // 	printM256i(nodes);
      // // 	printM256i(recs);
      // // }
      // cmpResult = _mm256_cmpeq_epi32(nodes, recs);
      // // if(i > this->connector->nx * 10) printM256i(cmpResult);

      // cmpResultf = _mm256_cvtepi32_ps(cmpResult);
      // cmpResultf = _mm256_add_ps(cmpResultf,p1);

      // tadd = _mm256_mul_ps(ttAs,cmpResultf);
      // // if(i > this->connector->nx * 10) printM256(cmpResultf);
      // // if(i > this->connector->nx * 10) printM256(tadd);
      // // if(i > this->connector->nx * 10) printM256(ttArecs);
      // ttArecs = _mm256_add_ps(tadd,ttArecs);
      // // if(i > this->connector->nx * 10) printM256(ttArecs);
      std::array<uint32_t, 8> ino{this->stack[i - 0], this->stack[i - 1],
                                  this->stack[i - 2], this->stack[i - 3],
                                  this->stack[i - 4], this->stack[i - 5],
                                  this->stack[i - 6], this->stack[i - 7]};
      std::array<int, 8> recs{this->connector->Sreceivers(ino[0]),
                              this->connector->Sreceivers(ino[1]),
                              this->connector->Sreceivers(ino[2]),
                              this->connector->Sreceivers(ino[3]),
                              this->connector->Sreceivers(ino[4]),
                              this->connector->Sreceivers(ino[5]),
                              this->connector->Sreceivers(ino[6]),
                              this->connector->Sreceivers(ino[7])};

      ttAs = _mm256_setr_ps(tA[ino[0]], tA[ino[1]], tA[ino[2]], tA[ino[3]],
                            tA[ino[4]], tA[ino[5]], tA[ino[6]], tA[ino[7]]);
      ttArecs =
          _mm256_setr_ps(tA[recs[0]], tA[recs[1]], tA[recs[2]], tA[recs[3]],
                         tA[recs[4]], tA[recs[5]], tA[recs[6]], tA[recs[7]]);
      ttArecs = _mm256_add_ps(ttAs, ttArecs);

      // if(i > this->connector->nx * 10) std::cout << std::endl;

      alignas(32) float rtArec[8];
      // alignas(32) int rnodes[8], rrecs[8];

      // _mm256_store_ps(rtA, ttAs);
      _mm256_store_ps(rtArec, ttArecs);
      // _mm256_storeu_si256((__m256i*)rnodes,nodes);
      // _mm256_storeu_si256((__m256i*)rrecs,recs);
      // std::cout << "J" << std::endl;

      // printM256i(nodes);
      // printM256i(recs);

      // printM256(ttAs);
      // printM256(ttArecs);

      tA[recs[0]] = rtArec[0];
      tA[recs[1]] = rtArec[1];
      tA[recs[2]] = rtArec[2];
      tA[recs[3]] = rtArec[3];
      tA[recs[4]] = rtArec[4];
      tA[recs[5]] = rtArec[5];
      tA[recs[6]] = rtArec[6];
      tA[recs[7]] = rtArec[7];
      // tA[rrecs[0]] += rtArec[0];
      // tA[rrecs[1]] += rtArec[1];
      // tA[rrecs[2]] += rtArec[2];
      // tA[rrecs[3]] += rtArec[3];
      // tA[rrecs[4]] += rtArec[4];
      // tA[rrecs[5]] += rtArec[5];
      // tA[rrecs[6]] += rtArec[6];
      // tA[rrecs[7]] += rtArec[7];

      // tA[node] += this->connector->get_cell_area(node);
      // int rec = this->connector->Sreceivers(node);
      // if(rec != node)
      // {
      // 	tA[rec] += tA[node];
      // }
    }
    return tA;
  }

  std::vector<float> _compute_drainage_area_beta() {

    std::vector<float> tDA(this->stack.size(),
                           this->connector->get_cell_area(0));

    // Perform vectorized addition A[i + 8] += A[i] using AVX2
    for (int i = 0; i < this->stack.size() - 8; i += 8) {
      __m256 vecA =
          _mm256_loadu_ps(&tDA[i]); // Load 8 elements from A[i] to A[i + 7]
      __m256 vecB = _mm256_loadu_ps(
          &tDA[i + 8]); // Load 8 elements from A[i + 8] to A[i + 15]
      __m256 result = _mm256_add_ps(vecA, vecB); // Perform vectorized addition
      _mm256_storeu_ps(
          &tDA[i + 8],
          result); // Store the result back to A[i + 8] to A[i + 15]
    }

    // return tDA;

    std::vector<float> DA(this->connector->nxy);
    reorderVector_classics_R_r(this->stack, DA, tDA);
    return DA;
  }

  std::vector<float> _compute_drainage_area_beta1() {
    std::vector<float> tA(this->connector->nxy,
                          this->connector->get_cell_area(0));
    float *ptrA = tA.data();
    int nremaining = this->connector->nxy % 8;
    __m256 DA = _mm256_set1_ps(this->connector->get_cell_area(0));
    __m256 m1 = _mm256_set1_ps(-1);
    __m256 p1 = _mm256_set1_ps(1);
    __m256i nodes;
    __m256i recs;
    __m256 ttAs;
    __m256 ttArecs;
    __m256i cmpResult;
    __m256 cmpResultf;
    __m256 tadd;

    // std::vector<int> argrec(this->connector->nxy);
    // for(int i=0;i<this->connector->nxy;++i)
    // {
    // 	this->argrec[i] =
    // this->argstack[this->connector->Sreceivers(this->stack[i])];
    // }

    // for(int i=this->connector->nxy-1;i>7; i -= 8)
    // {
    // 	// int node = this->stack[i];
    // 	// __m256i nodes = _mm256_loadu_si256(reinterpret_cast<const
    // __m256i*>(&this->stack[i+7])); 	nodes =
    // _mm256_setr_epi32(i-0,i-1,i-2,i-3,i-4,i-5,i-6,i-7); 	recs =
    // _mm256_setr_epi32(argrec[i-0],argrec[i-1],argrec[i-2],argrec[i-3],argrec[i-4],argrec[i-5],argrec[i-6],argrec[i-7]);
    // 	// printM256i(nodes);
    // 	// printM256i(recs);
    // 	ttAs = _mm256_i32gather_ps(ptrA, nodes, sizeof(float));
    // 	ttArecs = _mm256_set1_ps(0.);
    // 	ttAs = _mm256_add_ps(ttAs, DA);

    // 	// if(i > this->connector->nx * 10)
    // 	// {
    // 	// 	printM256i(nodes);
    // 	// 	printM256i(recs);
    // 	// }
    // 	cmpResult = _mm256_cmpeq_epi32(nodes, recs);
    // 	// if(i > this->connector->nx * 10) printM256i(cmpResult);

    // 	cmpResultf = _mm256_cvtepi32_ps(cmpResult);
    // 	cmpResultf = _mm256_add_ps(cmpResultf,p1);

    // 	tadd = _mm256_mul_ps(ttAs,cmpResultf);
    // 	// if(i > this->connector->nx * 10) printM256(cmpResultf);
    // 	// if(i > this->connector->nx * 10) printM256(tadd);
    // 	// if(i > this->connector->nx * 10) printM256(ttArecs);
    // 	ttArecs = _mm256_add_ps(tadd,ttArecs);
    // 	// if(i > this->connector->nx * 10) printM256(ttArecs);

    // 	// if(i > this->connector->nx * 10) std::cout << std::endl;

    // 	alignas(32) float rtA[8], rtArec[8];
    // 	alignas(32) int rnodes[8], rrecs[8];

    // 	_mm256_store_ps(rtA, ttAs);
    // 	_mm256_store_ps(rtArec, ttArecs);
    // 	_mm256_storeu_si256((__m256i*)rnodes,nodes);
    // 	_mm256_storeu_si256((__m256i*)rrecs,recs);
    // 	// std::cout << "J" << std::endl;

    // 	// printM256i(nodes);
    // 	// printM256i(recs);

    // 	// printM256(ttAs);
    // 	// printM256(ttArecs);

    // 	tA[rnodes[0]] = rtA[0];
    // 	tA[rnodes[1]] = rtA[1];
    // 	tA[rnodes[2]] = rtA[2];
    // 	tA[rnodes[3]] = rtA[3];
    // 	tA[rnodes[4]] = rtA[4];
    // 	tA[rnodes[5]] = rtA[5];
    // 	tA[rnodes[6]] = rtA[6];
    // 	tA[rnodes[7]] = rtA[7];
    // 	tA[rrecs[0]] += rtArec[0];
    // 	tA[rrecs[1]] += rtArec[1];
    // 	tA[rrecs[2]] += rtArec[2];
    // 	tA[rrecs[3]] += rtArec[3];
    // 	tA[rrecs[4]] += rtArec[4];
    // 	tA[rrecs[5]] += rtArec[5];
    // 	tA[rrecs[6]] += rtArec[6];
    // 	tA[rrecs[7]] += rtArec[7];

    // 	// tA[node] += this->connector->get_cell_area(node);
    // 	// int rec = this->connector->Sreceivers(node);
    // 	// if(rec != node)
    // 	// {
    // 	// 	tA[rec] += tA[node];
    // 	// }

    // }
    // for(int i=0;i<this->connector->nxy;++i)
    for (int i = this->connector->nxy - 1; i > 7; i -= 8) {
      // if (i != argrec[i])
      ttAs = _mm256_setr_ps(tA[i - 0], tA[i - 1], tA[i - 2], tA[i - 3],
                            tA[i - 4], tA[i - 5], tA[i - 6], tA[i - 7]);
      ttArecs = _mm256_setr_ps(tA[argrec[i - 0]], tA[argrec[i - 1]],
                               tA[argrec[i - 2]], tA[argrec[i - 3]],
                               tA[argrec[i - 4]], tA[argrec[i - 5]],
                               tA[argrec[i - 6]], tA[argrec[i - 7]]);
      ttArecs = _mm256_add_ps(ttAs, ttArecs);

      tA[argrec[i]] += tA[i];
      alignas(32) float rtA[8];
      _mm256_store_ps(rtA, ttArecs);
      tA[argrec[i - 0]] = rtA[0];
      tA[argrec[i - 1]] = rtA[1];
      tA[argrec[i - 2]] = rtA[2];
      tA[argrec[i - 3]] = rtA[3];
      tA[argrec[i - 4]] = rtA[4];
      tA[argrec[i - 5]] = rtA[5];
      tA[argrec[i - 6]] = rtA[6];
      tA[argrec[i - 7]] = rtA[7];
    }

    // std::vector<float> ttA(this->connector->nxy,0.);
    // for(int i=0; i<this->connector->nxy;++i)
    // {
    // 	ttA[this->argstack[i]] = tA[this->stack[i]];
    // }

    // return ttA;
    return tA;
  }

  std::vector<float> _compute_drainage_area_classic() {
    std::vector<float> tA(this->connector->nxy, 0.);
    for (int i = this->connector->nxy - 1; i >= 0; --i) {
      int node = this->stack[i];
      tA[node] += this->connector->get_cell_area(node);
      int rec = this->connector->Sreceivers(node);
      if (rec != node) {
        tA[rec] += tA[node];
      }
    }
    return tA;
  }

  void save_drainage_area(std::string filename) {
    auto tA = this->_compute_drainage_area_classic();
    const std::vector<std::uint64_t> shape{std::uint64_t(this->connector->ny),
                                           std::uint64_t(this->connector->nx)};
    const bool fortran_order{false};
    const std::string path{filename};
    npy::SaveArrayAsNumpy(path, fortran_order, shape.size(), shape.data(), tA);
  }
};

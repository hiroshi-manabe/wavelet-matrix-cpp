/* 
 *  Copyright (c) 2010 Daisuke Okanohara
  * 
 *   Redistribution and use in source and binary forms, with or without
 *   modification, are permitted provided that the following conditions
 *   are met:
 * 
 *   1. Redistributions of source code must retain the above Copyright
 *      notice, this list of conditions and the following disclaimer.
 *
 *   2. Redistributions in binary form must reproduce the above Copyright
 *      notice, this list of conditions and the following disclaimer in the
 *      documentation and/or other materials provided with the distribution.
 *
 *   3. Neither the name of the authors nor the names of its contributors
 *      may be used to endorse or promote products derived from this
 *      software without specific prior written permission.
 */

#include <gtest/gtest.h>
#include <vector>
#include <string>
#include <sstream>
#include "../src/bit_array.hpp"

using namespace std;
using namespace wat_array;

TEST(bitvec, trivial){
  BitArray ba;
  ASSERT_EQ(0, ba.length());
  ASSERT_EQ(0, ba.one_num());
  
  ostringstream oss;
  ba.Save(oss);
  ASSERT_EQ(false, !oss);
  
  istringstream iss(oss.str());
  BitArray ba_load;
  ba_load.Load(iss);
  ASSERT_EQ(false, !iss);

  ASSERT_EQ(0, ba_load.length());
  ASSERT_EQ(0, ba_load.one_num());
}

TEST(bitvec, selectblock){
  uint64_t x = 0;

  for (uint64_t i = 0; i < 64; ++i){
    ASSERT_EQ(i, BitArray::SelectInBlock(~x, i+1));
  }

  for (uint64_t i = 0; i < 64; ++i){
    x |= (1LLU << i);
  }

  for (uint64_t i = 0; i < 64; ++i){
    ASSERT_EQ(i, BitArray::SelectInBlock(x, i+1));
  }
}

TEST(bitvec, trivial_zero){
  const int N = 100;
  BitArray ba(N);
  
  ba.Build();
  ASSERT_EQ(N, ba.length());
  for (size_t i = 0; i < ba.length(); ++i){
    ASSERT_EQ(0  , ba.Lookup(i));
    ASSERT_EQ(i  , ba.Rank(0, i));
    ASSERT_EQ(i  , ba.Select(0, i+1));
  }
}

TEST(bitvec, trivial_one){
  const int N = 1000;
  BitArray ba(N);
  for (int i = 0; i < 1000; ++i){
    ba.SetBit(1, i);
  }

  ba.Build();
  ASSERT_EQ(N, ba.length());
  for (size_t i = 0; i < ba.length(); ++i){
    ASSERT_EQ(1  , ba.Lookup(i));
    ASSERT_EQ(i  , ba.Rank(1, i));
    ASSERT_EQ(i  , ba.Select(1, i+1));
  }
}

TEST(bitvec, random){
  const int N = 100000;
  BitArray ba(N);
  vector<int> B;
  for (int i = 0; i < N; ++i){
    int b = rand() % 2;
    ba.SetBit(b, i);
    B.push_back(b);
  }
  
  ba.Build();
  ASSERT_EQ(N, ba.length());
  int sum = 0;
  for (size_t i = 0; i < ba.length(); ++i){
    ASSERT_EQ(B[i], ba.Lookup(i));
    if (B[i]){
      ASSERT_EQ(sum, ba.Rank(1, i));
      EXPECT_EQ(i,   ba.Select(1, sum+1));
    } else {
      ASSERT_EQ(i - sum, ba.Rank(0, i));
      EXPECT_EQ(i,       ba.Select(0, i-sum+1));
    }
    sum += B[i];
  }
}

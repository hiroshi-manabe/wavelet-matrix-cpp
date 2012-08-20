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
#include <sstream>
#include <algorithm>
#include "../src/wat_array.hpp"

using namespace std;

struct RandomQuery{
  RandomQuery(int n){
    for (;;){
      beg = rand() % n;
      end = rand() % (n+1);
      if (beg != end) break;
    }
    if (beg > end) swap(beg, end);
  }
  uint64_t beg;
  uint64_t end;
};

TEST(wat_array, trivial){
  wat_array::WatArray wa;
  ASSERT_EQ(0, wa.alphabet_num());
  ASSERT_EQ(0, wa.length());
  ASSERT_EQ(wat_array::NOTFOUND, wa.Rank(0, 0));
  ASSERT_EQ(wat_array::NOTFOUND, wa.Select(0, 0));
  ASSERT_EQ(wat_array::NOTFOUND, wa.Lookup(0));
  ASSERT_EQ(wat_array::NOTFOUND, wa.Freq(0));
  ASSERT_EQ(wat_array::NOTFOUND, wa.FreqSum(0, 1));
  
  uint64_t rank = 0;
  uint64_t rank_less_than = 0;
  uint64_t rank_more_than = 0;
  wa.RankAll(0, 0, rank, rank_less_than, rank_more_than);
  ASSERT_EQ(wat_array::NOTFOUND, rank);
  ASSERT_EQ(wat_array::NOTFOUND, rank_more_than);
  ASSERT_EQ(wat_array::NOTFOUND, rank_less_than);

  uint64_t pos = 0;
  uint64_t val = 0;
  wa.MaxRange(0, 0, pos, val);
  ASSERT_EQ(wat_array::NOTFOUND, pos);
  wa.MinRange(0, 0, pos, val);
  ASSERT_EQ(wat_array::NOTFOUND, pos);

  ostringstream oss;
  wa.Save(oss);
  ASSERT_EQ(false, !oss);
  
  istringstream iss(oss.str());
  wat_array::WatArray ws_load;
  ws_load.Load(iss);
  ASSERT_EQ(false, !iss);
  
  ASSERT_EQ(0, ws_load.alphabet_num());
  ASSERT_EQ(0, ws_load.length());
  ASSERT_EQ(wat_array::NOTFOUND, ws_load.Rank(0, 0));
  ASSERT_EQ(wat_array::NOTFOUND, ws_load.Select(0, 0));
  ASSERT_EQ(wat_array::NOTFOUND, ws_load.Lookup(0));
}

TEST(wat_array, alphanum_one){
  vector<uint64_t> A;
  A.push_back(0);
  A.push_back(0);
  A.push_back(0);
  A.push_back(0);
  A.push_back(0);



  wat_array::WatArray wa;
  wa.Init(A);
  ASSERT_EQ(5, wa.length());
  ASSERT_EQ(1, wa.alphabet_num());
  ASSERT_EQ(5, wa.Freq(0));
  ASSERT_EQ(5, wa.FreqSum(0, 1));
  for (uint64_t i = 0; i < wa.length(); ++i){
    ASSERT_EQ(i, wa.Rank(0, i));
    ASSERT_EQ(0, wa.RankLessThan(0, i));
    ASSERT_EQ(0, wa.RankMoreThan(0, i));
    ASSERT_EQ(i, wa.Select(0, i+1));
    for (uint64_t j = i+1; j <= wa.length(); ++j){
      uint64_t pos = 0;
      uint64_t val = 0;
      wa.MaxRange(i, j, pos, val);
      ASSERT_EQ(i, pos);
      ASSERT_EQ(0, val);

      pos = 0;
      val = 0;
      wa.MinRange(i, j, pos, val);
      ASSERT_EQ(i, pos);
      ASSERT_EQ(0, val);

      vector<wat_array::ListResult> lrs;
      wa.ListMinRange(0, 1, i, j, j-i, lrs);
      ASSERT_EQ(1, lrs.size());
      ASSERT_EQ(0, lrs[0].c);
      ASSERT_EQ(j-i, lrs[0].freq);

      wa.ListModeRange(0, 1, i, j, j-i, lrs);
      sort(lrs.begin(), lrs.end());
      ASSERT_EQ(1, lrs.size());
      ASSERT_EQ(0, lrs[0].c);
      ASSERT_EQ(j-i, lrs[0].freq);
    }
  }
}

TEST(wat_array, save){
  vector<uint64_t> A;
  A.push_back(1);
  A.push_back(1);
  A.push_back(1);
  A.push_back(1);
  A.push_back(3);
  A.push_back(0);
  A.push_back(1);
  A.push_back(1);
  A.push_back(1);
  A.push_back(1);
  A.push_back(3);
  A.push_back(1);
  A.push_back(2);
  A.push_back(2);
  wat_array::WatArray wa;
  wa.Init(A);
  ASSERT_EQ(2, wa.Rank(3, 14));
  ostringstream os;
  wa.Save(os);
  istringstream is(os.str());
  wat_array::WatArray wa_load;
  wa_load.Load(is);
  ASSERT_EQ(2, wa_load.Rank(3, 14));
}
  
TEST(wat_array, small){

  vector<uint64_t> array;
  const uint64_t alphabet_num = 200;
  for (uint64_t i = 0; i < alphabet_num; ++i){
    array.push_back(i);
  }
  const uint64_t length = array.size();

  wat_array::WatArray wa;
  wa.Init(array);
  ASSERT_EQ(alphabet_num, wa.alphabet_num());
  ASSERT_EQ(length, wa.length());

  for (uint64_t i = 0; i < alphabet_num; ++i){
    ASSERT_EQ(1, wa.Freq(i));
    for (uint64_t j = i; j < alphabet_num; ++j){
      ASSERT_EQ(j-i, wa.FreqSum(i, j));
    }
  }

  ASSERT_EQ(1, wa.Rank(wa.alphabet_num()-1, wa.length()));
  ASSERT_EQ(wat_array::NOTFOUND, wa.Rank(wa.alphabet_num()  , wa.length()));
  ASSERT_EQ(wat_array::NOTFOUND, wa.Rank(wa.alphabet_num()+1, wa.length()));

  vector<uint64_t> counts(alphabet_num);
  for (uint64_t i = 0; i < length; ++i){
    uint64_t c = array[i];
    ASSERT_EQ(c, wa.Lookup(i));
    uint64_t sum = 0;
    for (uint64_t j = 0; j < alphabet_num; ++j){
      ASSERT_EQ(counts[j],                      wa.Rank(j, i));
      ASSERT_EQ(sum,                            wa.RankLessThan(j, i));
      ASSERT_EQ(i - sum - counts[j],            wa.RankMoreThan(j, i));
      sum += counts[j];
    }
    counts[c]++;
    ASSERT_EQ(i, wa.Select(c, counts[c]));

    for (uint64_t j = i+1; j <= wa.length(); ++j){
      uint64_t pos = 0;
      uint64_t val = 0;
      wa.MaxRange(i, j, pos, val);
      ASSERT_EQ(j - 1, pos);
      ASSERT_EQ(j - 1, val);

      pos = 0;
      val = 0;
      wa.MinRange(i, j, pos, val);
      ASSERT_EQ(i, pos);
      ASSERT_EQ(i, val);

      vector<wat_array::ListResult> lrs;
      wa.ListMinRange(0, alphabet_num, i, j, j-i, lrs);
      
      for (size_t k = 0; k < lrs.size(); ++k){
	ASSERT_EQ(i+k, lrs[k].c);
	ASSERT_EQ(1,   lrs[k].freq);
      }

      wa.ListModeRange(0, alphabet_num, i, j, j-i, lrs);
      sort(lrs.begin(), lrs.end());
      for (size_t k = 0; k < lrs.size(); ++k){
	ASSERT_EQ(i+k, lrs[k].c);
	ASSERT_EQ(1,   lrs[k].freq);
      }
    }
  }
}

TEST(wat_array, random){
  vector<uint64_t> array;

  uint64_t alphabet_num = 100;
  uint64_t n = 10000;
  vector<uint64_t> freq(alphabet_num);
  for (uint64_t i = 0; i < n; ++i){
    uint64_t c = rand() % alphabet_num;
    array.push_back(c);
    freq[c]++;
  }

  wat_array::WatArray wa;
  wa.Init(array);

  ASSERT_EQ(alphabet_num, wa.alphabet_num());
  ASSERT_EQ(n, wa.length());
  for (uint64_t i = 0; i < alphabet_num; ++i){
    ASSERT_EQ(freq[i], wa.Freq(i));
  }

  vector<uint64_t> counts(alphabet_num);
  for (uint64_t i = 0; i < array.size(); ++i){
    uint64_t c = array[i];

    ASSERT_EQ(c, wa.Lookup(i));
    uint64_t sum = 0;
    for (uint64_t j = 0; j < alphabet_num; ++j){
      if ((rand() % 100) == 0){
	ASSERT_EQ(counts[j],           wa.Rank(j, i));
	ASSERT_EQ(sum,                 wa.RankLessThan(j, i));
	ASSERT_EQ(i - sum - counts[j], wa.RankMoreThan(j, i));
      }
      sum += counts[j];
    }
    counts[c]++;
    ASSERT_EQ(i, wa.Select(c, counts[c]));
  }
}


void SetVals(const RandomQuery& rq,
	     const vector<uint64_t>& array,
	     vector<pair<uint64_t, size_t> >& vals){
  for (size_t i = rq.beg; i < rq.end; ++i){
    vals.push_back(make_pair(array[i], i));
  }
  sort(vals.begin(), vals.end());
}

void UniqCount(const vector<pair<uint64_t, size_t> >& vals,
	       vector<pair<uint64_t, uint64_t> >& ret){
  if (vals.size() == 0) return;
  uint64_t prev = vals[0].first;
  uint64_t count = 1;
  for (size_t i = 1; i < vals.size(); i++){
    if (prev != vals[i].first){
      ret.push_back(make_pair(prev, count));
      prev = vals[i].first;
      count = 1;
    } else {
      count++;
    }
  }
  ret.push_back(make_pair(prev, count));
}

void FilterRange(const RandomQuery& char_range, vector<pair<uint64_t, uint64_t> >& uniq_counts){
  vector<pair<uint64_t, uint64_t> > new_uniq_counts;
  for (size_t i = 0; i < uniq_counts.size(); ++i){
    uint64_t val = uniq_counts[i].first;
    if (char_range.beg <= val && char_range.end > val){
      new_uniq_counts.push_back(uniq_counts[i]);
    }
  }
  uniq_counts.swap(new_uniq_counts);
}

class FreqComp{
public:
  bool operator () (const pair<uint64_t, uint64_t>& left,
		    const pair<uint64_t, uint64_t>& right) const{
    if (left.second != right.second) return left.second > right.second;
    return left.first < right.first;
  }
};

class FreqCompLR{
public:
  bool operator () (const wat_array::ListResult& left,
		    const wat_array::ListResult& right) const{
    if (left.freq != right.freq) return left.freq > right.freq;
    return left.c < right.c;
  }
};

void WatRandomInitialize(wat_array::WatArray& wa,
			 vector<uint64_t>& array,
			 uint64_t alphabet_num, 
			 uint64_t n){
  vector<uint64_t> freq(alphabet_num);
  for (uint64_t i = 0; i < n; ++i){
    uint64_t c = rand() % alphabet_num;
    array.push_back(c);
    freq[c]++;
  }
  wa.Init(array);
}


TEST(wat_array, min_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 100, 1000);

  for (size_t iter = 0; iter < 10; ++iter){
    RandomQuery rq(wa.length());
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);

    uint64_t min_pos = 0;
    uint64_t min_val = 0;
    wa.MinRange(rq.beg, rq.end, min_pos, min_val);
    ASSERT_EQ(vals.front().first , min_val); 
    ASSERT_EQ(vals.front().second, min_pos);
  }
}

TEST(wat_array, quantile_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 100, 1000);

  for (size_t iter = 0; iter < 10; ++iter){
    RandomQuery rq(wa.length());
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);

    uint64_t kth_pos = 0;
    uint64_t kth_val = 0;
    uint64_t k = rand() % (rq.end - rq.beg);
    wa.QuantileRange(rq.beg, rq.end, k, kth_pos, kth_val);
    ASSERT_EQ(vals[k].first, kth_val);
  }
}

TEST(wat_array, max_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 10, 1000);

  for (size_t iter = 0; iter < 10; ++iter){
    RandomQuery rq(wa.length());
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);

    uint64_t max_pos = 0;
    uint64_t max_val = 0;
    wa.MaxRange(rq.beg, rq.end, max_pos, max_val);
    ASSERT_EQ(vals.back().first , max_val);
  }
}

TEST(wat_array, freq_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 10, 1000);

  for (size_t iter = 0; iter < 10; ++iter){
    RandomQuery rq(wa.length());
    RandomQuery arq(wa.alphabet_num());
    uint64_t count = 0;
    for (uint64_t i = rq.beg; i < rq.end; ++i){
      if (arq.beg <= array[i] &&
	  array[i] < arq.end){
	++count;
      }
    }
    ASSERT_EQ(count, wa.FreqRange(arq.beg, arq.end, rq.beg, rq.end));
  }
}


TEST(wat_array, list_mode_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 100, 100);

  for (size_t iter = 0; iter < 10; ++iter){
    RandomQuery rq(wa.length());
    RandomQuery arq(wa.alphabet_num());
    
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);

    vector<pair<uint64_t, uint64_t> > uniq_counts;
    UniqCount(vals, uniq_counts);
    sort(uniq_counts.begin(), uniq_counts.end(), FreqComp());
    FilterRange(arq, uniq_counts);

    vector<wat_array::ListResult> lrs;
    uint64_t num = rq.end - rq.beg;
    wa.ListModeRange(arq.beg, arq.end, rq.beg, rq.end, num, lrs);
    sort(lrs.begin(), lrs.end(), FreqCompLR());

    for (size_t i = 0; i < lrs.size() && i < uniq_counts.size(); ++i){
      ASSERT_EQ(uniq_counts[i].first, lrs[i].c);
      ASSERT_EQ(uniq_counts[i].second, lrs[i].freq);
    }
  }
}

TEST(wat_array, list_min_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 100, 1000);

  for (size_t iter = 0; iter < 10; ++iter){ 
    RandomQuery rq(wa.length());
    RandomQuery arq(wa.alphabet_num());
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);
    vector<pair<uint64_t, uint64_t> > uniq_counts;
    UniqCount(vals, uniq_counts);
    FilterRange(arq, uniq_counts);

    vector<wat_array::ListResult> lrs;
    uint64_t num = rq.end - rq.beg;
    wa.ListMinRange(arq.beg, arq.end, rq.beg, rq.end, num, lrs);
    for (size_t i = 0; i < lrs.size() && i < uniq_counts.size(); ++i){
      ASSERT_EQ(uniq_counts[i].first, lrs[i].c);
      ASSERT_EQ(uniq_counts[i].second, lrs[i].freq);
    }
  }
}

TEST(wat_array, list_max_range){
  wat_array::WatArray wa;
  vector<uint64_t> array;
  WatRandomInitialize(wa, array, 100, 1000);

  for (size_t iter = 0; iter < 10; ++iter){ 
    RandomQuery rq(wa.length());
    RandomQuery arq(wa.alphabet_num());
    vector<pair<uint64_t, size_t> > vals;
    SetVals(rq, array, vals);
    reverse(vals.begin(), vals.end());
    vector<pair<uint64_t, uint64_t> > uniq_counts;
    UniqCount(vals, uniq_counts);
    FilterRange(arq, uniq_counts);

    vector<wat_array::ListResult> lrs;
    uint64_t num = rq.end - rq.beg;
    wa.ListMaxRange(arq.beg, arq.end, rq.beg, rq.end, num, lrs);
    for (size_t i = 0; i < lrs.size() && i < uniq_counts.size(); ++i){
      ASSERT_EQ(uniq_counts[i].first, lrs[i].c);
      ASSERT_EQ(uniq_counts[i].second, lrs[i].freq);
    }
  }
}




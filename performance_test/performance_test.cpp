#include <iomanip>
#include <stdlib.h>
#include <time.h>
#include <algorithm>
#include <sys/time.h>
#include <stdio.h>

#include "../src/wavelet_matrix.hpp"

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

double gettimeofday_sec() {
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

void SortUniqCount(const vector<uint64_t>& array, 
		   uint64_t c_min, uint64_t c_max,
		   uint64_t beg, uint64_t end,
		   vector<pair<uint64_t, uint64_t> >& res){
  if (array.size() == 0) return;
  vector<uint64_t> vals; 
  for (uint64_t i = beg; i < end; ++i){
    if (array[i] >= c_min &&
	array[i] <  c_max){
      vals.push_back(array[i]);
    }
  }
  sort(vals.begin(), vals.end());

  uint64_t prev = vals[0];
  uint64_t count = 1;
  for (size_t i = 0; i < vals.size(); ++i){
    if (vals[i] != prev){
      res.push_back(make_pair(count, prev));
      prev = vals[i];
      count = 1;
    } else {
      count++;
    }
  }
  res.push_back(make_pair(count, prev));
}


struct QuerySet {
  QuerySet(int iter_num, uint64_t length, uint64_t alphabet_num) : 
    iter_num(iter_num),
    length(length),
    alphabet_num(alphabet_num),
    non_zero_num(0),
    array(length),
    freqs(alphabet_num), 
    pos_queries(iter_num),
    char_queries(iter_num),
    select_char_queries(iter_num),
    select_rank_queries(iter_num),
    range_queries(iter_num, length),
    range_alpha_queries(iter_num, alphabet_num) {

    for (uint64_t i = 0; i < length; ++i){
      array[i] = rand() % alphabet_num;
      freqs[array[i]]++;
    }

    for (size_t i = 0; i < freqs.size(); ++i){
      if (freqs[i] > 0) non_zeros.push_back(i);
    }
    non_zero_num = non_zeros.size();
    
    for (int i = 0; i < iter_num; ++i){
      pos_queries[i] = rand() % (length+1);
      char_queries[i] = rand() % alphabet_num;
      uint64_t select_char = non_zeros[rand() % non_zero_num];
      select_char_queries[i] = select_char;
      select_rank_queries[i] = (rand() % freqs[select_char]) + 1; 
    }
  }

  int iter_num;  
  uint64_t length;
  uint64_t alphabet_num;
  uint64_t non_zero_num;
  vector<uint64_t> array;
  vector<uint64_t> freqs;
  vector<uint64_t> non_zeros;


  vector<uint64_t> pos_queries;
  vector<uint64_t> char_queries;
  vector<uint64_t> select_char_queries;
  vector<uint64_t> select_rank_queries;
  vector<RandomQuery> range_queries;
  vector<RandomQuery> range_alpha_queries;
};

void TestBaseline(QuerySet& qs){
  double begin_time = 0.0;

  double init_time = 0.0;
  
  uint64_t dummy = 0;
  int iter_num = qs.iter_num;
  vector<uint64_t>& array = qs.array;

  // lookup
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    dummy += array[qs.pos_queries[i]];
  }
  double lookup_time = gettimeofday_sec() - begin_time;

  // rank
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    uint64_t rank = 0;
    uint64_t c = qs.char_queries[i]; 
    for (uint64_t j = 0; j < qs.pos_queries[i]; ++j){
      if (array[j] == c) ++rank;
    }
    dummy += rank;
  }
  double rank_time = gettimeofday_sec() - begin_time;

  // select
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    uint64_t rank = qs.select_rank_queries[i];
    uint64_t c = qs.select_char_queries[i];
    uint64_t pos = 0;
    for (uint64_t j = 0; j < array.size(); ++j){
      if (array[j] == c){
	--rank;
	if (rank == 0) {
	  pos = j;
	  break;
	} 
      }
    }
    dummy += pos;
  }
  double select_time = gettimeofday_sec() - begin_time;

  // max_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    uint64_t max_val = 0;
    uint64_t max_pos = 0;
    for (uint64_t j = qs.range_queries[i].beg; j < qs.range_queries[i].end; ++j){
      if (array[j] > max_val){
	max_val = array[j];
	max_pos = j;
      }
    }
    dummy += max_val;
  }
  double max_range_time = gettimeofday_sec() - begin_time;

  // quantile_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    RandomQuery& rq = qs.range_queries[i];
    vector<pair<uint64_t, uint64_t> > vals;
    for (uint64_t j = rq.beg; j < rq.end; ++j){
      vals.push_back(make_pair(array[j], j));
    }
    partial_sort(vals.begin(), vals.begin() + (rq.end - rq.beg) / 2, vals.end());
    dummy += vals[(rq.end - rq.beg) / 2].first;
  }
  double quantile_range_time = gettimeofday_sec() - begin_time;

  // freq_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    RandomQuery& rq = qs.range_queries[i];
    RandomQuery& arq = qs.range_alpha_queries[i];
    uint64_t freq = 0;
    for (uint64_t j = rq.beg; j < rq.end; ++j){
      if (array[j] >= arq.beg && array[j] < arq.end){
	++freq;
      }
    }
    dummy += freq;
  }
  double freq_range_time = gettimeofday_sec() - begin_time;

  double ratio_micro = 1.0 / qs.iter_num * 1000000.0;
  cerr  << scientific<< qs.length   << "\t"
        << scientific<< qs.alphabet_num  << "\t"
        << scientific<< init_time  << "\t"
        << scientific<< lookup_time  * ratio_micro << "\t"
        << scientific<< rank_time  * ratio_micro << "\t"
        << scientific<< select_time  * ratio_micro << "\t"
        << scientific<< max_range_time  * ratio_micro << "\t"
        << scientific<< quantile_range_time  * ratio_micro << "\t" 
        << scientific<< freq_range_time  * ratio_micro << endl;

  if (dummy == 7777) {
    cerr << ""; // remove optimization
  }
}


void TestWaveletMatrix(QuerySet& qs){
  wavelet_matrix::WaveletMatrix ws;
  double begin_time = 0.0;

  begin_time = gettimeofday_sec();
  ws.Init(qs.array);
  double init_time = gettimeofday_sec() - begin_time;
  

  uint64_t dummy = 0;
  int iter_num = qs.iter_num;
  // lookup
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    dummy += ws.Lookup(qs.pos_queries[i]);
  }
  double lookup_time = gettimeofday_sec() - begin_time;

  // rank
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    dummy += ws.Rank(qs.char_queries[i], qs.pos_queries[i]);
  }
  double rank_time = gettimeofday_sec() - begin_time;

  // select
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    dummy += ws.Select(qs.select_char_queries[i], qs.select_rank_queries[i]);
  }
  double select_time = gettimeofday_sec() - begin_time;

  // max_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    uint64_t pos = 0;
    uint64_t val = 0;
    ws.MaxRange(qs.range_queries[i].beg, qs.range_queries[i].end, pos, val);
    dummy += pos;
  }
  double max_range_time = gettimeofday_sec() - begin_time;

  // quantile_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    RandomQuery& rq = qs.range_queries[i];
    uint64_t pos = 0;
    uint64_t val = 0;
    ws.QuantileRange(rq.beg, rq.end, (rq.end-rq.beg)/2, pos, val);
    dummy += pos;
  }
  double quantile_range_time = gettimeofday_sec() - begin_time;

  // freq_range
  begin_time = gettimeofday_sec();
  for (int i = 0; i < iter_num; ++i){
    RandomQuery& rq = qs.range_queries[i];
    RandomQuery& arq = qs.range_alpha_queries[i];
    dummy += ws.FreqRange(arq.beg, arq.end, rq.beg, rq.end);
  }
  double freq_range_time = gettimeofday_sec() - begin_time;

  double ratio_micro = 1.0 / qs.iter_num * 1000000.0;
  cerr  << scientific<< qs.length  << "\t"
        << scientific<< qs.alphabet_num  << "\t"
        << scientific<< init_time  << "\t"
        << scientific<< lookup_time  * ratio_micro << "\t"
        << scientific<< rank_time  * ratio_micro << "\t"
        << scientific<< select_time  * ratio_micro << "\t"
        << scientific<< max_range_time  * ratio_micro << "\t"
        << scientific<< quantile_range_time  * ratio_micro << "\t" 
        << scientific<< freq_range_time  * ratio_micro << endl;
  if (dummy == 7777) cerr << "";
}

int main(int argc, char* argv[]){
  Test();

  cerr << "Performance Test init=total_time(sec.) other=avg_time(micro sec.) " << endl;
  cerr  << "method"  << "\t"
	<< "length"  << "\t"
        << "alnum"  << "\t"
        << "init"  << "\t"
        << "lookup"  << "\t"
        << "rank"  << "\t"
        << "select"  << "\t"
        << "max_range"  << "\t"
        << "quan_range"  << "\t" 
        << "freq_range"  << endl;

  for (uint64_t length = 1000; length <= 100000000; length *= 10){
    for (uint64_t alphabet_num = 10; alphabet_num <= length ; alphabet_num *= 100){
      QuerySet qs(100, length, alphabet_num);
      cerr << "ws "; TestWaveletMatrix(qs);
    }
  }

  return 0;

  for (uint64_t length = 1000; length <= 100000000; length *= 10){
    for (uint64_t alphabet_num = 10; alphabet_num <= length ; alphabet_num *= 100){
      QuerySet qs(100, length, alphabet_num);
      cerr << "bs "; TestBaseline(qs);
    }
  }


      
  return 0;
}

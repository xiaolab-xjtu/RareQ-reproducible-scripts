#include <Rcpp.h>
using namespace Rcpp;

// 辅助函数：聚类标签转为0-based
IntegerVector cluster_to_0based(IntegerVector clusters) {
  int n = clusters.size();
  IntegerVector unique_cl = unique(clusters);
  int k = unique_cl.size();
  
  std::map<int, int> cl_map;
  for (int i = 0; i < k; ++i) {
    cl_map[unique_cl[i]] = i;
  }
  
  IntegerVector res(n);
  for (int i = 0; i < n; ++i) {
    res[i] = cl_map[clusters[i]];
  }
  return res;
}

// [[Rcpp::export]]
NumericVector silhouette_rcpp_final(NumericMatrix data, IntegerVector clusters, int block_size = 1000) {
  int n = data.rows();
  int p = data.cols();
  
  // 处理聚类标签
  IntegerVector cl = cluster_to_0based(clusters);
  IntegerVector unique_cl = unique(cl);
  int k = unique_cl.size();
  if (k == 1) return NumericVector(n, 0.0);
  
  // 分块计算轮廓系数
  NumericVector sil_width(n);
  for (int b = 0; b < n; b += block_size) {
    int end = std::min(b + block_size, n);
    for (int i = b; i < end; ++i) {
      int xi = cl[i];
      
      // 计算a(i)
      double a = 0.0;
      int m = 0;
      for (int j = 0; j < n; ++j) {
        if (j == i || cl[j] != xi) continue;
        double dist = 0.0;
        for (int d = 0; d < p; ++d) {
          double diff = data(i, d) - data(j, d);
          dist += diff * diff;
        }
        a += sqrt(dist);
        m++;
      }
      if (m == 0) {
        sil_width[i] = 0.0;
        continue;
      }
      a /= m;
      
      // 计算b(i)
      double min_b = INFINITY;
      for (int c = 0; c < k; ++c) {
        int curr_cl = unique_cl[c];
        if (curr_cl == xi) continue;
        
        double b_sum = 0.0;
        int cnt = 0;
        for (int j = 0; j < n; ++j) {
          if (cl[j] != curr_cl) continue;
          double dist = 0.0;
          for (int d = 0; d < p; ++d) {
            double diff = data(i, d) - data(j, d);
            dist += diff * diff;
          }
          b_sum += sqrt(dist);
          cnt++;
        }
        if (cnt == 0) continue;
        double b = b_sum / cnt;
        if (b < min_b) min_b = b;
      }
      
      sil_width[i] = (min_b - a) / std::max(a, min_b);
    }
    checkUserInterrupt();
  }
  
  return sil_width;
}
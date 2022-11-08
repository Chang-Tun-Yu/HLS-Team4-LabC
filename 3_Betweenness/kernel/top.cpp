#include "top.hpp"
#include "container.hpp"
#include <list>
#include <queue>
#include <stack>
#define __SYNTHESIS__
void BFS(
	unsigned numVert,
	float sigma[N],
	Stack<data_t, flag_t> &stack,
    unsigned offset[N+1],
    unsigned column[E],
	unsigned int p[N][L],
  unsigned i_source,
  unsigned &o_source,
	unsigned int cnt[N]
)
{
  // for (unsigned i=0; i < numVert; i++) {
    int   dist[N];
    static Queue<data_t, flag_t> queue;
    unsigned source = i_source;

    stack.reset();
    queue.reset();
    each_vert: for (int j = 0; j < numVert; j++)
    {
      sigma[j] = 0;
      dist[j] = -1;
      cnt[j] = 0;
    }
    sigma[source] = 1;
    dist[source] = 0;

    queue.push_back(source);
    each_child: while (!queue.empty())
    {
      unsigned v = queue.front();

      stack.push_back(v);
      unsigned start = offset[v];
      unsigned end = offset[v+1];

      int dist_v = dist[v];
      float sigma_v= sigma[v];
	  
      for (unsigned j = start; j < end; j++)
      {
#pragma HLS pipeline
#pragma HLS dependence variable=sigma type=inter false
#pragma HLS dependence variable=cnt type=inter false
#pragma HLS dependence variable=dist type=inter false
#pragma HLS dependence variable=p type=inter false

    	  unsigned w = column[j];
        int dist_w = dist[w];
        bool flag = false;
        if (dist_w < 0)
        {
          flag = true;
          queue.push_back(w);
          dist[w] = dist_v + 1;
        }
        if (dist_w == dist_v + 1 || flag)
        {
          unsigned cnt_tmp = cnt[w];
          sigma[w] = sigma[w] + sigma_v;
          p[w][cnt_tmp] = v;
          cnt[w] = cnt_tmp + 1;

        }
      }
      queue.pop_front();

    }
    o_source = source;
  // }
}
void BACK(
	unsigned numVert,
	float delta[N],
	float sigma[N],
	Stack<data_t, flag_t> &stack,
	float btwn[N+1],
  unsigned i,
	unsigned int p[N][L],
	unsigned int cnt[N]
)
{
  // for (unsigned i=0; i < numVert; i++) {
    unsigned source = i;
    init: for (int j = 0; j < numVert; j++)
    {
      delta[j] = 0;
    }

    each_vert: while (!stack.empty())
    {
      unsigned w = stack.back();

      if (source != w)
      {
        btwn[w] = btwn[w] + delta[w];
      }

      float sigma_w = sigma[w];
      float delta_w = delta[w];

      each_parents: for (int k = 0; k < cnt[w]; k++)
      {
#pragma HLS pipeline
#pragma HLS dependence variable=delta type=inter false
#pragma HLS dependence variable=sigma type=inter false
        unsigned v = p[w][k];

        delta[v] = delta[v] + (sigma[v] / sigma_w) * (1 + delta_w);
        // if (source != w) {
        //     btwn[w] = btwn[w] + delta[w];
        // }
      }
      //return;
      stack.pop_back();

    }
  // }
}

void wrapper(
    unsigned numVert,
    unsigned i,
    unsigned offset[N+1],
    unsigned column[E],
    float btwn[N]
 ) {
    #pragma HLS allocation instances=BFS limit=1 function
    #pragma HLS allocation instances=BACK limit=1 function
    #pragma HLS DATAFLOW

	float delta1[N]; //, delta2[N];
	float sigma1[N]; //, sigma2[N];
	Stack<data_t, flag_t> stack1; //, stack2;
	unsigned int p1[N][L]; //, p2[N][L];
	unsigned int cnt1[N]; //, cnt2[N];
	unsigned source_channel;
    BFS(numVert, sigma1, stack1, offset, column, p1, i, source_channel, cnt1);
	  BACK(numVert, delta1, sigma1, stack1, btwn, source_channel, p1, cnt1);
  // }
}

extern "C" void dut(
    unsigned numVert,
    unsigned numEdge,
    unsigned *offset,
    unsigned *column,
    float *btwn,
    unsigned *tmp0,
    unsigned *tmp1,
    unsigned *tmp2,
    unsigned *tmp3)
{
  // clang-format off

    const unsigned MEMSIZE=INTERFACE_MEMSIZE;
#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem0 port = offset depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem1 port = column depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 256 max_read_burst_length = 256 bundle = gmem3 port = btwn depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem6 port = tmp0 depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem7 port = tmp1 depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem8 port = tmp2 depth = MEMSIZE

#pragma HLS INTERFACE m_axi offset = slave latency = 32 num_write_outstanding = 1 num_read_outstanding = \
    16 max_write_burst_length = 2 max_read_burst_length = 256 bundle = gmem9 port = tmp3 depth = MEMSIZE

// clang-format on
#ifndef __SYNTHESIS__
  for (int i = 0; i < numVert; i++)
  {
    btwn[i] = 0;
  }
  for (int i = 0; i < numVert; i++)
  {
    std::stack<unsigned> s;
    std::vector<std::list<unsigned>> p(numVert);
    std::vector<float> sigma(numVert);
    std::vector<int> dist(numVert);
    std::queue<unsigned> q;
    unsigned source = i;

    for (int j = 0; j < numVert; j++)
    {
      sigma[j] = 0;
      dist[j] = -1;
    }
    sigma[source] = 1;
    dist[source] = 0;

    q.push(source);
    while (!q.empty())
    {
      unsigned v = q.front();
      s.push(v);
      for (int j = offset[v]; j < offset[v + 1]; j++)
      {
        unsigned w = column[j];
        if (dist[w] < 0)
        {
          q.push(w);
          dist[w] = dist[v] + 1;
        }
        if (dist[w] == dist[v] + 1)
        {
          sigma[w] = sigma[w] + sigma[v];
          p[w].push_back(v);
        }
      }
      q.pop();
    }

    std::vector<float> delta(numVert);
    for (int j = 0; j < numVert; j++)
    {
      delta[j] = 0;
    }
    while (!s.empty())
    {
      unsigned w = s.top();
      if (source != w)
      {
        btwn[w] = btwn[w] + delta[w];
      }
      for (std::list<unsigned>::iterator it = p[w].begin(); it != p[w].end();
           it++)
      {
        unsigned v = *it;
        delta[v] = delta[v] + (sigma[v] / sigma[w]) * (1 + delta[w]);
        // if (source != w) {
        //     btwn[w] = btwn[w] + delta[w];
        // }
      }
      s.pop();
    }
  }
#else

  float _btwn[N+1];
  unsigned _offset[N+1];
  unsigned _column[E];
  read_data: for (int i = 0; i < numEdge; i++)
  {
    #pragma HLS PIPELINE II=1
    if (i <= numVert ) {
      _offset[i] = offset[i];
      _btwn[i] = 0;
      _column[i] = column[i];
    }
    else {
      _column[i] = column[i];
    }    
  }

	float delta1[N], delta2[N];
	float sigma1[N], sigma2[N];
	static Stack<data_t, flag_t> stack1, stack2;
	unsigned int p1[N][L], p2[N][L];
	unsigned int cnt1[N], cnt2[N];

    BFS(numVert, sigma1, stack1, _offset, _column, 0, p1, cnt1);

  each_source: for (int i = 0; i < numVert; i++)
  {
    //wrapper(numVert, i,_offset, _column, _btwn);
	if(i%2 == 0)
	{
		BACK(numVert, delta1, sigma1, stack1, _btwn, i, p1, cnt1);
        BFS(numVert, sigma2, stack2, _offset, _column, i+1, p2, cnt2);
	}
	else
	{
        BACK(numVert, delta2, sigma2, stack2, _btwn, i, p2, cnt2);
		BFS(numVert, sigma1, stack1, _offset, _column, i+1, p1, cnt1);
	}

  }
  
  if(numVert % 2 == 0)
    BACK(numVert, delta2, sigma2, stack2, _btwn, numVert-1, p2, cnt2);
  else
    BACK(numVert, delta1, sigma1, stack1, _btwn, numVert-1, p1, cnt1);


  write_btwn: for (int i = 0; i < numVert; i++)
  {
    #pragma HLS PIPELINE II=1
    btwn[i] = _btwn[i];
  }
#endif
}

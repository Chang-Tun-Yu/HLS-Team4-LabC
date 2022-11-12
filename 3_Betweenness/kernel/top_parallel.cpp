#include "top.hpp"
#include "container.hpp"
#include <list>
#include <queue>
#include <stack>
#include <math.h>
#define UNROLL_FACTOR 4
void BFS(
	unsigned numVert,
	float sigma[N],
	Stack<data_t, flag_t> &stack,
    unsigned offset[N+1],
    unsigned column[E],
	unsigned i,
	unsigned int p[N][L],
	unsigned int cnt[N]
)
{
    int   dist[N];
    static Queue<data_t, flag_t> queue;
    unsigned source = i;

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
}
void BACK(
	unsigned numVert,
	float delta[N],
	float sigma[N],
	Stack<data_t, flag_t> &stack,
	float btwn[N+1],
	unsigned source,
	unsigned int p[N][L],
	unsigned int cnt[N]
)
{
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
}
void Wrapper(
	unsigned i,
	unsigned j,
	unsigned numVert,
	float sigma1[N],
	float sigma2[N],
	Stack<data_t, flag_t> &stack1,
	Stack<data_t, flag_t> &stack2,
  unsigned _offset[N+1],
	unsigned _column[E],
	unsigned int p1[N][L],
  unsigned int p2[N][L],
	unsigned int cnt1[N],
  unsigned int cnt2[N],
  float delta1[N],
  float delta2[N],
  float _btwn[N]
  )
{
  
	if(i%2 == 0)
	{
		BACK(numVert, delta1, sigma1, stack1, _btwn, i*UNROLL_FACTOR+j, p1, cnt1);
	  if((i+1)*UNROLL_FACTOR+j < numVert) BFS(numVert, sigma2, stack2, _offset, _column, (i+1)*UNROLL_FACTOR+j, p2, cnt2);
	}
	else
	{
	  BACK(numVert, delta2, sigma2, stack2, _btwn, i*UNROLL_FACTOR+j, p2, cnt2);
		if((i+1)*UNROLL_FACTOR+j < numVert) BFS(numVert, sigma1, stack1, _offset, _column, (i+1)*UNROLL_FACTOR+j, p1, cnt1);
	}	
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


const unsigned UNROLL_FACTOR_tmp=UNROLL_FACTOR;
//#pragma HLS allocation instances=BFS limit=UNROLL_FACTOR_tmp function
//#pragma HLS allocation instances=BACK limit=UNROLL_FACTOR_tmp function


  float _btwn[UNROLL_FACTOR][N+1];
  unsigned _offset[UNROLL_FACTOR][N+1];
  unsigned _column[UNROLL_FACTOR][E];
#pragma HLS array_partition variable=_btwn  dim=1
#pragma HLS array_partition variable=_column  dim=1
#pragma HLS array_partition variable=_offset  dim=1
  read_data: for (int i = 0; i < numEdge; i++)
  {
#pragma HLS UNROLL factor=UNROLL_FACTOR_tmp
    for(int j = 0; j < UNROLL_FACTOR; j++){
#pragma HLS PIPELINE II=1
        if (i <= numVert ) {
          _offset[j][i] = offset[i];
          _btwn[j][i] = 0;
          _column[j][i] = column[i];
        }
        else {
          _column[j][i] = column[i];
        }
    }
  }


	float delta1[UNROLL_FACTOR][N], delta2[UNROLL_FACTOR][N];
	float sigma1[UNROLL_FACTOR][N], sigma2[UNROLL_FACTOR][N];
	static Stack<data_t, flag_t> stack1[UNROLL_FACTOR], stack2[UNROLL_FACTOR];
	unsigned int p1[UNROLL_FACTOR][N][L], p2[UNROLL_FACTOR][N][L];
	unsigned int cnt1[UNROLL_FACTOR][N], cnt2[UNROLL_FACTOR][N];
//#pragma HLS bind_storage variable=stack1[0].storage type=RAM_2P
#pragma HLS array_partition variable=stack1  dim=1
#pragma HLS array_partition variable=stack2  dim=1
#pragma HLS array_partition variable=delta1  dim=1
#pragma HLS array_partition variable=delta2  dim=1
#pragma HLS array_partition variable=sigma1  dim=1
#pragma HLS array_partition variable=sigma2  dim=1
#pragma HLS array_partition variable=p1  dim=1
#pragma HLS array_partition variable=p2  dim=1
#pragma HLS array_partition variable=cnt1 dim=1
#pragma HLS array_partition variable=cnt2 dim=1
///#pragma HLS array_partition variable=stack1[0].storage  dim=1
    for(int j=0; j<UNROLL_FACTOR; j++)
    {
#pragma HLS UNROLL factor=UNROLL_FACTOR_tmp
        BFS(numVert, sigma1[j], stack1[j], _offset[j], _column[j], j, p1[j], cnt1[j]);
    }
  int bound = ceil(numVert/float(UNROLL_FACTOR));

  each_source: for (int i = 0; i < bound-1; i++)
  { 
    for(int j=0; j<UNROLL_FACTOR; j++)
    {
#pragma HLS LOOP_FLATTEN
    Wrapper(i, j, numVert, sigma1[j], sigma2[j], stack1[j], stack2[j], _offset[j], _column[j], p1[j], p2[j], cnt1[j], cnt2[j], delta1[j], delta2[j], _btwn[j]);
    }
  }

for(int j=0; j<UNROLL_FACTOR; j++)
{
#pragma HLS UNROLL factor=UNROLL_FACTOR_tmp
  if((bound-1)*UNROLL_FACTOR+j < numVert)
  {
    if(bound % 2 == 0)
      BACK(numVert, delta2[j], sigma2[j], stack2[j], _btwn[j], (bound-1)*UNROLL_FACTOR+j, p2[j], cnt2[j]);
    else
      BACK(numVert, delta1[j], sigma1[j], stack1[j], _btwn[j], (bound-1)*UNROLL_FACTOR+j, p1[j], cnt1[j]);
  }

}
  write_btwn: for (int i = 0; i < numVert; i++)
  {
#pragma HLS UNROLL factor=UNROLL_FACTOR_tmp
    for(int j = 1; j < UNROLL_FACTOR; j++)
    {
    #pragma HLS PIPELINE II=1
      _btwn[0][i] += _btwn[j][i];
    }
    btwn[i] = _btwn[0][i];
  }

}

double astar(long start_i, long start_j, long end_i, long end_j, double *grid,
              double *cost_so_far, long i_range, long j_range);
double astar_bw(long start_i, long start_j, long end_i, long end_j, double *grid,
              double *cost_so_far, long i_range, long j_range);
void get_neighbours( long i, long j, long i_size, long j_size,  long neighbours[8]);

#include "astar.h"
#include "linkedlist.h"
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <time.h>
#include <string.h>

double astar(long start_i, long start_j, long end_i, long end_j, double *grid,
  double *cost_so_far, long i_range, long j_range) {
    double result = 0.0 ;
    bool (*visited)[j_range] = malloc(sizeof(bool[i_range][j_range]));
    memset(visited, 0, sizeof(visited[0][0]) * i_range * j_range);

    /* time to check loop lenght */
    clock_t start, end;
    double cpu_time_used;

    start = clock();


    NODE list_head;
    NODE * current_list;
    list_head.next = 0;
    list_head.data = 0;
    long list_count = 0; 
    long neighbours[8];


    InsertOrdered(&list_head, 0, start_i, start_j);
    list_count++;
    //cost_so_far[start_i][start_j] = 0;
    cost_so_far[start_i*i_range+start_j] = 0;

    long current_i, current_j;
    visited[start_i][start_j] = 1;
    while (list_count>0) {

       current_list = list_head.next;
       list_count--;
       Delete_second(&list_head);

       current_i = current_list->i;
       current_j = current_list->j;

       if (current_i == end_i &&  current_j == end_j) break;

       end = clock(); /*check for how long have we run a search */
       cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
       if (cpu_time_used > 5) {
         result = 1;
         break;}
 
       get_neighbours(current_i, current_j, i_range, j_range, neighbours);

       int k;
       double new_cost;
       long next_i, next_j;
       double priority, heuristic;

       for (k=0; k<4; k++) {
         next_i = neighbours[k*2];
         next_j = neighbours[k*2+1];

         if (next_i <0 || next_j <0) continue;
         
         new_cost = cost_so_far[current_i*i_range+current_j] + grid[next_i*i_range + next_j] + 1;
         //if (new_cost < cost_so_far[next_i][next_j] || visited[next_i][next_j]==0) {
         if (new_cost < cost_so_far[next_i*i_range+next_j] || visited[next_i][next_j]==0) {
            //cost_so_far[next_i][next_j] = new_cost;
            cost_so_far[next_i*i_range+next_j] = new_cost;
            visited[next_i][next_j] = 1;
            //heuristic =  (end_i-next_i)^2 + (end_j - next_j)^2;
            heuristic =  abs(end_i-next_i) + abs(end_j - next_j);
            priority = new_cost + heuristic;
            InsertOrdered(&list_head, priority, next_i, next_j);
            list_count++;
         }
       }
    }
    free(visited);
    return result;
}

void get_neighbours( long i, long j, long i_size, long j_size, long neighbours[8]) {
    neighbours[0] = -1;
    neighbours[1] = -1;
    neighbours[2] = -1;
    neighbours[3] = -1;
    neighbours[4] = -1;
    neighbours[5] = -1;
    neighbours[6] = -1;
    neighbours[7] = -1;
    
    if (i>0) {
       neighbours[0] = i-1;
       neighbours[1] = j;}
    if (i<(i_size-1)) {
       neighbours[2] = i+1;
       neighbours[3] = j;}
    if (j>0) {
       neighbours[4] = i;
       neighbours[5] = j-1;}
    if (j<(j_size-1)) {
       neighbours[6] = i;
       neighbours[7] = j+1;}
    return;
}

double astar_bw(long start_i, long start_j, long end_i, long end_j, double *grid,
  double *cost_so_far, long i_range, long j_range) {
    double result = 0.0 ;
    bool (*visited)[j_range] = malloc(sizeof(bool[i_range][j_range]));
    memset(visited, 0, sizeof(visited[0][0]) * i_range * j_range);

    /* time to check loop lenght */
    clock_t start, end;
    double cpu_time_used;

    start = clock();


    NODE list_head;
    NODE * current_list;
    list_head.next = 0;
    list_head.data = 0;
    long list_count = 0; 
    long neighbours[8];


    InsertOrdered(&list_head, 0, start_i, start_j);
    list_count++;
    //cost_so_far[start_i][start_j] = 0;
    cost_so_far[start_i*i_range+start_j] = 0;

    long current_i, current_j;
    visited[start_i][start_j] = 1;
    while (list_count>0) {

       current_list = list_head.next;
       list_count--;
       Delete_second(&list_head);

       current_i = current_list->i;
       current_j = current_list->j;

       if (current_i == end_i &&  current_j == end_j) break;

       end = clock(); /*check for how long have we run a search */
       cpu_time_used = ((double) (end - start)) / CLOCKS_PER_SEC;
       if (cpu_time_used > 3) {
         result = 1;
         break;}
 
       get_neighbours(current_i, current_j, i_range, j_range, neighbours);

       int k;
       double new_cost;
       long next_i, next_j;
       double priority, heuristic;

       for (k=0; k<4; k++) {
         next_i = neighbours[k*2];
         next_j = neighbours[k*2+1];

         if (next_i <0 || next_j <0) continue;
         
         
         new_cost = cost_so_far[current_i*i_range+current_j] +  1;
         //if (new_cost < cost_so_far[next_i][next_j] || visited[next_i][next_j]==0) {
         if (new_cost < cost_so_far[next_i*i_range+next_j] || visited[next_i][next_j]==0
                                                     || grid[next_i*i_range + next_j]==0) {
            //cost_so_far[next_i][next_j] = new_cost;
            cost_so_far[next_i*i_range+next_j] = new_cost;
            visited[next_i][next_j] = 1;
            //heuristic =  (end_i-next_i)^2 + (end_j - next_j)^2;
            heuristic =  abs(end_i-next_i) + abs(end_j - next_j);
            priority = new_cost + heuristic;
            InsertOrdered(&list_head, priority, next_i, next_j);
            list_count++;
         }
       }
    }
    free(visited);
    return result;
}


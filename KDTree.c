#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>


double data_mean_dim(
    double data[],
    int start_index,
    int data_len,
    int chosen_dim,
    int dim)
    {
        int i;
        double array_sum =0;
        double array_average;
        for (i = start_index + chosen_dim; i< start_index + data_len; i+=dim)
        {
            array_sum += data[i];
        }
        array_average = array_sum/(data_len/dim);
        return array_average;
    }


double data_max_dim(
    double data[],
    int start_index,
    int data_len,
    int chosen_dim,
    int dim)
    {
        int i;
        double array_sum =0;
        double array_max =-9999;
        for (i = chosen_dim; i<data_len; i+=dim)
        {
            if (array_max < data[i+start_index])
            {
                array_max = data[i+start_index];
            }
        }
        return array_max;
    }

double data_min_dim(
    double data[],
    int start_index,
    int data_len,
    int chosen_dim,
    int dim)
    {
        int i;
        double array_sum =0;
        double array_min = 9999;
        for (i = chosen_dim; i<data_len; i+=dim)
        {
            if (array_min > data[i + start_index])
            {
                array_min = data[i + start_index];
            }
        }
        return array_min;
    }

void cluster_assign_func(
    double data[],
    int start_index,
    int data_len,
    short* cluster_assign,
    int cluster_assign_start_index,
    int cluster_assign_len,
    int chosen_dim,
    int dim,
    double mean)
    {
        int j=0;
        int i;
        // printf("cluster_val: %d\n", cluster_k);
        for(i = start_index+chosen_dim; i< start_index+data_len; i+=dim)
        {
            if (*(data+i) <= mean)
            {
                *(cluster_assign + cluster_assign_start_index) = 0;
            // printf("cluster assign 0: %hd, %f, %f\n",*(cluster_assign+j), *(data+i), mean);
                cluster_assign_start_index += 1;
            }
            else
            {
                *(cluster_assign+cluster_assign_start_index) = 0 +1 ;
                // printf("cluster assign 1: %hd, %f, %f\n",*(cluster_assign+j), *(data+i), mean);
                cluster_assign_start_index += 1;
                }
        }
            // cluster_k += 2;
        
        
        }


void insertion_sort_for_cluster_assign_and_data(
    double* data,
    int start_index,
    int data_len,
    short* cluster_assign,
    int cluster_assign_start_index,
    int cluster_assign_len,
    int dim)
    {
        int i;
        int temp_index;
        double temp_val;
        double temp_val_data[dim];


        for(i=cluster_assign_start_index + 1; i < cluster_assign_start_index + cluster_assign_len; i +=1)
        {
            temp_index = i;
            temp_val = cluster_assign[temp_index];
            int j, k;
            k=0;
            for(j=temp_index*dim; j< temp_index*dim + dim; j+=1)
            {
                temp_val_data[k] = data[j];
                k+=1;
            }
            while((temp_val < cluster_assign[temp_index - 1]) && (temp_index > cluster_assign_start_index))
            {
                cluster_assign[temp_index] = cluster_assign[temp_index - 1];
                for(j= temp_index*dim; j< temp_index*dim + dim; j+=1)
                {
                    data[j] = data[j-dim];
                }

                temp_index -=1;
            }
            cluster_assign[temp_index] = temp_val;
            k =0;
                for(j= temp_index*dim; j< temp_index*dim + dim; j+=1)
                {
                    data[j] = temp_val_data[k];
                    k+=1;
                }
        }


    }


void cluster_start_size_init(
    short* cluster_assign,
    int cluster_assign_start_index,
    int cluster_assign_len,
    int dim,
    int* cluster_start,
    int* cluster_size)
    {
        int i;
        cluster_start[0] = dim*cluster_assign_start_index;
        for(i=cluster_assign_start_index+1; i< cluster_assign_start_index+cluster_assign_len; i+=1)
        {
            if(cluster_assign[i]==cluster_assign[i-1]+1)
            {
                cluster_start[0+1] = i*dim;
                cluster_size[0] = (i-cluster_assign_start_index)*dim;
                // printf("i val %d, cluster_assign_len %d, cluster_assign_start_ind %d \n\n", i, cluster_assign_len, cluster_assign_start_index);
                cluster_size[0+1] = (cluster_assign_start_index + cluster_assign_len - i) *dim;
            }
        }

    }


double* cluster_boundary(
    double* data,
    int start_index,
    int data_len,
    int dim)
    {
        int chosen_dim;
        int j=0;

        double* tempdata1;
            tempdata1 = (double*)malloc(sizeof(double) * dim*2);


        for (chosen_dim =0; chosen_dim < dim; chosen_dim +=1)
        {
            tempdata1[j] = data_min_dim(data, start_index, data_len, chosen_dim, dim);
            tempdata1[j+1] = data_max_dim(data, start_index, data_len, chosen_dim, dim);
            // printf("startind, datalen %d, %d, %f, %f\n", start_index, data_len, tempdata1[j], tempdata1[j+1]);

            j+=2;

        }
        // printf("\n");


        return tempdata1;
    }

double* cluster_centroid_calc(
    double* cluster_min_max,
    int dim)
    {
        double* tempCentroid;
        int i, j;
        tempCentroid = (double*)malloc(sizeof(double)*dim);

        
        j=0;
        for (i = 0; i<dim+1; i+=2)
        {
            tempCentroid[j] = data_mean_dim(cluster_min_max, i, 2, 0, 1);
            // printf("hero %f\n", *(tempCentroid + j));
            j +=1;
        }
        return tempCentroid;
    }

int max_variation_of_dim(
    double* data,
    int start_ind,
    int data_len,
    int dim)
    {
        int i, j;
        double mean[dim];
        double var[dim];
        double max_var;
        int max_var_dim;
        max_var = 0;

        for (i=0; i<dim; i+=1)
        {

        mean[i] = data_mean_dim(data, start_ind, data_len, i, dim);
        printf("max_var_mean: %f\n", mean[i]);
        }

        for (j=0; j<dim; j+=1)
        {
            for (i = 0; i<data_len; i+=dim)
            {
                var[j] += pow((*(data+i) - mean[j]), 2);
            }
            printf("var[j] %f\n",var[j]);
            if (max_var < var[j])
            {
                max_var_dim = j;
            }
        }
            printf("max_var_dim %d\n",max_var_dim);
        return max_var_dim;
    }

double dim_distance(
        double query_pt,
        double boundary_low,
        double boundary_high)
    {
        double result;
        printf("QQ %f, %f, %f\n", query_pt, boundary_low, boundary_high);
        // printf("QQ %f, %f\n", fabs(query_pt-boundary_low) ,fabs(query_pt-boundary_high));
        if (boundary_low < query_pt && query_pt <= boundary_high)
        {
            printf("HELLO gonna pass result, result %f\n", 0.000);
            return 0.0;
        }
        else
        {
            if(fabs(query_pt-boundary_low) <= fabs(query_pt-boundary_high))
            {
                
                result = fabs(query_pt-boundary_low);
                // printf("HELLsds%f", result);
                printf("HELLO gonna pass result, result %f\n", fabs(query_pt-boundary_low));
                return result;
            }
            else
            {
            // printf("HELLO");
                result = fabs(query_pt-boundary_high);
                printf("HELLO gonna pass result, result %f\n", fabs(query_pt-boundary_high));
                return result;
            }
        }
    }

double* get_closest_pt_in_cluster(
    double* query_pt,
    double* data,
    int* cluster_start_for_kk,
    int* cluster_size_for_kk,
    double* result_pt,
    int closest_cluster,
    int dim)
{
    int data_index_to_start = cluster_start_for_kk[closest_cluster];
    int data_size_to_explore = cluster_size_for_kk[closest_cluster];

    printf( "TEST TEST%f, %d\n\n", data[8], closest_cluster);


    double sum_for_point, shortest_pt_dist;
    double sum_for_point_total[data_size_to_explore];
    int i, j,l, k=0;
    for (j=0; j<data_size_to_explore; j+=dim)
    {
        for( i= data_index_to_start; i <data_index_to_start + dim ; i++)
        {
            sum_for_point += pow(fabs(query_pt[i - data_index_to_start] -data[i+j]), 2);
            printf("%f, %f, %f, %f, querypt and datapt\n", query_pt[i - data_index_to_start], data[i+j], query_pt[i - data_index_to_start]-data[i+j], sum_for_point);
        }
        sum_for_point_total[k] =pow(sum_for_point, 0.5);

        k+=1;
        sum_for_point = 0;

    }
    shortest_pt_dist = sum_for_point_total[0];
    l = 0;
    for (k=0; k<data_size_to_explore/dim; k++)
    {
        printf("distance of query from points within chosen cluster %f, %d\n", sum_for_point_total[k], k);
        if ((shortest_pt_dist > sum_for_point_total[k+1])  && k+1< data_size_to_explore/dim)
        {
            shortest_pt_dist = sum_for_point_total[k+1];
            l =k+1;
        }
        printf("SHORTEST distance of query from points within chosen cluster %f, %d\n", shortest_pt_dist, l);
        result_pt[0] = shortest_pt_dist;
        printf("HELLO");
        int m =1;
        for(i = data_index_to_start + l*dim; i< data_index_to_start + l*dim +dim; i++)
        {
            printf("dist: %f, data %f\n", result_pt[0], data[i]);
            result_pt[m] = data[i];
            m++;
        }

        return result_pt;

    }

}

int bipartition(
    int dim, // 4 bytes: Number of dimensions of each data point
    int i0, // 4 bytes: Start INDEX of array of data points
    int data_len, // 4 bytes: End Ind of array of data points
    double* data, // 8 byte Pointer to Data: Reads & Traverses 8bytes at a time; Size = dim*(im-i0)
    int chosen_dim, // 4 bytes: Dimension to bipartition on
    int cluster_start[2], // 8 bytes: INDEX of cluster 0 start in data, INDEX of cluster 1 start in data
    int cluster_size[2], // 8 bytes : Number of elements in cluster 0, Number of elements in cluster 1
    double* cluster_bdry[2], // 16 bytes: 1 pointer to array of [min_dim1, max_dim1, min_dim2, max_dim2] in cluster 0; 1 pointer to array of [min_dim1, max_dim1, min_dim2, max_dim2] in cluster 1
    double* cluster_centroid[2], // 16 bytes: 1 pointer to array of [center_of_dim1, center_of_dim2] in cluster 0; 1 pointer to array of [center_of_dim1, center_of_dim2] in cluster 1 
    short* cluster_assign) // 8 bytes pointer to array on which cluster each datum belongs to
    {
        // Get the data on the chosen_dim
        double a = 0;
        int i, j, k;
        double mean;
        int cluster_assign_len;
        int m = 0;
        short b = 0;
        double* cluster_bdry_int;
        int cluster_assign_start_index;
        cluster_assign_start_index = i0/dim;
        cluster_assign = &b;
        cluster_assign_len = data_len/dim;


        mean = data_mean_dim(data, i0, data_len, chosen_dim, dim);
        printf("\nmean of chosen_dim: %f\n", mean);


        cluster_assign_func(data, i0, data_len, cluster_assign, cluster_assign_start_index, cluster_assign_len, chosen_dim, dim, mean);

        insertion_sort_for_cluster_assign_and_data(data, i0, data_len, cluster_assign, cluster_assign_start_index, cluster_assign_len, dim);
        
        
        for(i=i0/dim; i< i0/dim +cluster_assign_len; i+=1)
        {
            printf("cluster_assign %d \n", cluster_assign[i]);
        }

        for(i=i0; i<i0+data_len; i+=1)
        {
            printf("data %f \n", data[i]);
        }


        cluster_start_size_init(cluster_assign, cluster_assign_start_index, cluster_assign_len, dim, cluster_start, cluster_size);


        return 0;

    }


int kdtree(
    int dim, // 4 bytes: Number of dimensions of each data point
    int ndata,  // 4 bytes: Number of data points
    double* data,  // 8 byte pointer -> 8 byte double: Number of dimensions of each data point
    int kk, // 4 bytes: Number of clusters
    int* cluster_start, // 8 bytes pointer -> 4 bytes int: Number of dimensions of each data point
    int *cluster_size, // 
    double **cluster_bdry,
    double **cluster_centroid,
    int* cluster_start_for_kk,
    int* cluster_size_for_kk,
    short *cluster_assign)
    {
        cluster_start[0] = 0;
        cluster_size[0] = ndata*dim;
        int chosen_dim;
        int i; 

        // cluster_start_for_kk  = (int*)malloc(1*sizeof(int));
        cluster_start_for_kk[0] = 0;

        // cluster_size_for_kk  = (int*)malloc(1*sizeof(int));
        *cluster_size_for_kk = ndata*dim;

        // double** cluster_bdry_1  = (double**)malloc(2*sizeof(double*));


        for(i=2; i <kk+1; i+=1)
        {
         chosen_dim = max_variation_of_dim(data, cluster_start_for_kk[0], cluster_size_for_kk[0], dim); // always 0
         bipartition(dim, cluster_start_for_kk[0], cluster_size_for_kk[0], data, chosen_dim, cluster_start, cluster_size, cluster_bdry , cluster_centroid, cluster_assign);
         
         cluster_start_for_kk = (int*)realloc(cluster_start_for_kk, i*sizeof(int));

         for(int j=1; j< i; j+=1)
        {
            *(cluster_start_for_kk + j-1) = *(cluster_start_for_kk+ j);
        }

         *(cluster_start_for_kk + i-2) = cluster_start[0];
         *(cluster_start_for_kk + i-1) = cluster_start[1];


        cluster_size_for_kk = (int*)realloc(cluster_size_for_kk, i*sizeof(int));

        for(int j=1; j< i; j+=1)
        {
            *(cluster_size_for_kk + j-1) = *(cluster_size_for_kk+ j);
        }
         *(cluster_size_for_kk + i-2) = cluster_size[0];
         *(cluster_size_for_kk + i-1) = cluster_size[1];

         printf("Cluster %d start_ind : %d and num of data_elements: %d, and num of ndata %d, mem add %d\n", i-1 ,*(cluster_start_for_kk + i-2), cluster_size_for_kk[i-2], cluster_size_for_kk[i-2]/dim, (&cluster_start_for_kk[i-2]));
         printf("Cluster %d start_ind : %d and num of data_elements: %d, and num of ndata %d, mem add %d\n\n", i, cluster_start_for_kk[i-1], cluster_size_for_kk[i-1], cluster_size_for_kk[i-1]/dim, (&cluster_start_for_kk[i-1]));
        }

        // double** cluster_bdry_1 = (double**)malloc(2*dim*(sizeof(double*)));
        
        for (int j = 0; j<kk; j+=1)
        {
            printf("cli %d, %d\n\n", cluster_start_for_kk[j], cluster_size_for_kk[j]);
            cluster_bdry[j] = cluster_boundary(data, cluster_start_for_kk[j], cluster_size_for_kk[j],dim);
            cluster_centroid[j] = cluster_centroid_calc(cluster_bdry[j], dim);

            for (i =0; i <dim*2; i++)
            {
                printf("cluster_bdry: %f, j val= %d\n", cluster_bdry[j][i], j);
            }
            for (i=0;i<dim;i++)
            {
                printf("cluster_centroid: %f, j val= %d\n", cluster_centroid[j][i], j);
            }

        }
        return 0;
    }






double* search_kdtree(
    int dim,
    int ndata,
    double *data,
    int kk,
    int *cluster_start,
    int *cluster_size,
    double **cluster_bdry,
    double *query_pt,
    double *result_pt,
    int* cluster_start_for_kk,
    int* cluster_size_for_kk,
    double* global_shortest)
    {
        // double** dim_dist = (double**)malloc(kk*dim*2*sizeof(double*));
        double dim_dist[kk][dim];
        double* euc_dist = (double*)malloc(kk*sizeof(double));
        double sum_for_euc = 0;
        double euclidean_dist[kk];
        double closest_cluster_dist;
        int k=0;
        int i, j, l, m;
        double* result_pt_tot[dim+1];
        // double test;
        // test = dim_distance(query_pt[0], cluster_bdry[2][0], cluster_bdry[2][1]);
        // printf("%f passed\n", test);
            l=0;
        for(i=0; i<2*dim; i+=2)
        {
            for(j=0; j<kk; j++)
            {
                printf("He %f, %f, %f, %d, %d\n",query_pt[k], cluster_bdry[j][i], cluster_bdry[j][i+1], j,l);
                dim_dist[j][l] = dim_distance(query_pt[k], cluster_bdry[j][i], cluster_bdry[j][i+1]);
                printf("OK PASSED\n");
            }
            l+=1;
            k+=1;
        }
        for (j=0; j<kk; j++)
        {
            sum_for_euc =0;
            for (i = 0; i<dim; i++)
            {
                sum_for_euc += pow(dim_dist[j][i], 2);
            }
        euclidean_dist[j] = pow(sum_for_euc, 0.5);
        }

        // if get_closest_point_in_cluster
        closest_cluster_dist = euclidean_dist[0];
        m=0;
        for (j=0; j<kk; j++)
        {
            printf("distance of query point to cluster %d, is %f\n", j, euclidean_dist[j]);
            if (closest_cluster_dist > euclidean_dist[j+1] && j+1 <kk) 
            {
                closest_cluster_dist = euclidean_dist[j+1];
                m = j+1;
            }
        }
                printf("closest distance of query point to cluster %d, is %f\n", m, closest_cluster_dist);

        result_pt_tot[0] = get_closest_pt_in_cluster(query_pt, data, cluster_start_for_kk, cluster_size_for_kk,result_pt, m, dim);



        // double global_shortest[dim+1];
        for (i=0; i<dim+1; i++)
        {
            global_shortest[i] = result_pt_tot[0][i];
        }
        printf("result_pt_tot %f\n", global_shortest[0]);

        printf("global_shortest %f, %f", global_shortest[1], global_shortest[2]);


        int n=1;
        for (j=0; j<kk; j++)
        {
        if (euclidean_dist[j] < global_shortest[0] && j!=m)
        {
            printf("other cluster %d, and its distance is %f\n", j, euclidean_dist[j]);

            result_pt_tot[n] = get_closest_pt_in_cluster(query_pt, data, cluster_start_for_kk, cluster_size_for_kk,result_pt, j, dim);
            
            if (result_pt_tot[n][0]<global_shortest[0])
            {
                        for (i=0; i<dim+1; i++)
                        {
                            global_shortest[i] = result_pt_tot[n][i];
                        }
            }
            n+=1;
        }
        }
        printf("\nDONEEEE\n");
        return global_shortest;
    }




        
        // int i, j, l;
        // i=0;
        // l=0;
        // for (j=0; j<kk; j+=1)
        // {
        //     printf("%f, %f, %f, j=%d: check\n",cluster_bdry[j][l] , query_pt[i] , cluster_bdry[j][l+1], j);
        //     if (cluster_bdry[j][l] < query_pt[i] && query_pt[i] < cluster_bdry[j][l+1] && i<dim)
        //     {

        //     }
        //     while (cluster_bdry[j][l] < query_pt[i] && query_pt[i] < cluster_bdry[j][l+1] && i<dim)
        //  {

        //     printf("pass on dim %d for cluster %d\n", i, j);
        //     if (i == (dim -1))
        //     {
        //         printf("PASSED FULLY ON CLUSTER %d HOORAHHHH \n", j);
        //     }


        //     i+=1;
        //     l+=2;
        //  }
        //  i=0;
        //  l=0;
        // }





        // if (cluster_bdry[j][i] < query_pt[i] < cluster_bdry[j][i+1])
        // {
        //     if (cluster_bdry[j][i+2] < query_pt[i+1] < cluster_bdry[j][i+3])
        //     {
        //         if (cluster_bdry[j][i+4] < query_pt[i+2] < cluster_bdry[j][i+5])
        //         {
        //             if (cluster_bdry[j][i+4] < query_pt[dim] < cluster_bdry[j][i+5])
        //         }
        //     }
        // }

int main()
{
    int d, e;
    // d = kD_Tree(2, 6, data1, 2, );
    int cluster_start[2]; // 8 bytes: INDEX of cluster 0 start in data, INDEX of cluster 1 start in data
    int cluster_size[2]; // 8 bytes : Number of elements in cluster 0, Number of elements in cluster 1
    short* cluster_assign; // 8 bytes pointer to array on which cluster each datum belongs to
    int* cluster_start_for_kk = (int*)malloc(1*sizeof(int));;
    int* cluster_size_for_kk = (int*)malloc(1*sizeof(int));
    int dim = 2;
    double** cluster_bdry= (double**)malloc(2*dim*(sizeof(double*)));
    double** cluster_centroid = (double**)malloc(dim*(sizeof(double*)));
    int ndata = 6;
    double data1[] = {2, 3, 5, 4, 9, 6, 4, 7, 8, 1, 7, 2}; //[(2,3), (5,4), (9,6), (4,7), (8,1), (7,2)] -> [(2,3), (8,1), (7,2), ||(5,4), (9,6), (4,7)]  ->  [(8,1), (7,2) || (2,3)] || [(5,4) || (9,6), (4,7) ]
    int kk =4;
    int data_len = ndata*dim;
    double result_pt[dim+1];
    double global_shortest[dim+1];
    double* global_shortest_final;



    d = kdtree(dim, ndata, data1, kk, cluster_start, cluster_size, cluster_bdry , cluster_centroid, cluster_start_for_kk, cluster_size_for_kk, cluster_assign); // 8 bytes pointer to array on which cluster each datum belongs to
 
    for (int j = 0; j<kk; j+=1)
    {
        printf("cli main %d, %d\n\n", cluster_start_for_kk[j], cluster_size_for_kk[j]);
    }

    
    double query_pt[] = {7.8, 1.2};


    global_shortest_final = search_kdtree(dim, ndata, data1, kk, cluster_start, cluster_size, cluster_bdry, query_pt, result_pt, cluster_start_for_kk, cluster_size_for_kk, global_shortest);
    
    for(int i=1; i<dim+1; i++)
    {
        printf("YESSSSS %f", global_shortest_final[i]);
    }

    return 0;
}






    
    // jj = 1 // jj -- the number of clusters the dataset is partitioned into.
    //     // At the beginning, jj=1, meaing only one cluster, the whole dataset.
    //     while (jj < kk ) { // there are jj clusters at this moment
    //     for (j=0; j<jj; j++) { // j loops through indices of all jj clusters
    //     find the dimension of largest variance for the j-th cluster
    //     call bipartition() to partition the j-th cluster into 2 clusters
    //     }
    //     jj = 2*jj
    //     }

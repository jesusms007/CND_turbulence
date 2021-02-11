///LAST MOD .. NOV 9, 2016

////x0 - point below x
////x1 - poin above x
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "../allvars.h"
#include "../proto.h"

int ind_below(double coord)
{
    double delta2; double rad; int zin;
    delta2 = (double) delta;
    rad = (double) (-cloud_radius);
    zin = (fabs( (rad - coord) / delta2) ) ;
    
    return zin ;
    
}

int which_file(int number_of_files, double r) ///RETURNS A RANDOM NUMBER FROM 0 TO NUMBER_OF_FILES
{
    int result;
    int h = (int)(r * number_of_files);
    
    if (number_of_files == 1)
    {
        result = 1;
    }
    if (h == number_of_files)
    {
        result = (number_of_files-1) ;
    }
    else
    {
        result = h;
    }
        
    return result;
}

double trilinear(double x, double y, double z, double grid[grid_size], int coord, double *vel_grid)
                 
{
    double c;
    
    if       (x < grid[0] || x > grid[grid_size-1])
    {
        c = 0.0;
    }
    else if  (y < grid[0] || y > grid[grid_size-1])
    {
        c = 0.0;
    }
    else if  (z < grid[0] || z > grid[grid_size-1] )
    {
        c = 0.0;
    }
    else
    {
    double delta2 = (double) delta;
    ///indices
    int index_x0 = ind_below( x );
    int index_x1 = index_x0 + 1;
    int index_y0 = ind_below( y );
    int index_y1 = index_y0 + 1;
    int index_z0 = ind_below( z );
    int index_z1 = index_z0 + 1;
        
   
        if (index_x0 > grid_size-1 || index_y0 > grid_size-1 || index_z0 > grid_size-1)
        {
            printf("My_error, check here\n");
        }
        if (index_x0 < 0 || index_y0 < 0 || index_z0 < 0)
        {
            printf("My_error, check here. Negative index\n");
        }

    ///values
    double x0 = grid[index_x0];
    double y0 = grid[index_y0];
    double z0 = grid[index_z0];
    
    double pt_x0_y0_z0 = vel_grid[coord + index_x0 * dim + index_y0 * dim * grid_size + index_z0 * dim * grid_size * grid_size];
        
    double pt_x1_y0_z0 = vel_grid[coord + index_x1 * dim + index_y0 * dim * grid_size + index_z0 * dim * grid_size * grid_size];
    
    double pt_x0_y1_z0 = vel_grid[coord + index_x0 * dim + index_y1 * dim * grid_size + index_z0 * dim * grid_size * grid_size];
        
    double pt_x0_y0_z1 = vel_grid[coord + index_x0 * dim + index_y0 * dim * grid_size + index_z1 * dim * grid_size * grid_size];
        
    double pt_x1_y1_z0 = vel_grid[coord + index_x1 * dim + index_y1 * dim * grid_size + index_z0 * dim * grid_size * grid_size];
        
    double pt_x1_y0_z1 = vel_grid[coord + index_x1 * dim + index_y0 * dim * grid_size + index_z1 * dim * grid_size * grid_size];
        
    double pt_x0_y1_z1 = vel_grid[coord + index_x0 * dim + index_y1 * dim * grid_size + index_z1 * dim * grid_size * grid_size];
        
    double pt_x1_y1_z1 = vel_grid[coord + index_x1 * dim + index_y1 * dim * grid_size + index_z1 * dim * grid_size * grid_size];
        
        
    ////trilinear algorithm
    double xd = (x - x0)/ delta2;
    double yd = (y - y0)/ delta2;
    double zd = (z - z0)/ delta2;
        /*
    if (x < x0 || y < y0 || z < z0)
        {
            printf("Error, check here\n");
            printf("x=%f, x0=%f, index = %d, index1 = %d\n", x, x0, index_x0, index_x1);
            printf("y=%f, y0=%f, index = %d, index1 = %d\n", y, y0, index_y0, index_y1);
            printf("z=%f, z0=%f, index = %d, index1 = %d\n", z, z0, index_z0, index_z1);
            exit(0);
        }
         */
    double one_minus_xd = (1.0 - xd);

    double c00 = (pt_x0_y0_z0 * one_minus_xd) + (pt_x1_y0_z0 * xd);
    double c10 = (pt_x0_y1_z0 * one_minus_xd) + (pt_x1_y1_z0 * xd);
    double c01 = (pt_x0_y0_z1 * one_minus_xd) + (pt_x1_y0_z1 * xd);
    double c11 = (pt_x0_y1_z1 * one_minus_xd) + (pt_x1_y1_z1 * xd);

    double one_minus_yd = 1.0 - yd;

    double c0 = (c00*one_minus_yd) + (c10*yd);
    double c1 = (c01*one_minus_yd) + (c11*yd);

    c = (c0*(1.0-zd)) + (c1*zd);
    
    }
    
    return c;
}


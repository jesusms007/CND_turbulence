
////LAST MOD.. NOV 23, 2016
#include "../allvars.h"
#include <stdlib.h>
#include "../proto.h"

int check_octree(double x, double y, double z)
{
    int octree = -1;
    if (x <= 0.0)
    {
        if (y <= 0.0)
        {
            if (z > 0.0)
            {
                octree = 0;
            }
            if (z <= 0.0)
            {
                octree = 4;
            }
        }
        if (y > 0.0)
        {
            if (z > 0.0)
            {
                octree = 2;
            }
            if (z <= 0.0)
            {
                octree = 6;
            }
        }
    }
    
    if (x > 0)
    {
        if (y <= 0.0)
        {
            if (z > 0.0)
            {
                octree = 1;
            }
            if (z <= 0.0)
            {
                octree = 5;
            }
            
        }
        
        if (y > 0.0)
        {
            if (z > 0.0)
            {
                octree = 3;
            }
            if (z <= 0.0)
            {
                octree = 7;
            }
        }
        
    }
 
    if (octree == -1)
    {   printf("My_error:FATAL ERROR. Negative octree?\n");
        exit(0);
    }
    return octree;
}

int check_quadrant(double x, double y, double z)
{
    int quad = 0;
    if (x <= 0.0)
    {
        if (y <= 0.0)
        {
            quad = 2;
        }
        if (y > 0.0)
        {
            quad = 0;
        }
    }

    if (x > 0)
    {
        if (y <= 0.0)
        {
            quad = 3;
        }
        
        if (y > 0.0)
        {
            quad = 1;
        }
 
    }
    return quad;
}







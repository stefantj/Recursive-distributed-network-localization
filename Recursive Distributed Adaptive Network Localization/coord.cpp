/*
 *  coord.cpp
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */


#include "coord.h"
#include <math.h>

float coord_mult(coord c1, coord c2)
{
	return (c1[X]*c2[X]+c1[Y]*c2[Y]);
}

void coord_diff(coord c1, coord c2, coord res)
{
	res[X] = c1[X] - c2[X];
	res[Y] = c1[Y] - c2[Y];
}

void coord_add(coord c1, coord c2, coord res)
{
	res[X] = c1[X] + c2[X];
	res[Y] = c1[Y] + c2[Y];
}

void coord_scale(float a, coord c1)
{
	c1[X] = a*c1[X];
	c1[Y] = a*c1[Y];
}

void coord_rotate(float angle, coord c1)
{
	//Get magnitude:
	//assumed to be unit vector
	
	//Get current angle:
	float theta = atan2(c1[Y],c1[X]);
	theta += angle;
	
	c1[X] = cos(theta);
	c1[Y] = sin(theta);
}
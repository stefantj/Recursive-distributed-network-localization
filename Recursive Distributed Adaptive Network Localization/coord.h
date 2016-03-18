/*
 *  coord.h
 *
 *  Created by Stefan Jorgensen for Stanford EE373A project (Winter 2016).
 *  based on earlier work for EE359.
 *
 *  Code released under the MIT License (MIT)
 *
 */

#ifndef _COORD_H__
#define _COORD_H__

typedef float coord[2]; 
#define X 0
#define Y 1

//vector to scalar:
float coord_mult(coord c1, coord c2);

//vector to vector:
void coord_diff(coord c1, coord c2, coord res);
void coord_add(coord c1, coord c2, coord res);

void coord_scale(float a, coord c1);
void coord_rotate(float angle, coord c1);

#endif //_COORD_H__

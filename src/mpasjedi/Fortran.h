/*
 * (C) Copyright 2017-2023 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#pragma once

namespace mpas {

// Geometry key type
typedef int F90geom;
// GeometryIterator key type
typedef int F90iter;
// Model key type
typedef int F90model;
// Variables key type
typedef int F90vars;
// Locations key type
typedef int F90locs;
// Goms key type
typedef int F90goms;
// State key type
typedef int F90state;
// Increment key type
typedef int F90inc;
// Trajectory key type
typedef int F90traj;
// Background error covariance key type
typedef int F90bmat;
// Observation vector key type
typedef int F90ovec;
// Obs operator key type
typedef int F90hop;
// Observation data base type
typedef int F90odb;
// GetValues key
typedef int F90getvalues;
// LinearGetValues key
typedef int F90lineargetvalues;

typedef double real_type;
/*
 * Above, double is consistent with OOPS interface.
 * typedef float real_type will not work here, however,
 * branch feature/sp_t1 gives an example on how this works.
 */

}  // namespace mpas

/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef MPAS_MODEL_RANDOM_F_H_
#define MPAS_MODEL_RANDOM_F_H_

namespace mpas {
extern "C" {
  void random_f(const int &, double *);
}
}  // namespace mpas

#endif  // MPAS_MODEL_RANDOM_F_H_

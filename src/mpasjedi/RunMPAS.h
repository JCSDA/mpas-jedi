/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MPASJEDI_RUNMPAS_H_
#define MPASJEDI_RUNMPAS_H_

#include <oops/runs/Run.h>

namespace mpas {

/*!
 *  RunMPAS encapsulates one MPAS/OOPS run.
 */

// -----------------------------------------------------------------------------

class RunMPAS : public oops::Run {
 public:
  RunMPAS(int, char **);
  ~RunMPAS();
};

// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPASJEDI_RUNMPAS_H_

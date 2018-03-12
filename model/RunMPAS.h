#ifndef MPAS_RUNMPAS_H_
#define MPAS_RUNMPAS_H_

#include "oops/runs/Run.h"

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
#endif  // MPAS_RUNMPAS_H_

#ifndef MPAS_RUNMPAS_H_
#define MPAS_RUNMPAS_H_

#include "oops/runs/Run.h"

namespace mpas {

/*!
 *  RunLfric encapsulates one LFRic/OOPS run.
 */

// -----------------------------------------------------------------------------

class RunLfric : public oops::Run {
 public:
  RunLfric(int, char **);
  ~RunLfric();
};

// -----------------------------------------------------------------------------

}  // namespace mpas
#endif  // MPAS_RUNMPAS_H_

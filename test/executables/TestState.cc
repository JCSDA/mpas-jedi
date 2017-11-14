
#include "model/Traits.h"
#include "oops/runs/Run.h"
#include "test/interface/State.h"

int main(int argc,  char ** argv) {
  oops::Run run(argc, argv);
  test::State<mpas::Traits> tests;
  run.execute(tests);
  return 0;
};


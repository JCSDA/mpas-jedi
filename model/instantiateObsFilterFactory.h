/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef MODEL_INSTANTIATEOBSFILTERFACTORY_H_
#define MODEL_INSTANTIATEOBSFILTERFACTORY_H_

#include "MPASTraits.h"
#include "oops/base/instantiateObsFilterFactory.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/ObsFilter.h"
#include "ufo/BackgroundCheck.h"
#include "ufo/BlackList.h"
#include "ufo/ObsBoundsCheck.h"
#include "ufo/ObsDomainCheck.h"
#include "ufo/ObsPreQC.h"

namespace mpas {

void instantiateObsFilterFactory() {
  oops::instantiateObsFilterFactory<MPASTraits>();
  static oops::FilterMaker<MPASTraits,
                           oops::ObsFilter<MPASTraits, ufo::ObsPreQC>
                          > makerChk1_("PreQC");
  static oops::FilterMaker<MPASTraits,
                 oops::ObsFilter<MPASTraits, ufo::ObsDomainCheck>
                          > makerChk2_("Domain Check");
  static oops::FilterMaker<MPASTraits,
                 oops::ObsFilter<MPASTraits, ufo::ObsBoundsCheck>
                          > makerChk3_("Bounds Check");
  static oops::FilterMaker<MPASTraits,
                 oops::ObsFilter<MPASTraits, ufo::BlackList>
                          > makerChk4_("BlackList");
  static oops::FilterMaker<MPASTraits,
                 oops::ObsFilter<MPASTraits, ufo::BackgroundCheck>
                          > makerChk5_("Background Check");
}

}  // namespace mpas

#endif  // MODEL_INSTANTIATEOBSFILTERFACTORY_H_

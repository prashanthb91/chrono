// =============================================================================
// PROJECT CHRONO - http://projectchrono.org
//
// Copyright (c) 2014 projectchrono.org
// All right reserved.
//
// Use of this source code is governed by a BSD-style license that can be found
// in the LICENSE file at the top level of the distribution and at
// http://projectchrono.org/license-chrono.txt.
//
// =============================================================================
// Authors: Radu Serban
// =============================================================================
//
// Simple driveline model for a tracked vehicle. This template splits the input
// torque to the left and right tracks using a simple analytical model of a
// Torsen limited-slip differential and the given driver steering input.
//
// =============================================================================

#ifndef CH_SIMPLE_TRACK_DRIVELINE_H
#define CH_SIMPLE_TRACK_DRIVELINE_H

#include "chrono_vehicle/ChApiVehicle.h"
#include "chrono_vehicle/tracked_vehicle/ChTrackDriveline.h"

namespace chrono {
namespace vehicle {

///
///
///
class CH_VEHICLE_API ChSimpleTrackDriveline : public ChTrackDriveline {
  public:
    ChSimpleTrackDriveline(const std::string& name);
    virtual ~ChSimpleTrackDriveline() {}

    /// Initialize the driveline subsystem.
    /// This function connects this driveline subsystem to the sprockets of the
    /// two track assembly subsystems.
    virtual void Initialize(ChSharedPtr<ChBody> chassis,              ///< handle to the chassis body
                            ChSharedPtr<ChTrackAssembly> track_left,  ///< handle to the left track assembly
                            ChSharedPtr<ChTrackAssembly> track_right  ///< handle to the right track assembly
                            ) override;

    /// Get the angular speed of the driveshaft.
    /// This represents the output from the driveline subsystem that is passed to
    /// the powertrain system.
    virtual double GetDriveshaftSpeed() const override;

    /// Update the driveline subsystem.
    /// The motor torque represents the input to the driveline subsystem from the
    /// powertrain system.
    virtual void Update(double steering, double torque) override;

    /// Get the motor torque to be applied to the specified sprocket.
    virtual double GetSprocketTorque(VehicleSide side) const override;

  protected:
    /// Return the torque bias ratio for the differential.
    /// This is a simple model of a Torsen limited-slip differential.
    virtual double GetDifferentialMaxBias() const = 0;

  private:
      ChSharedPtr<ChShaft> m_shaft_left;
      ChSharedPtr<ChShaft> m_shaft_right;
};

}  // end namespace vehicle
}  // end namespace chrono

#endif
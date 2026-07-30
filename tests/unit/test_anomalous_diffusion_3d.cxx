
#include "gtest/gtest.h"

#include "fake_mesh_fixture.hxx"
#include "test_extras.hxx" // FakeMesh                                                                                                                                                                                                                                            

#include "../../include/anomalous_diffusion_3d.hxx"

/// Global mesh                                                                                                                                                                                                                                                                   
namespace bout {
namespace globals {
extern Mesh* mesh;
} // namespace globals                                                                                                                                                                                                                                                            
} // namespace bout                                                                                                                                                                                                                                                               

// The unit tests use the global mesh                                                                                                                                                                                                                                             
using namespace bout::globals;

#include <bout/field_factory.hxx> // For generating functions                                                                                                                                                                                                                     

// Reuse the "standard" fixture for FakeMesh                                                                                                                                                                                                                                      
using AnomalousDiffusion3DTest = FakeMeshFixture;

TEST_F(AnomalousDiffusion3DTest, NoDiffusion) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;

  AnomalousDiffusion3D component("h", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);
  mesh->communicate(N);

  Options state;
  state["species"]["h"]["density"] = N;

  // If D is not set, then the diffusion should not be calculated                                                                                                                                                                                                                 
  component.transform(state);

  ASSERT_FALSE(state["species"]["h"].isSet("density_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("momentum_source"));
  ASSERT_FALSE(state["species"]["h"].isSet("energy_source"));
}

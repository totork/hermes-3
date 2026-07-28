#include "gtest/gtest.h"
#include "test_extras.hxx" // FakeMesh                                                                                                               
#include "../../include/kmodel.hxx"

/// Global mesh                                                                                                                                    
                                                                                                                                                   
namespace bout{
  namespace globals{
    extern Mesh *mesh;
  }
}                                                                                                                                                     
using namespace bout::globals;

#include <bout/field_factory.hxx>                                                                                                                                                        
using KmodelTest = FakeMeshFixture;

TEST_F(KmodelTest, CreateComponent) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;
  mesh->getCoordinates()->Bxy = 1.0;

  Kmodel component("test", options, nullptr);
}

/*
TEST_F(KmodelTest, NoCharge) {
  Options options;
  options["units"]["meters"] = 1.0;
  options["units"]["seconds"] = 1.0;
  options["units"]["Tesla"] = 1.0;

  mesh->getCoordinates()->Bxy = FieldFactory::get()->create2D("1 + x", &options, mesh);

  
  Kmodel component("h+", options, nullptr);

  Field3D N = FieldFactory::get()->create3D("1 + y * (x - 0.5)", &options, mesh);

  Field3D T = FieldFactory::get()->create3D("1 + 0.2 * x + sin(z)", &options, mesh);

  mesh->communicate(N,T);

  Options state;
  state["species"]["h+"]["density"] = N;
  state["species"]["h+"]["temperature"] = T;
  state["species"]["h+"]["pressure"] = N*T;
  state["species"]["h+"]["charge"] = 1;

}
*/

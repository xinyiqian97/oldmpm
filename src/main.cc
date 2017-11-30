#include <array>
#include <iostream>
#include <memory>
#include <mpi.h>

#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
#include <boost/serialization/array.hpp>
#include <boost/serialization/vector.hpp>
#define EIGEN_DENSEBASE_PLUGIN "EigenDenseBaseAddons.h"

#include "container.h"
#include "handler.h"
#include "node.h"
#include "particle.h"

#include "Eigen/Dense"

int main(int argc, char** argv) {
  unsigned long long id = 42;
  const unsigned Dim = 3;
  Eigen::Matrix<double, Dim, 1> coord;
  coord.setZero();
  coord[0] = 0.1;
  coord[1] = 0.2;
  coord[2] = 0.3;
  /*

  auto node = std::make_shared<mpm::Node<Dim>>(id, coord);
  std::cout << "Node id: " << node->id() << '\n';

  auto nodehandler = std::make_shared<mpm::Handler<mpm::Node<Dim>>>();
  nodehandler->insert(node);

  for (auto itr = nodehandler->begin(); itr != nodehandler->end(); ++itr)
    std::cout << ((*itr).second)->id() << '\n';

  auto nodecontainer = std::make_shared<mpm::Container<mpm::Node<Dim>>>();
  nodecontainer->insert(node);
  nodecontainer->insert(node);


  */

  MPI_Init(&argc, &argv);
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  std::stringstream ss;
  if (rank == 0) {
    mpm::Particle<Dim> p(id, coord);
    std::cout << "master particle id " << p.id() << std::endl;

    boost::archive::text_oarchive oa(ss);

    oa << p;
  std::cout << ss.str();
  }


  const int buf_size = 2000;
  char buf[buf_size];
  strncpy(buf, ss.str().c_str(), buf_size);

  MPI_Bcast(buf, buf_size, MPI_BYTE, 0, MPI_COMM_WORLD);

  if (rank != 1) {
    ss << buf;

    mpm::Particle<Dim> q;

    boost::archive::text_iarchive(ss) >> q;
    std::cout << rank << " particle " << q.id() 
      << " " << q.coordinates()[0]
      << " " << q.coordinates()[1]
      << " " << q.coordinates()[2] << std::endl;
  }

  MPI_Finalize();
}

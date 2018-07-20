#ifndef MPM_MESH_H_
#define MPM_MESH_H_

#include <iostream>
#include <limits>
#include <memory>
#include <vector>

#include "Eigen/Dense"
#include <tbb/parallel_for.h>
#include <tbb/parallel_for_each.h>

#include "cell.h"
#include "container.h"
#include "factory.h"
#include "node.h"
#include "particle.h"
#include "particle_base.h"

namespace mpm {

//! Mesh class
//! \brief Base class that stores the information about meshes
//! \details Mesh class which stores the particles, nodes, cells and neighbours
//! \tparam Tdim Dimension
template <unsigned Tdim>
class Mesh {
 public:
  //! Define a vector of size dimension
  using VectorDim = Eigen::Matrix<double, Tdim, 1>;

  // Construct a mesh with a global unique id
  //! \param[in] id Global mesh id
  Mesh(unsigned id);

  //! Default destructor
  ~Mesh() = default;

  //! Delete copy constructor
  Mesh(const Mesh<Tdim>&) = delete;

  //! Delete assignement operator
  Mesh& operator=(const Mesh<Tdim>&) = delete;

  //! Return id of the mesh
  unsigned id() const { return id_; }

  //! Create nodes from coordinates
  //! \param[in] gnid Global node id
  //! \param[in] node_type Node type
  //! \param[in] coordinates Nodal coordinates
  //! \retval status Create node status
  bool create_nodes(mpm::Index gnid, const std::string& node_type,
                    const std::vector<VectorDim>& coordinates);

  //! Add a node to the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool add_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Remove a node from the mesh
  //! \param[in] node A shared pointer to node
  //! \retval insertion_status Return the successful addition of a node
  bool remove_node(const std::shared_ptr<mpm::NodeBase<Tdim>>& node);

  //! Return the number of nodes
  mpm::Index nnodes() const { return nodes_.size(); }

  //! Iterate over nodes
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_nodes(Toper oper);

  //! Create cells from list of nodes
  //! \param[in] gcid Global cell id
  //! \param[in] shapefn Shape function of the cell
  //! \param[in] cells Node ids of cells
  //! \retval status Create cells status
  bool create_cells(mpm::Index gnid,
                    const std::shared_ptr<mpm::ShapeFn<Tdim>>& shapefn,
                    const std::vector<std::vector<mpm::Index>>& cells);

  //! Add a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool add_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Remove a cell from the mesh
  //! \param[in] cell A shared pointer to cell
  //! \retval insertion_status Return the successful addition of a cell
  bool remove_cell(const std::shared_ptr<mpm::Cell<Tdim>>& cell);

  //! Number of cells in the mesh
  mpm::Index ncells() const { return cells_.size(); }

  //! Iterate over cells
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_cells(Toper oper);

  //! Create particles from coordinates
  //! \param[in] gpid Global particle id
  //! \param[in] particle_type Particle type
  //! \param[in] coordinates Nodal coordinates
  //! \retval status Create particle status
  bool create_particles(mpm::Index gpid, const std::string& particle_type,
                        const std::vector<VectorDim>& coordinates);

  //! Add a particle to the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool add_particle(const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Remove a particle from the mesh
  //! \param[in] particle A shared pointer to particle
  //! \retval insertion_status Return the successful addition of a particle
  bool remove_particle(
      const std::shared_ptr<mpm::ParticleBase<Tdim>>& particle);

  //! Number of particles in the mesh
  mpm::Index nparticles() const { return particles_.size(); }

  //! Locate particles in a cell
  //! Iterate over all cells in a mesh to find the cell in which particles
  //! are located.
  //! \retval particles Particles which cannot be located in the mesh
  std::vector<std::shared_ptr<mpm::ParticleBase<Tdim>>> locate_particles_mesh();

  //! Iterate over particles
  //! \tparam Toper Callable object typically a baseclass functor
  template <typename Toper>
  void iterate_over_particles(Toper oper);

  //! Return status of the mesh. A mesh is active, if at least one particle is
  //! present
  bool status() const { return particles_.size(); }

  //! Add a neighbour mesh, using the local id for the new mesh and a mesh
  //! pointer
  //! \param[in] local_id local id of the mesh
  //! \param[in] neighbour A shared pointer to the neighbouring mesh
  //! \retval insertion_status Return the successful addition of a node
  bool add_neighbour(unsigned local_id,
                     const std::shared_ptr<Mesh<Tdim>>& neighbour);

  //! Return the number of neighbouring meshes
  unsigned nneighbours() const { return neighbour_meshes_.size(); }

  /*
  std::vector<uint64_t> pack_particles(Container<std::shared_ptr<mpm::ParticleBase<Tdim>>> &particles) {
    std::list<std::vector<uint64_t>> buf;
    for (auto &p : particles) {
      buf.push( p->pack(buf) );
      remove_particle(p);
    }

    unsigned size = 0;
    for (auto &b : buf) {
      size += b.size();
    }
    std::vector<uint64_t> actual_buf(size);
    auto abuf = actual_buf.data();

    for (auto &b : buf) {
      memcpy(abuf, b.data(), b.size());
      abuf += b.size();
    }

    return actual_buf;
  }
  void unpack_particles(std::vector<uint64_t> buf) {
    std::list<std::vector<uint64_t>> bufs;

    auto b = buf.data();
    while (b < buf.cend()) {
      int len = pack_size(b);
      std::vector<uint64_t> bb(len);
      memcpy(bb.data(), b, len);
      bufs.push( bb );
    }

    for (auto &b : bufs) {
      auto mytype = typenmap[b[0]];
      auto p = particle_factory(mytype);
      p.unpack(b);
      add_particle(p);
    }

  }

  std::map<uint64_t, std::string> typemap;

  void sendrecv_particles(std::vector<std::shared_ptr<> &particles) {
    for (p : particles) {
      auto b = p->pack();
      mpi_request req;
      MPI_Ibsend(b.data(), p.mesh.rank, &req);
      reqs.push(req);
    }

    mpi_recv(&npart, 1);
    for (unsigned i = 0; i < npart; i++) {

      MPI_recv(&sz, )

      MPI_recv(b, sz)
    }

  }

  void sendrecv_particles(std::vector<std::shared_ptr<> &particles) {
    buffer_type b2[89];
    for (n : neighbour_meshes_) {
      auto b = pack_particles(filter(p for p in particles if p.mesh == n));
      reqs.push( mpi_isend(b, rank=n.rank) );
      reqs.push( mpi_irecv(b2, rank=n.rank) );
    }
    mpi_wait(reqs.data());

    for (bb : b2) {
      unpack_particles(bb);
    }
  }
  */

 protected:
  //! mesh id
  unsigned id_{std::numeric_limits<unsigned>::max()};

  //! Container of mesh neighbours
  Handler<Mesh<Tdim>> neighbour_meshes_;

  //! Container of particles
  Container<ParticleBase<Tdim>> particles_;

  //! Container of nodes
  Container<NodeBase<Tdim>> nodes_;

  //! Map of nodes for fast retrieval
  Handler<NodeBase<Tdim>> map_nodes_;

  //! Container of cells
  Container<Cell<Tdim>> cells_;
};  // Mesh class
}  // namespace mpm

#include "mesh.tcc"

#endif  // MPM_MESH_H_

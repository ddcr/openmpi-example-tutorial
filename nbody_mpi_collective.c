/*
* @Author: Domingos Rodrigues <ddcr@lcc.ufmg.br>
* @Date:   2019-06-28 15:00:08
* @Last Modified by:   Domingos Rodrigues
* @Last Modified time: 2019-07-01 07:09:26
* File imported from:
*        https://www.fz-juelich.de/SharedDocs/Downloads/IAS/JSC/EN/
*        slides/mpi/course-materials-mpi-openmp.zip?__blob=publicationFile
*/

#define _XOPEN_SOURCE 600

#include <assert.h>
#include <math.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mpi.h>

#ifdef TEST_SCRATCH
#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include "xstring.h"
#endif
#include "xdmf_write.h"

typedef struct {
  double* x;
  double* y;
  double* z;
} vectors;

typedef struct {
  size_t count;
  size_t capacity;
  vectors x;
  vectors v;
  vectors a;
  double* q;
} particles;

#define block_allocate(block, n) \
  if (posix_memalign((void**)&block, 64, n * sizeof(double))) { \
    perror("Could not allocate memory for particles"); \
    exit(EXIT_FAILURE); \
  } \

particles particles_allocate(size_t n) {
  particles p;

  block_allocate(p.x.x, n);
  block_allocate(p.x.y, n);
  block_allocate(p.x.z, n);

  block_allocate(p.v.x, n);
  block_allocate(p.v.y, n);
  block_allocate(p.v.z, n);

  block_allocate(p.a.x, n);
  block_allocate(p.a.y, n);
  block_allocate(p.a.z, n);

  block_allocate(p.q, n);

  for (size_t i = 0; i < n; ++i) {
    p.x.x[i] = p.x.y[i] = p.x.z[i] = 0;
    p.v.x[i] = p.v.y[i] = p.v.z[i] = 0;
    p.a.x[i] = p.a.y[i] = p.a.z[i] = 0;
    p.q[i] = 0;
  }

  p.count = 0;
  p.capacity = n;

  return p;
}

void sources_copy(particles source, particles* dest) {
  if (dest->capacity < source.count) {
    fprintf(stderr, "Destination particle container too small to hold source particles.");
    exit(EXIT_FAILURE);
  }

  for (size_t i = 0; i < source.count; ++i) {
    dest->x.x[i] = source.x.x[i]; dest->x.y[i] = source.x.y[i]; dest->x.z[i] = source.x.z[i];
    dest->q[i] = source.q[i];
  }

  dest->count = source.count;
}

#define block_bcast(block) MPI_Bcast(block, count, MPI_DOUBLE, root, MPI_COMM_WORLD);

void sources_bcast(particles* sources, int root) {
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);

  uint64_t count = sources->count;
  MPI_Bcast(&count, 1, MPI_UINT64_T, root, MPI_COMM_WORLD);

  block_bcast(sources->x.x);
  block_bcast(sources->x.y);
  block_bcast(sources->x.z);
  block_bcast(sources->q);

  sources->count = count;
}

#define block_read(block) MPI_File_read_ordered(f, block, nl, MPI_DOUBLE, MPI_STATUS_IGNORE);

particles particles_from_file(const char* filename) {
  MPI_File f;
  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_RDONLY, MPI_INFO_NULL, &f);
  MPI_File_set_errhandler(f, MPI_ERRORS_ARE_FATAL);

  uint64_t n_;
  MPI_File_read_all(f, &n_, 1, MPI_UINT64_T, MPI_STATUS_IGNORE);
  size_t n = n_;
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);
  size_t nl = n / s + ((r < (n % s)) ? 1 : 0);
  MPI_File_set_view(f, sizeof(uint64_t), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

  particles p = particles_allocate(nl);

  block_read(p.x.x);
  block_read(p.x.y);
  block_read(p.x.z);

  block_read(p.q);

  block_read(p.v.x);
  block_read(p.v.y);
  block_read(p.v.z);

  MPI_File_close(&f);

  p.count = nl;

  return p;
}

#define block_write(block) MPI_File_write_ordered(f, block, p.count, MPI_DOUBLE, MPI_STATUS_IGNORE);

#ifdef TEST_SCRATCH
#define block_write_proc(block) write(fd_proc, block, p.count*sizeof(double));

void particles_to_file_per_proc(const char* filename_p, particles p) {
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);
  int fd_proc = open(filename_p, O_WRONLY | O_CREAT | O_APPEND, 0664);

  if (!fd_proc){
    printf("Could not write to %s\n", filename_p);
    close(fd_proc);
    (void) unlink(filename_p);
  } else {
    printf("Let us write %ld particles to %s\n", p.count, filename_p);
    write(fd_proc, &p.count, sizeof(int));

    block_write_proc(p.x.x);
    block_write_proc(p.x.y);
    block_write_proc(p.x.z);

    block_write_proc(p.q);

    block_write_proc(p.v.x);
    block_write_proc(p.v.y);
    block_write_proc(p.v.z);

    close(fd_proc);
  }
}
#endif

void particles_to_file(const char* filename, particles p) {
  MPI_File f;

  MPI_File_open(MPI_COMM_WORLD, filename, MPI_MODE_WRONLY | MPI_MODE_CREATE, MPI_INFO_NULL, &f);
  MPI_File_set_errhandler(f, MPI_ERRORS_ARE_FATAL);

  uint64_t n;
  MPI_Reduce(&p.count, &n, 1, MPI_UINT64_T, MPI_SUM, 0, MPI_COMM_WORLD);

  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);

  if (r == 0)
    MPI_File_write(f, &n, 1, MPI_UINT64_T, MPI_STATUS_IGNORE);

  MPI_File_set_view(f, sizeof(uint64_t), MPI_DOUBLE, MPI_DOUBLE, "native", MPI_INFO_NULL);

  block_write(p.x.x);
  block_write(p.x.y);
  block_write(p.x.z);

  block_write(p.q);

  block_write(p.v.x);
  block_write(p.v.y);
  block_write(p.v.z);

  MPI_File_close(&f);
}

void particles_free(particles* p) {
  free(p->x.x); free(p->x.y); free(p->x.z);
  p->x.x = p->x.y = p->x.z = NULL;
  free(p->q);
  p->q = NULL;
  free(p->v.x); free(p->v.y); free(p->v.z);
  p->v.x = p->v.y = p->v.z = NULL;
  free(p->a.x); free(p->a.y); free(p->a.z);
  p->a.x = p->a.y = p->a.z = NULL;

  p->count = 0;
  p->capacity = 0;
}

void accelerations_calculate_self(
    size_t np,
    double ax[static restrict np],
    double ay[static restrict np],
    double az[static restrict np],
    const double x[static restrict np],
    const double y[static restrict np],
    const double z[static restrict np],
    const double q[static restrict np]
) {
  for (size_t i = 0; i < np; ++i) {
    double accx = 0.0;
    double accy = 0.0;
    double accz = 0.0;

    for (size_t j = 0; j < i; ++j) {
      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double dz = z[i] - z[j];
      double r = sqrt(dx * dx + dy * dy + dz * dz);
      double acc = q[i] * q[j] / (r * r * r);

      accx += acc * dx; accy += acc * dy; accz += acc * dz;
    }

    for (size_t j = i + 1; j < np; ++j) {
      double dx = x[i] - x[j];
      double dy = y[i] - y[j];
      double dz = z[i] - z[j];
      double r = sqrt(dx * dx + dy * dy + dz * dz);
      double acc = q[i] * q[j] / (r * r * r);

      accx += acc * dx; accy += acc * dy; accz += acc * dz;
    }

    ax[i] += accx; ay[i] += accy; az[i] += accz;
  }
}

void accelerations_calculate_other(
    size_t nt,
    double ax[static restrict nt],
    double ay[static restrict nt],
    double az[static restrict nt],
    const double xt[static restrict nt],
    const double yt[static restrict nt],
    const double zt[static restrict nt],
    const double qt[static restrict nt],
    size_t ns,
    const double xs[static restrict ns],
    const double ys[static restrict ns],
    const double zs[static restrict ns],
    const double qs[static restrict ns]
) {
  for (size_t i = 0; i < nt; ++i) {
    double accx = 0.0;
    double accy = 0.0;
    double accz = 0.0;

    for (size_t j = 0; j < ns; ++j) {
      double dx = xt[i] - xs[j];
      double dy = yt[i] - ys[j];
      double dz = zt[i] - zs[j];
      double r = sqrt(dx * dx + dy * dy + dz * dz);
      double acc = qt[i] * qs[j] / (r * r * r);

      accx += acc * dx; accy += acc * dy; accz += acc * dz;
    }

    ax[i] += accx; ay[i] += accy; az[i] += accz;
  }
}

void accelerations_calculate(particles targets) {
  for (size_t i = 0; i < targets.count; ++i)
    targets.a.x[i] = targets.a.y[i] = targets.a.z[i] = 0.0;

  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);

  particles sources = particles_allocate(targets.capacity + 1);

  for (int i = 0; i < s; ++ i) {
    if (i == r) {
      sources_copy(targets, &sources);
      sources_bcast(&sources, i);
      accelerations_calculate_self(targets.count, targets.a.x, targets.a.y, targets.a.z, targets.x.x, targets.x.y, targets.x.z, targets.q);
   } else {
      sources_bcast(&sources, i);
      accelerations_calculate_other(
          targets.count, targets.a.x, targets.a.y, targets.a.z, targets.x.x, targets.x.y, targets.x.z, targets.q,
          sources.count, sources.x.x, sources.x.y, sources.x.z, sources.q
      );
    }
  }

  particles_free(&sources);
}

void step(size_t np, double dt, double u[static restrict np], const double u_t[static restrict np]) {
  for (size_t ip = 0; ip < np; ++ip)
    u[ip] += dt * u_t[ip];
}

void velocity_step(particles p, double dt) {
  step(p.count, dt, p.v.x, p.a.x);
  step(p.count, dt, p.v.y, p.a.y);
  step(p.count, dt, p.v.z, p.a.z);
}

void position_step(particles p, double dt) {
  step(p.count, dt, p.x.x, p.v.x);
  step(p.count, dt, p.x.y, p.v.y);
  step(p.count, dt, p.x.z, p.v.z);
}

int main(int argc, char* argv[argc + 1]) {
  if (argc < 2) {
    fprintf(stderr, "Usage: nbody <input file>\n");
    exit(EXIT_FAILURE);
  }

  MPI_Init(&argc, &argv);
  int r; MPI_Comm_rank(MPI_COMM_WORLD, &r);
  int s; MPI_Comm_size(MPI_COMM_WORLD, &s);

  particles p = particles_from_file(argv[1]);

  accelerations_calculate(p);
  const size_t nt = 50;
  const double dt = 0.1;
  char filename[18];
#ifdef TEST_SCRATCH
  char *filename_rank = NULL;
#endif

  for (size_t it = 0; it < nt; ++it) {
    // Write current state to disk
    sprintf(filename, "particles_%03zu.bin", it);
    particles_to_file(filename, p);

#ifdef TEST_SCRATCH
    filename_rank = xstrdup_printf("particles_%03zu.bin.proc%02zu", it, r);
    particles_to_file_per_proc(filename_rank, p);
#endif

    if (r == 0) { fprintf(stderr, "Working on step %zu... ", it + 1); fflush(stderr); }
    // Velocity half step
    velocity_step(p, 0.5 * dt);

    // Position full step
    position_step(p, dt);

    // Calculate accelerations at new positions
    accelerations_calculate(p);

    // Velocity half step
    velocity_step(p, 0.5 * dt);
    if (r == 0) fprintf(stderr, "done.\n");
  }

  // Write final state to disk
  sprintf(filename, "particles_%03zu.bin", nt);
  particles_to_file(filename, p);

  if (r == 0) {
    xdmf_write(nt, dt);
  }

  particles_free(&p);

  MPI_Finalize();

  return EXIT_SUCCESS;
}

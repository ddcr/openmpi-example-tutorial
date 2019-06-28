#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

void xdmf_write(size_t nt, double dt) {
  FILE* particle_file = fopen("particles_000.bin", "rb");
  if (!particle_file) {
    perror("Could not open particle file particles_000.bin for inspection");
    exit(EXIT_FAILURE);
  }

  uint64_t np_;
  size_t nread = fread(&np_, sizeof(uint64_t), 1, particle_file);
  if (nread != 1) {
    fprintf(stderr, "Particle file too short or error reading.\n");
    exit(EXIT_FAILURE);
  }
  size_t np = np_;

  fclose(particle_file);

  FILE* file = fopen("particles.xdmf", "w");
  if (!file) {
    perror("Could not open particle output XDMF file");
    exit(EXIT_FAILURE);
  }

  fprintf(file, "<?xml version=\"1.0\" ?>\n<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n<Xdmf Version=\"2.0\">\n<Domain>\n<Grid Name=\"TimeSeries\" GridType=\"Collection\" CollectionType=\"Temporal\">\n");

  char filename[21];
  for (size_t it = 0; it <= nt; ++it) {
    sprintf(filename, "particles_%03zu.bin", it);
    // Grid
    fprintf(file, "<Grid Name=\"Particles\" GridType=\"Uniform\">");

    // Time
    fprintf(file, "<Time Value=\"%f\" />\n", it * dt);

    // Topology
    fprintf(file, "<Topology TopologyType=\"Polyvertex\" NodesPerElement=\"%zu\" />\n", np);

    // Geometry
    fprintf(file, "<Geometry GeometryType=\"X_Y_Z\">\n");
    for (size_t i = 0; i < 3; ++i) {
      fprintf(file, "<DataItem ItemType=\"Uniform\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Seek=\"%zu\" Dimensions=\"%zu\">\n", sizeof(uint64_t) + i * np * sizeof(double), np);
      fprintf(file, "%s\n", filename);
      fprintf(file, "</DataItem>\n");
    }
    fprintf(file, "</Geometry>\n");

    // Charge
    fprintf(file, "<Attribute Name=\"Charge\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(file, "<DataItem ItemType=\"Uniform\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Seek=\"%zu\" Dimensions=\"%zu\">%s</DataItem>\n", sizeof(uint64_t) + 3 * np * sizeof(double), np, filename);
    fprintf(file, "</Attribute>");

    // Velocity
    fprintf(file, "<Attribute Name=\"Velocity\" AttributeType=\"Vector\" Center=\"Node\">\n");
    fprintf(file, "<DataItem ItemType=\"Function\" DataType=\"Float\" Precision=\"8\" Dimensions=\"%zu 3\" Function=\"JOIN($0, $1, $2)\">\n", np);
    for (size_t i = 4; i < 7; ++i) {
      fprintf(file, "<DataItem ItemType=\"Uniform\" Format=\"Binary\" NumberType=\"Float\" Precision=\"8\" Seek=\"%zu\" Dimensions=\"%zu\">\n", sizeof(uint64_t) + i * np * sizeof(double), np);
      fprintf(file, "%s\n", filename);
      fprintf(file, "</DataItem>\n");
    }
    fprintf(file, "</DataItem>");
    fprintf(file, "</Attribute>");

    fprintf(file, "</Grid>\n");
  }

  fprintf(file, "</Grid>\n</Domain>\n</Xdmf>\n");

  if (fclose(file)) {
    perror("Could not close particle output XDMF file");
    exit(EXIT_FAILURE);
  }
}

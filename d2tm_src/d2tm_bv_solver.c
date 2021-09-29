#include "bem2_emf_tm.h"

int main(int argc,char **argv)
{
  DOMD md;
  
  read_data(argc,argv,&md);
  print_data(&md);
  //print_data_MKSA(&md);
  initialize_domd(&md);
  output_node_particles(argv[4],&md);

  solve_bieq(&md);
  dat_write(argv[4],&md);

  mfree_domd(&md);
  return 0;
}

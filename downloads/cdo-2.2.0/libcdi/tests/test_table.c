#include <assert.h>
#include <stdlib.h>

#include "cdi.h"

int
main()
{
  tableInqEntry(-1, -1, -1, NULL, NULL, NULL);
  // assert(tableInqEntry(-1, -1, -1, NULL, NULL, NULL) != 0);
  return EXIT_SUCCESS;
}

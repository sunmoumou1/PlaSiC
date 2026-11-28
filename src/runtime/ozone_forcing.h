#ifndef PLASIC_OZONE_FORCING_H
#define PLASIC_OZONE_FORCING_H

#include <stdint.h>

int runtime_ozone_forcing_open(const char *path);
/* Check coverage and values for steps starting at plasic_nstep. */
int runtime_ozone_forcing_validate(int32_t steps);
int runtime_ozone_forcing_apply(void);
void runtime_ozone_forcing_close(void);
const char *runtime_ozone_forcing_source(void);

#endif

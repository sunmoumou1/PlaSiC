#ifndef PLASIC_OZONE_FORCING_H
#define PLASIC_OZONE_FORCING_H

int runtime_ozone_forcing_open(const char *path);
int runtime_ozone_forcing_apply(void);
void runtime_ozone_forcing_close(void);
const char *runtime_ozone_forcing_source(void);

#endif

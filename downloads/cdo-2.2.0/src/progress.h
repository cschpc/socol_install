#ifndef PROGRESS_H
#define PROGRESS_H

namespace progress
{

extern bool stdoutIsTerminal;
extern bool silentMode;

void init();
void init(const char *p_context);
void update(double offset, double refval, double curval);
void set_context_function(const char *(*func)(void) );

}  // namespace progress

#endif

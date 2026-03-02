#ifndef TIMER_H
#define TIMER_H

extern int timer_read, timer_write;

int timer_new(const char *text);
double timer_val(int it);
void timer_report(void);
void timer_start(int it);
void timer_stop(int it);

#endif

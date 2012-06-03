#ifndef TIMER_H
#define TIMER_H

/* Start timer */
void timer_start (void);
/* CPU time from last call with MESSAGE printed */
void print_time (const char *message);
/* CPU and Astronomic time from last call with MESSAGE printed */
void print_full_time (const char *message);

/* Same as timer_start */
void TimerStart (void);
/* CPU time from last call with MESSAGE printed
   returns CPU time */
long PrintTime (const char *message);
/* CPU and Astronomic time from last call with MESSAGE printed
   returns CPU time, Astronomic time stored in pTotalTime */
long int PrintTimeT (const char *message, long int * pTotalTime);
/* Milliseconds from start of the process */
long TimerGet (void);

/* Print current time in string BUFFER */
void sprint_time (char *buffer);

/* CPU time for pthread in nanosec */
double get_time_pthread(void);

#endif // TIMER_H

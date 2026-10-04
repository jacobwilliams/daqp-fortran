/*
 * Helpers for the comparison with upstream DAQP: the parts of the C split
 * interface that need the size or the fields of DAQPWorkspace, which the
 * Fortran side does not mirror.
 */
#include "api.h"
#include <stdlib.h>

/* A new workspace that uses the given settings (owned by the caller) */
void* cmp_ws_new(DAQPSettings* settings){
    DAQPWorkspace* work = calloc(1, sizeof(DAQPWorkspace));
    work->settings = settings;
    return work;
}

/* Free a workspace from cmp_ws_new (and set up by setup_daqp) */
void cmp_ws_free(void* p){
    DAQPWorkspace* work = (DAQPWorkspace*)p;
    if(work == NULL) return;
    work->settings = NULL; /* owned by the caller */
    free_daqp_workspace(work);
    free_daqp_ldp(work);
    free(work);
}

int cmp_ws_n_active(void* p){
    return ((DAQPWorkspace*)p)->n_active;
}

/* The working set (0-based indices) and whether each is at its lower bound */
void cmp_ws_working_set(void* p, int* ws, int* lower){
    DAQPWorkspace* work = (DAQPWorkspace*)p;
    int i;
    for(i = 0; i < work->n_active; i++){
        ws[i] = work->WS[i];
        lower[i] = (work->sense[work->WS[i]] & DAQP_LOWER) != 0;
    }
}

/* Set the working set as daqp_set_working_set of the Fortran port does */
int cmp_ws_set_working_set(void* p, int na, const int* ids, const int* lower){
    DAQPWorkspace* work = (DAQPWorkspace*)p;
    int i, id;
    for(i = 0; i < work->m; i++)
        if(!(work->sense[i] & DAQP_IMMUTABLE)) work->sense[i] &= ~DAQP_ACTIVE;
    for(i = 0; i < na; i++){
        id = ids[i];
        if(work->sense[id] & DAQP_IMMUTABLE) continue;
        work->sense[id] |= DAQP_ACTIVE;
        if(lower[i]) work->sense[id] |= DAQP_LOWER;
        else work->sense[id] &= ~DAQP_LOWER;
    }
    reset_daqp_workspace(work);
    return daqp_activate_constraints(work);
}

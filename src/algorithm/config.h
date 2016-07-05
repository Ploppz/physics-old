/* Only for Resolution alg. Config. */


/*** Physical ***/
float RESTITUTION = 1;


/*** Iterating ***/
int MAX_ITERATIONS = 2;

/** Depth resolution **/
int MAX_DEPTH_ITER = 10;
float INSIGNIFICANT_DEPTH = 0.1f;

/** Time resolution **/
int MAX_TIME_ITER = 6;
float MIN_REWIND_TIME = 0.01f;
float MIN_REL_VEL = 1;
/* In numerically approximating the point in time of contact: */
int MAX_REWIND_ITER = 100;
float ERROR_THRESHOLD = 0.1f; // Acceptable distance that counts as contact

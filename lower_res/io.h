//--- Declare funtions

int load_snapshot(char *fname, int files);
int allocate_memory(void);
int savepositions_ioformat1(char *fname);
size_t my_fwrite(void *ptr, size_t size, size_t nmemb, FILE *stream);
size_t my_fread(void *ptr, size_t size, size_t nmemb, FILE *stream);

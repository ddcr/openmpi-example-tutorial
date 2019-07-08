/*
* @Author: Domingos Rodrigues <ddcr@lcc.ufmg.br>
* @Date:   2019-06-28 11:32:45
* @Last Modified by:   Domingos Rodrigues
* @Last Modified time: 2019-07-07 23:54:05
*
* Routines imported from SLURM sources
*  <https://www.schedmd.com>
*/

#include <sys/types.h>

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <assert.h>
#include <limits.h>

#include "xstring.h"

#define xrealloc(__p, __sz) new_xrealloc((void **)&(__p), __sz)
void *xmalloc(size_t);
void *new_xrealloc(void **, size_t);
int xsize(void *);

#define XMALLOC_MAGIC 0x42
#define XFGETS_CHUNKSIZE 64

static void makespace(char **str, int needed)
{
	int used;

    if (*str == NULL)
    	*str = xmalloc(needed + 1);
    else {
    	used = strlen(*str) + 1;
        while (used + needed > xsize(*str)) {
            int newsize = xsize(*str) + XFGETS_CHUNKSIZE;
            int actualsize;

            xrealloc(*str, newsize);
            actualsize = xsize(*str);

            assert(actualsize == newsize);
        }
    }
}

void *xmalloc(size_t size)
{
    void *new;
    int *p;

    assert(size >= 0 && size <= INT_MAX);
    p = (int *)malloc(size + 2*sizeof(int));
    if (!p)
        exit(1);
    p[0] = XMALLOC_MAGIC;   /* add "secret" magic cookie */
    p[1] = (int)size;       /* store size in buffer */

    new = &p[2];
    memset(new, 0, size);
    return new;
}

void * new_xrealloc(void **item, size_t newsize)
{
    int *p = NULL;

    /* xmalloc_assert(*item != NULL, file, line, func); */
    assert(newsize >= 0 && (int)newsize <= INT_MAX);

    if (*item != NULL) {
        int old_size;
        p = (int *)*item - 2;

        /* magic cookie still there? */
        assert(p[0] == XMALLOC_MAGIC);
        old_size = p[1];

        p = (int *)realloc(p, newsize + 2*sizeof(int));

        if (p == NULL)
            goto error;

        if (old_size < newsize) {
            char *p_new = (char *)(&p[2]) + old_size;
            memset(p_new, 0, (int)(newsize-old_size));
        }
        assert(p[0] == XMALLOC_MAGIC);

    } else {
        /* Initalize new memory */
        p = (int *)malloc(newsize + 2*sizeof(int));
        if (p == NULL)
            goto error;

        memset(&p[2], 0, newsize);
        p[0] = XMALLOC_MAGIC;
    }

    p[1] = (int)newsize;
    *item = &p[2];
    return *item;
    error:
    abort();
}


int xsize(void *item)
{
    int *p = (int *)item - 2;
    assert(item != NULL);
    assert(p[0] == XMALLOC_MAGIC);
    return p[1];
}




char *xstrdup(const char *str)
{
	size_t siz, rsiz;
    char   *result;

    if (str == NULL) {
        return NULL;
    }
    siz = strlen(str) + 1;
    result = (char *)xmalloc(siz);

    rsiz = strlcpy(result, str, siz);

    assert(rsiz == siz-1);

    return result;
}


char *xstrdup_printf(const char *fmt, ...)
{
    /* Start out with a size of 100 bytes. */
    int n, size = 100;
    char *p = NULL;
    va_list ap;

    if((p = xmalloc(size)) == NULL)
        return NULL;
    while(1) {
        /* Try to print in the allocated space. */
        va_start(ap, fmt);
        n = vsnprintf(p, size, fmt, ap);
        va_end (ap);
        /* If that worked, return the string. */
        if (n > -1 && n < size)
            return p;
        /* Else try again with more space. */
        if (n > -1)               /* glibc 2.1 */
            size = n + 1;           /* precisely what is needed */
        else                      /* glibc 2.0 */
            size *= 2;              /* twice the old size */
        if ((p = xrealloc(p, size)) == NULL)
            return NULL;
    }
}


size_t strlcpy(dst, src, siz)
	char *dst;
    const char *src;
    size_t siz;
{
    register char *d = dst;
    register const char *s = src;
    register size_t n = siz;

    /* Copy as many bytes as will fit */
    if (n != 0 && --n != 0) {
        do {
            if ((*d++ = *s++) == 0)
                break;
        } while (--n != 0);
    }

    /* Not enough room in dst, add NUL and traverse rest of src */
    if (n == 0) {
        if (siz != 0)
            *d = '\0';              /* NUL-terminate dst */
        while (*s++)
            ;
    }

    return(s - src - 1);    /* count does not include NUL */
}


void _xstrcat(char **str1, const char *str2)
{
    if (str2 == NULL)
        str2 = "(null)";

    makespace(str1, strlen(str2));
    strcat(*str1, str2);
}


#ifdef TESTBED
#include <sys/param.h>
#include <unistd.h>


int main(int argc, char const *argv[])
{
    char *newroot;
    char *buf = NULL;
    char *filename = NULL;
    char *env_val = NULL;
    char *tmp;

    newroot = xstrdup("root");
    printf("newroot = %s\n", newroot);

    xstrcat(buf, "CPU_Bind");
    xstrcat(buf, ",");
    xstrcat(buf, " Memory_Bind");
    printf("buf= %s\n", buf);

    if ((env_val = getenv("TMPDIR"))) {
        printf("TMPDIR has been set to: %s\n", env_val);
        filename = xstrdup_printf("%s/testfile%03d.dat", env_val, 2);
    } else {
        printf("TMPDIR is not set!\n");
        filename = xstrdup_printf("testfile%03d.dat", 2);
    }
    printf("%s\n", filename);

    tmp = malloc(PATH_MAX);
    if ((tmp = getcwd(tmp, PATH_MAX)))
        printf("getcwd set to :=>%s<=\n", tmp);
    else
        printf("getcwd failed\n");

    // if ((tmp = get_current_dir_name()))
    //     printf("getcwd set to :=>%s<=\n", tmp);
    // else
    //     printf("getcwd failed\n");


    return 0;
}
#endif

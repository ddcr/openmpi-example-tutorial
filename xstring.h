/*
* @Author: Domingos Rodrigues <ddcr@lcc.ufmg.br>
* @Date:   2019-06-27 07:42:54
* @Last Modified by:   ddcr
* @Last Modified time: 2019-06-28 13:26:46
*/
#ifndef _XSTRING_H
#define _XSTRING_H 1

char *xstrdup(const char *);
char *xstrdup_printf(const char *, ...);
size_t strlcpy(char *dst, const char *src, size_t siz);

#define xstrcat(__p, __q) _xstrcat(&(__p), __q)
void _xstrcat(char **str1, const char *str2);

#endif
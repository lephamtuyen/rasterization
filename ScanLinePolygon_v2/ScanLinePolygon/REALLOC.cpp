#include <stdlib.h>

void *REALLOC(void *p, int n, int size)
{
	if (n > 0)
	{
		if (p)
			return (realloc(p, n * size));
		else
			return (calloc(n, size));
	}
	else
	{
		if (p) free(p);
		return (NULL);
	}
}

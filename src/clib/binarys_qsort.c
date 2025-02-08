#include <stdint.h>
#include "binarys_qsort.h"

int binsearch_range(uint64_t key, uint32_t *v, int64_t n,  int64_t *range, int8_t k_off)
{
    int64_t l=0, r=n-1, m;
    uint32_t tmp = 0;
    range[0] = range[1] = -1;

    //printf("k_off = %d\n", k_off);

    if (k_off == 0)
    {
        while(l <= r)
        {
            m = (l + r)/2;
            if (key < v[m])
            {
                r = m - 1;
            }
            else if (key > v[m])
            {
                l = m + 1;
            }
            else
            {
                range[0] = range[1] = m;
                return 1;
            }
        }
    }
    else
    {
        while (l <= r)
        {
            m = (l+r)/2;
            tmp = v[m] >> k_off;
            if (tmp == key)
            {
                range[0] = range[1] = m;
                
                //run low bound
                int64_t sl=l, sr=m-1, sm;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[0] = sm;
                        sr = sm-1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }

                //run upper bound
                sl = m+1; sr = r;
                while (sl <= sr)
                {
                    sm = (sl+sr)/2;
                    tmp = v[sm] >> k_off;
                    if (tmp == key)
                    {
                        range[1] = sm;
                        sl = sm+1;
                    }
                    else if (tmp > key) sr = sm - 1;
                    else    sl = sm + 1;
                }
                return 1;
            }
            else if (tmp > key) r = m - 1;
            else l = m + 1;
        }
    }

    return -1;
}

#ifdef UNPIPATH_OFF_K20
int multi_binsearch_offset64(uint32_t x, uint32_t v[], int64_t n, uint64_t offset, int64_t seed_binary[], int8_t k_r)
#else
int multi_binsearch_offset(uint32_t x, uint32_t v[], int64_t n, uint32_t offset, int64_t seed_binary[], int8_t k_r)
#endif
{
    int64_t low, high, mid;
    uint32_t temp = 0;
#ifdef UNPIPATH_OFF_K20
    uint64_t current_offset = 0;
#else
    uint32_t current_offset = 0;
#endif

    int8_t k_2r = k_r << 1;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        temp = v[mid + offset] >> k_2r;
        if(x < temp)
        {
            high = mid - 1;
        }
        else if(x > temp)
        {
            low = mid + 1;
        }
        else  /*found match*/
        {   
            //找到当前match的kmer，然后向两侧search，找到match的上下界
            if (mid == 0)
            {
                seed_binary[0] = offset;

                current_offset = mid + offset + 1;
                while((current_offset < offset + n) &&  (x == (v[current_offset] >> k_2r)))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }
            else if (mid == (n - 1))
            {
                seed_binary[1] = n + offset- 1;

                current_offset = mid + offset - 1;
                while((current_offset >= offset) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset--;
                }
                seed_binary[0] = current_offset + 1;
            }
            else
            {
                //up
                current_offset = mid + offset - 1;
                while((current_offset >= offset) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset--;
                }
                seed_binary[0] = current_offset + 1;    
                
                //down
                current_offset = mid + offset + 1;
                while((current_offset < offset + n) && (x == (v[current_offset] >> k_2r)))
                {
                    current_offset++;
                }
                seed_binary[1] = current_offset - 1;
            }

            return 1;
        }
    }
    return -1;
}

int64_t binsearch_interval_unipath64(uint64_t x, uint64_t v[], uint64_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}

int64_t binsearch_interval_unipath(uint32_t x, uint32_t v[], uint32_t n)
{
    int64_t low, high, mid;

    low = 0;
    high = n - 1;

    while ( low <= high )
    {
        mid = ((int64_t )(low + high)) >> 1;
        if(x < v[mid])
        {
            high = mid - 1;
        }
        else if(x > v[mid])
        {
            low = mid + 1;
        }
        else  /*found match*/
        {
            return mid;
        }
    }

    return high;
}




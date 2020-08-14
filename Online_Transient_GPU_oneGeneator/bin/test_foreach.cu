#include <thrust/for_each.h>
#include <thrust/iterator/counting_iterator.h>

struct sum_Functor {
    int *sum;
    sum_Functor(int *s){sum = s;}
    __host__ __device__
    void operator()(int i)
    {
        *sum+=i;
        printf("In functor: i %d sum %d\n",i,*sum);
    }

};

int main(){

    thrust::counting_iterator<int> first(0);
    thrust::counting_iterator<int> last = first+10;
    int sum = 0;
    sum_Functor sf(&sum);
    printf("After constructor: value is %d\n", *(sf.sum));
    for(int i=0;i<5;i++){
        sf(i);
    }

    printf("Initiating for_each call - current value %d\n", (*(sf.sum)));
    thrust::for_each(first,last,sf);

    cudaDeviceSynchronize();
    printf("After for_each: value is %d\n",*(sf.sum));
}

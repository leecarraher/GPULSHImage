#include<stdio.h>

int MIN_CLOSE_DIST = 128*10;
int main(int l,char** args)
{
    printf("hello world\n");

}

//This function finds the nearest key greater than grtThan
int findMatch(int lenB,int r,int* keypointsA, int* keypointsB)
{
    int d,mind=100000000 , mind2 = 10000000;
    int minidx1 = 0;
    int j;

    for(j=0;j<lenB/128;j++)
    {
        
        d = DistSquared(r,j*128,keypointsA,keypointsB);
        if(d<mind)
        {
          mind2=mind;
          mind=d;
          minidx1 = j;
        }else if(d<mind2)mind2=d;
    }
    if ( 10*10 * mind < 6*6* mind2)return 1;
    /*second round
    for(j=0;j<lenB/128;j++)
    {
        if (keypointsB[j*128]!=-1)
            d = DistSquared(r,j*128,keypointsA,keypointsB);
        if(d<mind2 && d>mind)
        {
          mind2=d;

        }
    }
    //.6 metric from sift match demo code
    if(mind < (mind2*.6)){
        //invalidate the nearest key
        keypointsB[minidx1*128]=-1;
        return 1;
    }*/
    return 0;


}

int mock(){
    return 0;
}

int numMatches(int lenA,int lenB, int* keypointsA, int* keypointsB)
{
    int matches=0;
    int k;
    for (k=0;k<lenA/128;k++)
    {
        matches=matches+findMatch(lenB,k*128,keypointsA, keypointsB);
    }
    return matches;

}


/* Return squared distance between two keypoint descriptors.
*/
int DistSquared(int s1,int s2,int* keypointsA, int* keypointsB)
{
    int i, dif, distsq = 0;
    for (i = 0; i < 128; i++) {
      dif = keypointsA[s1+i] - keypointsB[s2+i];
      distsq += dif * dif;
    }
    return distsq;
}
